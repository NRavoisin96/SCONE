module distributedSource_inter

  use dictionary_class,           only : dictionary
  use errors_mod,                 only : fatalError
  use genericProcedures,          only : numToChar, rotateVector
  use geometry_inter,             only : geometry
  use neutronMaterial_inter,      only : neutronMaterial, neutronMaterial_CptrCast
  use nuclearDatabase_inter,      only : nuclearDatabase
  use nuclearDataReg_mod,         only : get
  use numPrecision
  use RNG_class,                  only : RNG
  use source_inter,               only : init_super => init, kill_super => kill, source
  use transportObjectState_class, only : buildTransportObjectStatePayload, transportObjectState

  implicit none
  private

  ! Public procedures.
  public :: init, kill

  !!
  !!
  !!
  type, public, abstract, extends(source) :: distributedSource
    private
    integer(shortInt)              :: nMaxIterations = 0
    real(defReal), dimension(3, 2) :: bounds = ZERO
  contains
    procedure(finalisePayload), deferred               :: finalisePayload
    procedure(getDefaultMaxIterationsNumber), deferred :: getDefaultMaxIterationsNumber
    procedure                                          :: init
    procedure(isMaterialInvalid), deferred             :: isMaterialInvalid
    procedure                                          :: kill
    procedure(printInfiniteLoopError), deferred        :: printInfiniteLoopError
    procedure                                          :: rejectMaterial
    procedure                                          :: sampleState
  end type distributedSource

  abstract interface
    !!
    !!
    !!
    subroutine finalisePayload(self, mat, database, temperature, rand, payload, mu, phi)
      import :: buildTransportObjectStatePayload, defReal, distributedSource, neutronMaterial, nuclearDatabase, RNG
      class(distributedSource), intent(in)                   :: self
      class(neutronMaterial), intent(in)                     :: mat
      class(nuclearDatabase), intent(in)                     :: database
      real(defReal), intent(in)                              :: temperature
      type(RNG), intent(inout)                               :: rand
      class(buildTransportObjectStatePayload), intent(inout) :: payload
      real(defReal), intent(out)                             :: mu, phi
    end subroutine finalisePayload

    !!
    !!
    !!
    elemental function getDefaultMaxIterationsNumber(self) result(nDefaultMaxIterations)
      import                               :: distributedSource, shortInt
      class(distributedSource), intent(in) :: self
      integer(shortInt)                    :: nDefaultMaxIterations
    end function getDefaultMaxIterationsNumber

    !!
    !!
    !!
    elemental function isMaterialInvalid(self, materialIdx) result(isInvalid)
      import                               :: defBool, distributedSource, shortInt
      class(distributedSource), intent(in) :: self
      integer(shortInt), intent(in)        :: materialIdx
      logical(defBool)                     :: isInvalid
    end function isMaterialInvalid

    !!
    !!
    !!
    subroutine printInfiniteLoopError(self, nMaxIterations)
      import                               :: distributedSource, shortInt
      class(distributedSource), intent(in) :: self
      integer(shortInt), intent(in)        :: nMaxIterations
    end subroutine printInfiniteLoopError

  end interface

contains
  !!
  !!
  !!
  subroutine init(self, dict, geom)
    class(distributedSource), intent(inout)  :: self
    class(dictionary), intent(in)            :: dict
    class(geometry), pointer, intent(in)     :: geom
    integer(shortInt)                        :: nBounds
    real(defReal), dimension(:), allocatable :: bounds
    real(defReal), dimension(6)              :: geometryBounds
    character(*), parameter                  :: here = 'init (distributedSource_inter.f90)'

    ! Initialise superclass.
    call init_super(self, dict, geom)

    ! Load the maximum number of iterations allowed and check that it is valid.
    call dict % getOrDefault(self % nMaxIterations, 'maxIterations', self % getDefaultMaxIterationsNumber())
    if (self % nMaxIterations < 1) &
    call fatalError(here, 'Maximum number of iterations allowed: '//numToChar(self % nMaxIterations)//' is negative.')

    ! Load bounding box if provided, else set as the geometry's default bounding box.
    geometryBounds = geom % bounds()
    self % bounds(:, 1) = geometryBounds(1:3)
    if (dict % isPresent('bottom')) then
      call dict % get(bounds, 'bottom')
      nBounds = size(bounds)
      if (nBounds /= 3) call fatalError(here, 'Source lower bounds must have dimension 3. Has: '//numToChar(nBounds)//'.')
      self % bounds(:, 1) = bounds

    end if

    self % bounds(:, 2) = geometryBounds(4:6)
    if (dict % isPresent('top')) then
      call dict % get(bounds, 'top')
      nBounds = size(bounds)
      if (nBounds /= 3) call fatalError(here, 'Source upper bounds must have dimension 3. Has: '//numToChar(nBounds)//'.')
      self % bounds(:, 2) = bounds

    end if

    ! Check that bounding box is valid.
    if (any(self % bounds(:, 2) < self % bounds(:, 1))) &
    call fatalError(here, 'Source upper bounds must be strictly greater than lower bounds.')

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(distributedSource), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % nMaxIterations = 0
    self % bounds = ZERO

  end subroutine kill

  !!
  !!
  !!
  elemental function rejectMaterial(self, material) result(isRejected)
    class(distributedSource), intent(in) :: self
    class(neutronMaterial), intent(in)   :: material
    logical(defBool)                     :: isRejected

    ! By default do nothing.
    isRejected = .false.

  end function rejectMaterial

  !!
  !!
  !!
  subroutine sampleState(self, rand, state)
    class(distributedSource), intent(inout)               :: self
    class(RNG), intent(inout)                             :: rand
    class(transportObjectState), allocatable, intent(out) :: state
    class(buildTransportObjectStatePayload), allocatable  :: payload
    class(geometry), pointer                              :: geometryPtr
    class(neutronMaterial), pointer                       :: mat
    class(nuclearDatabase), pointer                       :: nuclearDatabasePtr
    integer(shortInt)                                     :: i
    real(defReal)                                         :: mu, phi, temperature
    character(*), parameter                               :: here = 'sampleParticle (distributedSource_inter.f90)'

    ! Sample state then get pointer to appropriate nuclear database.
    call self % allocatePayloadAndState(payload, state)
    nuclearDatabasePtr => get(state)

    ! Get pointer to geometry.
    geometryPtr => self % getGeometryPtr()

    ! Sample particle.
    i = 0
    rejection : do
      ! Protect against infinite loop
      i = i + 1
      if (self % nMaxIterations < i) call self % printInfiniteLoopError(self % nMaxIterations)

      ! Sample initial position.
      call geometryPtr % sampleInitialPosition(self % bounds(:, 1), self % bounds(:, 2), rand, &
                                               payload % materialIdx, payload % uniqueId, payload % rGlobal, temperature)

      ! Check if material needs to be rejected and cycle if so.
      if (self % isMaterialInvalid(payload % materialIdx)) cycle rejection
      mat => neutronMaterial_CptrCast(nuclearDatabasePtr % getMaterial(payload % materialIdx))
      if (.not. associated(mat)) call fatalError(here, "Nuclear data did not return neutron material.")

      ! Cycle if material is rejected.
      if (self % rejectMaterial(mat)) cycle rejection

      ! Finalise state depending on specific class then initialise state.
      call self % finalisePayload(mat, nuclearDatabasePtr, temperature, rand, payload, mu, phi)
      payload % uGlobal = rotateVector([ONE, ZERO, ZERO], mu, phi)
      call state % init(payload)

      ! Exit loop.
      exit rejection

    end do rejection

  end subroutine sampleState

end module distributedSource_inter