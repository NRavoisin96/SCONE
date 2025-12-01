module uniFissSitesField_class

  use dictionary_class,           only : dictionary
  use errors_mod,                 only : fatalError
  use field_inter,                only : field
  use genericProcedures,          only : numToChar
  use geometry_inter,             only : geometry
  use neutronMaterial_inter,      only : neutronMaterial, neutronMaterial_CptrCast
  use nuclearDatabase_inter,      only : nuclearDatabase
  use nuclearDataReg_mod,         only : ndReg_getNeutronCE => getNeutronCE, ndReg_getNeutronMG => getNeutronMG
  use numPrecision
  use RNG_class,                  only : RNG
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap
  use transportObject_inter,      only : transportObject
  use transportObjectState_class, only : transportObjectState
  use universalVariables,         only : OUTSIDE_MAT, VOID_MAT, P_NEUTRON_CE
  use vectorField_inter,          only : vectorField

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: uniFissSitesField_TptrCast

  !!
  !! Uniform Fission Sites Field
  !!
  !! Returns a 3D vector with:
  !! - fraction of fissionable material volume occupied by the required cell
  !! - percentage of fission sites in the required cell
  !! - unused entry, filled with a contant. It is there to match the interface
  !!
  !! Sample Dictionary Input:
  !!   uniformFissionSites { type uniFissSitesField;
  !!                         #uniformVolMap 0;#         optional
  !!                         #popVolumes 1.0e7;#        optional
  !!                         map { <map definition> } }
  !!
  !! NOTE: If uniformVolMap is set to 0 (False), the map bins may contain different
  !! volumes of fissile material. If that's the case, the volume in each bin has to be
  !! estimated with a Monte Carlo calculation: we use 'popVolumes' test points
  !! distributed uniformly over the problem geometry to estimate the volume fractions.
  !! This may be associated with a significant stochastic error if the fissile
  !! volume in any cell (or in the geometry in general) is small
  !!
  !! Public Members:
  !!   map ->  map that lays over the geometry. It should be a spatial map, an
  !!           energy map wouldn't make much sense!
  !!   N   ->  total number of map bins
  !!   pop ->  particle population used for the volume estimation
  !!   uniformVolMap  -> flag to indicate whether the map has bins with uniform volumes
  !!   volFraction    -> array with the volume fraction of each bin
  !!   sourceFraction -> array with the percentage of fission sites in each bin
  !!   buildSource    -> array used to 'tally' fission sites
  !!
  !! Interface:
  !!   vectorField interface
  !!
  type, public, extends(vectorField) :: uniFissSitesField
    private
    class(tallyMap), allocatable             :: map
    integer(shortInt)                        :: N = 0, pop = 0
    logical(defBool)                         :: uniformVolMap = .false.
    real(defReal), dimension(:), allocatable :: buildSource, sourceFraction, volFraction
  contains
    ! Superclass interface
    procedure :: init
    procedure :: kill
    procedure :: estimateVol
    procedure :: at
    procedure :: storeFS
    procedure :: updateMap
  end type uniFissSitesField

contains

  !!
  !! Initialise from dictionary
  !!
  !! See field_inter for details
  !!
  subroutine init(self, dict)
    class(uniFissSitesField), intent(inout) :: self
    class(dictionary), intent(in)           :: dict
    integer(shortInt), parameter            :: ALL = 0

    ! Initialise overlay map
    call new_tallyMap(self % map, dict % getDictPtr('map'))
    self % N = self % map % bins(ALL)

    ! Allocate and initialise arrays
    allocate(self % sourceFraction(self % N), self % buildSource(self % N))
    self % sourceFraction = ONE / self % N
    self % buildSource = ZERO

    ! Settings for volume calculation
    call dict % getOrDefault(self % uniformVolMap, 'uniformVolMap', .false.)
    if (.not. self % uniformVolMap) call dict % getOrDefault(self % pop, 'popVolumes', 1000000)

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(uniFissSitesField), intent(inout) :: self

    call self % map % kill()
    deallocate(self % map)
    deallocate(self % sourceFraction)
    deallocate(self % buildSource)

    self % N = 0
    self % uniformVolMap = .false.

  end subroutine kill

  !!
  !! Generate random points to estimate the volume of the elements on the map
  !!
  subroutine estimateVol(self, geom, type, rand)
    class(uniFissSitesField), intent(inout) :: self
    class(geometry), pointer, intent(in)    :: geom
    integer(shortInt), intent(in)           :: type
    type(RNG), intent(inout), optional     :: rand
    real(defReal), dimension(6)             :: bounds
    real(defReal), dimension(3)             :: bottom, top
    real(defReal), dimension(3)             :: r
    type(transportObjectState)              :: state
    integer(shortInt)                       :: i
    integer(shortInt)                       :: j, binIdx, matIdx, uniqueID
    class(nuclearDatabase), pointer         :: nucData
    class(neutronMaterial), pointer         :: mat
    character(*), parameter                 :: Here = 'estimateVol (uniFissSitesField_class.f90)'

    allocate(self % volFraction(self % N))

    ! Check if volume estimation is needed or not
    if (self % uniformVolMap) then
      self % volFraction = ONE / self % N

    else
      self % volFraction = ZERO
      
      ! Get pointer to appropriate nuclear database
      if (type == P_NEUTRON_CE) then
        nucData => ndReg_getNeutronCE()

      else
        nucData => ndReg_getNeutronMG()

      end if
      if (.not. associated(nucData)) call fatalError(Here, 'Failed to retrieve Nuclear Database')

      ! Check that pointer to random number generator is associated.
      if (.not. present(rand)) call fatalError(Here, 'Random number generator was not provided.')

      ! Set bounding region
      bounds = geom % bounds()
      bottom = bounds(1:3)
      top = bounds(4:6)

      print *, "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
      print *, "VOLUME CALCULATION FOR UFS"

      ! Iterate over number of points desired
      !$omp parallel do private(r, state, j, binIdx, matIdx, uniqueID, mat)
      do i = 1, self % pop
        j = 0
        rejection : do
          ! Protect against infinite loop
          j = j + 1
          if (1000 < j) call fatalError(Here, 'Infinite loop in sampling of fission sites. Please check that&
                                              & defined volume contains fissile material.')

          ! Sample initial positionm.
          call geom % sampleInitialPosition(bottom, top, rand, matIdx, uniqueID, r)

          ! Reject if there is no material
          if (any([OUTSIDE_MAT, VOID_MAT] == matIdx)) cycle rejection

          mat => neutronMaterial_CptrCast(nucData % getMaterial(matIdx))
          if (.not. associated(mat)) call fatalError(Here, "Nuclear data did not return neutron material.")

          ! Resample position if material is not fissile
          if (.not. mat % isFissile()) cycle rejection
          call state % setGlobalPosition(r)
          call state % setMaterialIdx(matIdx)

          ! Read map bin index
          binIdx = self % map % map(state)

          ! Return if invalid bin index
          if (binIdx == 0) cycle rejection

          ! Add point to the volume fraction map
          !$omp atomic
          self % volFraction(binIdx) = self % volFraction(binIdx) + 1

          ! Exit the loop
          exit rejection

        end do rejection

      end do
      !$omp end parallel do

      ! Normalise the volume fraction map
      self % volFraction = self % volFraction / sum(self % volFraction)

      print *, "DONE!"
      print *, "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>"

    end if

  end subroutine estimateVol

  !!
  !! Get value of the vector field given the phase-space location of a particle
  !!
  !! See vectorField_inter for details
  !!
  function at(self, object) result(val)
    class(uniFissSitesField), intent(in)  :: self
    class(transportObject), intent(inout) :: object
    real(defReal), dimension(3)           :: val
    integer(shortInt)                     :: binIdx

    ! Read map bin index
    binIdx = self % map % map(object % updateAndGetCurrentStatePtr())

    ! Return if invalid bin index
    val = ONE
    if (binIdx == 0) return
    val = [self % volFraction(binIdx), self % sourceFraction(binIdx), ZERO]

  end function at

  !!
  !! Store the fission sites generated in a vector
  !!
  !! Args:
  !! state [in] -> particle state of the fission site
  !!
  subroutine storeFS(self, state)
    class(uniFissSitesField), intent(inout) :: self
    class(transportObjectState), intent(in) :: state
    integer(shortInt)                       :: idx

    idx = self % map % map(state)
    if (idx == 0) return
    ! Add fission sites where appropriate
    self % buildSource(idx) = self % buildSource(idx) + state % getWeight()

  end subroutine storeFS

  !!
  !! Calculates fission site probability distribution
  !! It accounts for possible ares of the map having zero events
  !!
  subroutine updateMap(self)
    class(uniFissSitesField), intent(inout) :: self
    integer(shortInt) :: i

    ! Eliminate zeros in the distribution
    do i = 1, self % N
      if (self % buildSource(i) == ZERO) self % buildSource(i) = ONE

    end do
    ! Normalise to calculate probability
    self % sourceFraction = self % buildSource / sum(self % buildSource)
    self % buildSource = ZERO

  end subroutine updateMap

  !!
  !! Cast field pointer to uniFissSitesField pointer
  !!
  !! Args:
  !!   source [in] -> source pointer of class field
  !!
  !! Result:
  !!   Null is source is not of uniFissSitesField
  !!   Pointer to source if source is uniFissSitesField type
  !!
  pure function uniFissSitesField_TptrCast(source) result(ptr)
    class(field), pointer, intent(in) :: source
    type(uniFissSitesField), pointer  :: ptr

    select type (source)
      type is (uniFissSitesField)
        ptr => source

      class default
        ptr => null()
    end select

  end function uniFissSitesField_TptrCast


end module uniFissSitesField_class
