module physicalParticleState_class

  use errors_mod,                 only : fatalError
  use numPrecision
  use transportObjectState_class, only : buildTransportObjectStatePayload, display_super => display, init_super => init, &
                                         kill_super => kill, preparePayload_super => preparePayload, transportObjectState

  implicit none
  private

  ! Public procedures.
  public :: castPhysicalParticleStatePtr, display, init, kill, preparePayload

  !!
  !!
  !!
  type, public, extends(buildTransportObjectStatePayload) :: buildPhysicalParticleStatePayload
    integer(shortInt) :: broodId = 0
    real(defReal)     :: mass = ZERO
  end type buildPhysicalParticleStatePayload

  !!
  !!
  !!
  type, public, extends(transportObjectState) :: physicalParticleState
    private
    integer(shortInt) :: broodId = 0, nCollisions = 0
    real(defReal)     :: mass = ZERO
  contains
    procedure :: display
    procedure :: getBroodId
    procedure :: getCollisionsNumber
    procedure :: getMass
    procedure :: init
    procedure :: kill
    procedure :: preparePayload
    procedure :: setBroodId
    procedure :: setCollisionsNumber
    procedure :: setMass
  end type physicalParticleState

  !!
  !!
  !!
  type, public :: physicalParticleStateBox
    class(physicalParticleState), pointer :: ptr => null()
  end type physicalParticleStateBox

contains
  !!
  !!
  !!
  function castPhysicalParticleStatePtr(source, fatal) result(ptr)
    class(transportObjectState), intent(in) :: source
    logical(defBool), intent(in), optional  :: fatal
    logical(defBool)                        :: throwError
    class(physicalParticleState), pointer   :: ptr
    character(*), parameter                 :: here = 'castPhysicalParticleStatePtr (physicalParticleState_class.f90)'

    select type(temp => source)
      class is(physicalParticleState)
        ptr => temp

      class default
        ptr => null()

    end select

    ! Throw error if requested.
    throwError = .false.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) &
    call fatalError(here, "Transport object state is not of class 'physicalParticleState'.")

  end function castPhysicalParticleStatePtr

  !!
  !!
  !!
  subroutine display(self)
    class(physicalParticleState), intent(in) :: self

    ! Superclass.
    call display_super(self)

    ! Local.
    print *, 'Brood id: ', self % broodId
    print *, 'Number of collisions: ', self % nCollisions

  end subroutine display

  !!
  !!
  !!
  function getBroodId(self) result(broodId)
    class(physicalParticleState), intent(in) :: self
    integer(shortInt)                        :: broodId

    broodId = self % broodId

  end function getBroodId

  !!
  !!
  !!
  elemental function getCollisionsNumber(self) result(nCollisions)
    class(physicalParticleState), intent(in) :: self
    integer(shortInt)                        :: nCollisions

    nCollisions = self % nCollisions

  end function getCollisionsNumber

  !!
  !!
  !!
  elemental function getMass(self) result(mass)
    class(physicalParticleState), intent(in) :: self
    real(defReal)                            :: mass

    mass = self % mass

  end function getMass

  !!
  !!
  !!
  subroutine init(self, payload)
    class(physicalParticleState), intent(inout)         :: self
    class(buildTransportObjectStatePayload), intent(in) :: payload
    class(buildPhysicalParticleStatePayload), pointer   :: payloadPtr
    character(*), parameter                             :: HERE = 'init (physicalParticleState_class.f90)'

    ! Initialise superclass.
    call init_super(self, payload)

    ! Downcast payload to correct class.
    select type(ptr => payload)
      class is(buildPhysicalParticleStatePayload)
        payloadPtr => ptr

      class default
        call fatalError(HERE, 'Invalid payload type.')

    end select

    ! Local.
    self % broodId = payloadPtr % broodId
    self % mass = payloadPtr % mass

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(physicalParticleState), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % broodId = 0
    self % nCollisions = 0
    self % mass = ZERO

  end subroutine kill

  !!
  !!
  !!
  subroutine preparePayload(self, payload)
    class(physicalParticleState), intent(in)               :: self
    class(buildTransportObjectStatePayload), intent(inout) :: payload
    class(buildPhysicalParticleStatePayload), pointer      :: payloadPtr
    character(*), parameter                                :: HERE = 'preparePayload (physicalParticleState_class.f90)'

    ! Superclass.
    call preparePayload_super(self, payload)

    ! Downcast payload to correct class.
    select type(ptr => payload)
      class is(buildPhysicalParticleStatePayload)
        payloadPtr => ptr

      class default
        call fatalError(HERE, 'Invalid payload type.')

    end select

    ! Local.
    payloadPtr % broodId = self % broodId
    payloadPtr % mass = self % mass

  end subroutine preparePayload

  !!
  !!
  !!
  elemental subroutine setBroodId(self, broodId)
    class(physicalParticleState), intent(inout) :: self
    integer(shortInt), intent(in)               :: broodId

    self % broodId = broodId

  end subroutine setBroodId

  !!
  !!
  !!
  elemental subroutine setCollisionsNumber(self, nCollisions)
    class(physicalParticleState), intent(inout) :: self
    integer(shortInt), intent(in)               :: nCollisions

    self % nCollisions = nCollisions

  end subroutine setCollisionsNumber

  !!
  !!
  !!
  elemental subroutine setMass(self, mass)
    class(physicalParticleState), intent(inout) :: self
    real(defReal), intent(in)                   :: mass

    self % mass = mass

  end subroutine setMass

end module physicalParticleState_class