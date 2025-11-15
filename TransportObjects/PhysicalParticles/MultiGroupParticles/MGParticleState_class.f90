module MGParticleState_class

  use errors_mod,                  only : fatalError
  use numPrecision
  use physicalParticleState_class, only : buildPhysicalParticleStatePayload, display_super => display, &
                                          init_fromPayload_super => init_fromPayload, kill_super => kill, &
                                          physicalParticleState, preparePayload_super => preparePayload
  use transportObjectState_class,  only : buildTransportObjectStatePayload, transportObjectState

  implicit none
  private

  ! Public procedures.
  public :: castBuildMGParticleStatePayloadPtr, castMGParticleStatePtr

  !!
  !!
  !!
  type, public, extends(buildPhysicalParticleStatePayload) :: buildMGParticleStatePayload
    integer(shortInt) :: energyGroup = 0
  end type buildMGParticleStatePayload

  !!
  !!
  !!
  type, public, extends(physicalParticleState) :: MGParticleState
    private
    integer(shortInt) :: energyGroup = 0
  contains
    procedure :: display
    procedure :: getEnergyGroup
    procedure :: init_fromPayload
    procedure :: kill
    procedure :: setEnergyGroup
  end type MGParticleState

contains
  !!
  !!
  !!
  function castBuildMGParticleStatePayloadPtr(source, fatal) result(ptr)
    class(buildTransportObjectStatePayload), intent(in) :: source
    logical(defBool), intent(in), optional              :: fatal
    class(buildMGParticleStatePayload), pointer         :: ptr
    logical(defBool)                                    :: throwError
    character(*), parameter :: HERE = 'castBuildMGParticleStatePayloadPtr (MGParticleState_class.f90)'

    select type(temp => source)
      type is(buildMGParticleStatePayload)
        ptr => temp

      class default
        ptr => null()

    end select

    ! Throw error if requested.
    throwError = .true.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) &
    call fatalError(HERE, "Payload is not of type 'buildMGParticleStatePayload'.")

  end function castBuildMGParticleStatePayloadPtr

  !!
  !!
  !!
  function castMGParticleStatePtr(source, fatal) result(ptr)
    class(transportObjectState), intent(in) :: source
    logical(defBool), intent(in), optional  :: fatal
    class(MGParticleState), pointer         :: ptr
    logical(defBool)                        :: throwError
    character(*), parameter                 :: HERE = 'castMGParticleStatePtr (MGParticleState_class.f90)'

    select type(temp => source)
      class is(MGParticleState)
        ptr => temp

      class default
        ptr => null()

    end select

    ! Throw error if requested.
    throwError = .false.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) &
    call fatalError(HERE, "Transport object state is not of type 'MGParticleState'.")

  end function castMGParticleStatePtr

  !!
  !!
  !!
  subroutine display(self)
    class(MGParticleState), intent(in) :: self

    ! Superclass.
    call display_super(self)

    ! Local.
    print *, 'Energy group: ', self % energyGroup

  end subroutine display

  !!
  !!
  !!
  elemental function getEnergyGroup(self) result(energyGroup)
    class(MGParticleState), intent(in) :: self
    integer(shortInt)                  :: energyGroup

    energyGroup = self % energyGroup

  end function getEnergyGroup

  !!
  !!
  !!
  subroutine init_fromPayload(self, payload)
    class(MGParticleState), intent(inout)               :: self
    class(buildTransportObjectStatePayload), intent(in) :: payload
    type(buildMGParticleStatePayload), pointer          :: payloadPtr

    ! Initialise superclass.
    call init_fromPayload_super(self, payload)

    ! Downcast payload to correct type.
    payloadPtr => castBuildMGParticleStatePayloadPtr(payload)

    ! Set energy group.
    self % energyGroup = payloadPtr % energyGroup

  end subroutine init_fromPayload

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(MGParticleState), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % energyGroup = 0

  end subroutine kill

  !!
  !!
  !!
  subroutine preparePayload(self, payload)
    class(MGParticleState), intent(in)                     :: self
    class(buildTransportObjectStatePayload), intent(inout) :: payload
    type(buildMGParticleStatePayload), pointer             :: payloadPtr

    ! Superclass.
    call preparePayload_super(self, payload)

    ! Downcast payload to correct type.
    payloadPtr => castBuildMGParticleStatePayloadPtr(payload)

    ! Local.
    payloadPtr % energyGroup = self % energyGroup

  end subroutine preparePayload

  !!
  !!
  !!
  elemental subroutine setEnergyGroup(self, energyGroup)
    class(MGParticleState), intent(inout) :: self
    integer(shortInt), intent(in)         :: energyGroup

    self % energyGroup = energyGroup

  end subroutine setEnergyGroup

end module MGParticleState_class