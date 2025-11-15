module CEParticleState_class

  use errors_mod,                  only : fatalError
  use numPrecision
  use physicalParticleState_class, only : buildPhysicalParticleStatePayload, display_super => display, &
                                          init_fromPayload_super => init_fromPayload, kill_super => kill, &
                                          physicalParticleState, preparePayload_super => preparePayload
  use transportObjectState_class,  only : buildTransportObjectStatePayload, transportObjectState

  implicit none
  private

  ! Public procedures.
  public :: castBuildCEParticleStatePayloadPtr, castCEParticleStatePtr

  !!
  !!
  !!
  type, public, extends(buildPhysicalParticleStatePayload) :: buildCEParticleStatePayload
    real(defReal) :: energy = ZERO
  end type buildCEParticleStatePayload

  !!
  !!
  !!
  type, public, extends(physicalParticleState) :: CEParticleState
    private
    real(defReal) :: energy = ZERO
  contains
    procedure :: display
    procedure :: getEnergy
    procedure :: init_fromPayload
    procedure :: kill
    procedure :: preparePayload
    procedure :: setEnergy
  end type CEParticleState

contains
  !!
  !!
  !!
  function castBuildCEParticleStatePayloadPtr(source, fatal) result(ptr)
    class(buildTransportObjectStatePayload), intent(in) :: source
    logical(defBool), intent(in), optional              :: fatal
    class(buildCEParticleStatePayload), pointer         :: ptr
    logical(defBool)                                    :: throwError
    character(*), parameter                             :: HERE = 'castBuildCEParticleStatePayloadPtr (CEParticleState_class.f90)'

    select type(temp => source)
      type is(buildCEParticleStatePayload)
        ptr => temp

      class default
        ptr => null()

    end select

    throwError = .true.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) call fatalError(HERE, "Payload is not of type 'buildCEParticleStatePayload'.")

  end function castBuildCEParticleStatePayloadPtr

  !!
  !!
  !!
  function castCEParticleStatePtr(source, fatal) result(ptr)
    class(transportObjectState), intent(in) :: source
    logical(defBool), intent(in), optional  :: fatal
    class(CEParticleState), pointer         :: ptr
    logical(defBool)                        :: throwError
    character(*), parameter                 :: here = 'castCEParticleStatePtr (CEParticleState_class.f90)'

    select type(temp => source)
      class is(CEParticleState)
        ptr => temp

      class default
        ptr => null()

    end select

    throwError = .false.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) call fatalError(here, "Transport object state is not of class 'CEParticleState'.")

  end function castCEParticleStatePtr

  !!
  !!
  !!
  subroutine display(self)
    class(CEParticleState), intent(in) :: self

    ! Superclass.
    call display_super(self)

    ! Local.
    print *, 'Energy: ', self % energy

  end subroutine display

  !!
  !!
  !!
  function getEnergy(self) result(energy)
    class(CEParticleState), intent(in) :: self
    real(defReal)                      :: energy

    energy = self % energy

  end function getEnergy

  !!
  !!
  !!
  subroutine init_fromPayload(self, payload)
    class(CEParticleState), intent(inout)               :: self
    class(buildTransportObjectStatePayload), intent(in) :: payload
    type(buildCEParticleStatePayload), pointer          :: payloadPtr

    ! Initialise superclass.
    call init_fromPayload_super(self, payload)

    ! Downcast payload to correct type.
    payloadPtr => castBuildCEParticleStatePayloadPtr(payload)

    ! Local.
    self % energy = payloadPtr % energy

  end subroutine init_fromPayload

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(CEParticleState), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % energy = ZERO

  end subroutine kill

  !!
  !!
  !!
  subroutine preparePayload(self, payload)
    class(CEParticleState), intent(in)                     :: self
    class(buildTransportObjectStatePayload), intent(inout) :: payload
    type(buildCEParticleStatePayload), pointer             :: payloadPtr

    ! Superclass.
    call preparePayload_super(self, payload)

    ! Downcast payload to correct type.
    payloadPtr => castBuildCEParticleStatePayloadPtr(payload)

    ! Local.
    payloadPtr % energy = self % energy

  end subroutine preparePayload

  !!
  !!
  !!
  elemental subroutine setEnergy(self, energy)
    class(CEParticleState), intent(inout) :: self
    real(defReal), intent(in)             :: energy

    self % energy = energy

  end subroutine setEnergy

end module CEParticleState_class