module CEParticleState_class

  use errors_mod,                  only : fatalError
  use numPrecision
  use physicalParticleState_class, only : buildPhysicalParticleStatePayload, display_super => display, init_super => init, &
                                          kill_super => kill, preparePayload_super => preparePayload, physicalParticleState
  use transportObjectState_class,  only : buildTransportObjectStatePayload, transportObjectState

  implicit none
  private

  ! Public procedures.
  public :: castCEParticleStatePtr

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
    procedure :: init
    procedure :: kill
    procedure :: preparePayload
    procedure :: setEnergy
  end type CEParticleState

contains
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
  subroutine init(self, payload)
    class(CEParticleState), intent(inout)               :: self
    class(buildTransportObjectStatePayload), intent(in) :: payload
    type(buildCEParticleStatePayload), pointer          :: payloadPtr
    character(*), parameter                             :: HERE = 'init (CEParticleState_class.f90)'

    ! Initialise superclass.
    call init_super(self, payload)

    ! Downcast payload to correct type.
    select type(ptr => payload)
      type is(buildCEParticleStatePayload)
        payloadPtr => ptr

      class default
        call fatalError(HERE, 'Invalid payload type.')

    end select

    ! Local.
    self % energy = payloadPtr % energy

  end subroutine init

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
    character(*), parameter                                :: HERE = 'preparePayload (CEParticleState_class.f90)'

    ! Superclass.
    call preparePayload_super(self, payload)

    ! Downcast payload to correct type.
    select type(ptr => payload)
      type is(buildCEParticleStatePayload)
        payloadPtr => ptr

      class default
        call fatalError(HERE, 'Invalid payload type.')

    end select

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