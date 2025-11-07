module MGParticleState_class

  use errors_mod,                  only : fatalError
  use numPrecision
  use physicalParticleState_class, only : display_super => display, kill_super => kill, physicalParticleState
  use transportObjectState_class,  only : transportObjectState

  implicit none
  private

  ! Public procedures.
  public :: castMGParticleStatePtr

  !!
  !!
  !!
  type, public, extends(physicalParticleState) :: MGParticleState
    private
    integer(shortInt) :: energyGroup = 0
  contains
    procedure :: display
    procedure :: getEnergyGroup
    procedure :: kill
    procedure :: setEnergyGroup
  end type MGParticleState

contains

  !!
  !!
  !!
  function castMGParticleStatePtr(source, fatal) result(ptr)
    class(transportObjectState), intent(in) :: source
    logical(defBool), intent(in), optional  :: fatal
    class(MGParticleState), pointer         :: ptr
    logical(defBool)                        :: throwError
    character(*), parameter                 :: here = 'castMGParticleStatePtr (MGParticleState_class.f90)'

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
    call fatalError(here, "Transport object state is not of type 'MGParticleState'.")

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
  elemental subroutine setEnergyGroup(self, energyGroup)
    class(MGParticleState), intent(inout) :: self
    integer(shortInt), intent(in)         :: energyGroup

    self % energyGroup = energyGroup

  end subroutine setEnergyGroup

end module MGParticleState_class