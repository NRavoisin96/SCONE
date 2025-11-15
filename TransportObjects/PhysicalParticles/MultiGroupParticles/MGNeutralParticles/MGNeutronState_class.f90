module MGNeutronState_class

  use errors_mod,                 only : fatalError
  use MGParticleState_class,      only : MGParticleState
  use numPrecision
  use transportObjectState_class, only : init_base_super => init_base, transportObjectState
  use universalVariables,         only : neutronMass

  implicit none
  private

  ! Public procedures.
  public :: castMGNeutronStatePtr

  !!
  !!
  !!
  type, public, extends(MGParticleState) :: MGNeutronState
    private
  contains
    procedure :: init_base
  end type MGNeutronState

contains
  !!
  !!
  !!
  function castMGNeutronStatePtr(source, fatal) result(ptr)
    class(transportObjectState), intent(in) :: source
    logical(defBool), intent(in), optional  :: fatal
    type(MGNeutronState), pointer           :: ptr
    logical(defBool)                        :: throwError
    character(*), parameter                 :: here = 'castMGNeutronStatePtr (MGNeutronState_class.f90)'

    select type(temp => source)
      type is(MGNeutronState)
        ptr => temp

      class default
        ptr => null()

    end select

    ! Throw error if requested.
    throwError = .false.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) &
    call fatalError(here, "Transport object state is not of type 'MGNeutronState'.")

  end function castMGNeutronStatePtr

  !!
  !!
  !!
  subroutine init_base(self)
    class(MGNeutronState), intent(inout) :: self

    ! Initialise superclass then set mass.
    call init_base_super(self)
    call self % setMass(neutronMass)

  end subroutine init_base

end module MGNeutronState_class