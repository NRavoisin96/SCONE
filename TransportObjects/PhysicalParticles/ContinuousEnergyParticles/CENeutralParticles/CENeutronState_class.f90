module CENeutronState_class

  use CEParticleState_class,      only : buildCEParticleStatePayload, CEParticleState
  use errors_mod,                 only : fatalError
  use numPrecision
  use transportObjectState_class, only : init_base_super => init_base, transportObjectState
  use universalVariables,         only : neutronMass

  implicit none
  private

  ! Public procedures.
  public :: castCENeutronStatePtr, newCENeutronState

  !!
  !!
  !!
  type, public, extends(CEParticleState) :: CENeutronState
    private
  contains
    procedure :: init_base
  end type CENeutronState

contains
  !!
  !!
  !!
  function castCENeutronStatePtr(source, fatal) result(ptr)
    class(transportObjectState), intent(in) :: source
    logical(defBool), intent(in), optional  :: fatal
    type(CENeutronState), pointer           :: ptr
    logical(defBool)                        :: throwError
    character(*), parameter                 :: here = 'castCENeutronStatePtr (CENeutronState_class.f90)'

    select type(temp => source)
      type is(CENeutronState)
        ptr => temp

      class default
        ptr => null()

    end select

    ! Throw error if requested.
    throwError = .false.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) &
    call fatalError(here, "Transport object state is not of type 'CENeutronState'.")

  end function castCENeutronStatePtr

  !!
  !!
  !!
  subroutine init_base(self)
    class(CENeutronState), intent(inout) :: self

    ! Initialise superclass then set mass.
    call init_base_super(self)
    call self % setMass(neutronMass)

  end subroutine init_base

  !!
  !!
  !!
  function newCENeutronState(payload) result(new)
    type(buildCEParticleStatePayload), intent(in) :: payload
    type(CENeutronState)                          :: new

    ! Initialise from payload.
    call new % init(payload)

  end function newCENeutronState

end module CENeutronState_class