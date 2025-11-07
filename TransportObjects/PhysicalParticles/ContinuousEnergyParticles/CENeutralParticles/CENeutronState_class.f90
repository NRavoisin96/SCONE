module CENeutronState_class

  use CEParticleState_class,      only : buildCEParticleStatePayload, CEParticleState
  use errors_mod,                 only : fatalError
  use numPrecision
  use transportObjectState_class, only : transportObjectState

  implicit none
  private

  ! Public procedures.
  public :: castCENeutronStatePtr, newCENeutronState

  !!
  !!
  !!
  type, public, extends(CEParticleState) :: CENeutronState
    private
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
  function newCENeutronState(payload) result(new)
    type(buildCEParticleStatePayload), intent(in) :: payload
    type(CENeutronState)                          :: new

    ! Initialise from payload.
    call new % init(payload)

  end function newCENeutronState

end module CENeutronState_class