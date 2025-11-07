module CEPhotonState_class

  use CEParticleState_class,      only : CEParticleState
  use errors_mod,                 only : fatalError
  use numPrecision
  use transportObjectState_class, only : transportObjectState

  implicit none
  private

  ! Public procedures.
  public :: castCEPhotonStatePtr

  !!
  !!
  !!
  type, public, extends(CEParticleState) :: CEPhotonState
    private
  end type CEPhotonState

contains
  !!
  !!
  !!
  function castCEPhotonStatePtr(source, fatal) result(ptr)
    class(transportObjectState), intent(in) :: source
    logical(defBool), intent(in), optional  :: fatal
    type(CEPhotonState), pointer            :: ptr
    logical(defBool)                        :: throwError
    character(*), parameter                 :: here = 'castCEPhotonStatePtr (CEPhotonState_class.f90)'

    select type(temp => source)
      type is(CEPhotonState)
        ptr => temp

      class default
        ptr => null()

    end select

    ! Throw error if requested.
    throwError = .false.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) &
    call fatalError(here, "Transport object state is not of type 'CEPhotonState'.")

  end function castCEPhotonStatePtr

end module CEPhotonState_class