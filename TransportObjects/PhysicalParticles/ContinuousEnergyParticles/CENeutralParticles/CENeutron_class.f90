module CENeutron_class

  use CENeutronState_class,       only : castCENeutronStatePtr, CENeutronState
  use CEParticle_inter,           only : CEParticle
  use errors_mod,                 only : fatalError
  use numPrecision
  use physicalParticle_inter,     only : init_base_super => init_base
  use transportObject_inter,      only : transportObject
  use transportObjectState_class, only : transportObjectState
  use universalVariables,         only : P_NEUTRON_CE

  implicit none
  private

  ! Public procedures.
  public :: castCENeutronPtr

  !!
  !!
  !!
  type, public, extends(CEParticle) :: CENeutron
    private
  contains
    procedure :: allocateState
    procedure :: getType
  end type CENeutron

contains
  !!
  !!
  !!
  subroutine allocateState(self, state)
    class(CENeutron), intent(inout)                       :: self
    class(transportObjectState), allocatable, intent(out) :: state

    allocate(CENeutronState :: state)

  end subroutine allocateState

  !!
  !!
  !!
  function castCENeutronPtr(source, fatal) result(ptr)
    class(transportObject), intent(in)     :: source
    logical(defBool), intent(in), optional :: fatal
    type(CENeutron), pointer               :: ptr
    logical(defBool)                       :: throwError
    character(*), parameter                :: HERE = 'castCENeutronPtr (CENeutron_class.f90)'

    select type(temp => source)
      type is(CENeutron)
        ptr => temp

      class default
        ptr => null()

    end select

    ! Throw error if requested.
    throwError = .true.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) call fatalError(HERE, "Transport object is not of type 'CENeutron'.")

  end function castCENeutronPtr

  !!
  !!
  !!
  elemental function getType(self) result(type)
    class(CENeutron), intent(in) :: self
    integer(shortInt)            :: type

    type = P_NEUTRON_CE

  end function getType

end module CENeutron_class