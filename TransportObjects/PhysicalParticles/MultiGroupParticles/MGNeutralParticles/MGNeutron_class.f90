module MGNeutron_class

  use errors_mod,                 only : fatalError
  use MGNeutronState_class,       only : MGNeutronState
  use MGParticle_inter,           only : MGParticle
  use numPrecision
  use physicalParticle_inter,     only : init_base_super => init_base
  use transportObject_inter,      only : transportObject
  use transportObjectState_class, only : transportObjectState
  use universalVariables,         only : P_NEUTRON_MG

  implicit none
  private

  ! Public procedures.
  public :: castMGNeutronPtr

  !!
  !!
  !!
  type, public, extends(MGParticle) :: MGNeutron
    private
  contains
    procedure :: allocateState
    procedure :: getType
  end type MGNeutron

contains
  !!
  !!
  !!
  subroutine allocateState(self, state)
    class(MGNeutron), intent(inout)                       :: self
    class(transportObjectState), allocatable, intent(out) :: state

    allocate(MGNeutronState :: state)

  end subroutine allocateState

  !!
  !!
  !!
  function castMGNeutronPtr(source, fatal) result(ptr)
    class(transportObject), intent(in)     :: source
    logical(defBool), intent(in), optional :: fatal
    logical(defBool)                       :: throwError
    type(MGNeutron), pointer               :: ptr
    character(*), parameter                :: here = 'castMGNeutronPtr (MGNeutron_class.f90)'

    select type(temp => source)
      type is(MGNeutron)
        ptr => temp

      class default
        ptr => null()

    end select

    throwError = .false.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) call fatalError(here, "Physical particle is not of type 'MGNeutron'.")

  end function castMGNeutronPtr

  !!
  !!
  !!
  elemental function getType(self) result(type)
    class(MGNeutron), intent(in) :: self
    integer(shortInt)            :: type

    type = P_NEUTRON_MG

  end function getType

end module MGNeutron_class