module randomWalker_class

  use errors_mod,                 only : fatalError
  use numPrecision
  use transportObject_inter,      only : transportObject
  use transportObjectState_class, only : transportObjectState
  use universalVariables,         only : P_RANDOM_WALKER
  use virtualParticle_inter,      only : virtualParticle

  implicit none
  private

  ! Public procedures.
  public :: castRandomWalkerPtr, newRandomWalker

  !!
  !!
  !!
  type, public, extends(virtualParticle) :: randomWalker
    private
    real(defReal) :: accumulatedValue = ZERO
  contains
    procedure :: accumulateValue
    procedure :: allocateState
    procedure :: getAccumulatedValue
    procedure :: getType
    procedure :: setAccumulatedValue
  end type randomWalker

contains
  !!
  !!
  !!
  subroutine allocateState(self, state)
    class(randomWalker), intent(inout)                    :: self
    class(transportObjectState), allocatable, intent(out) :: state

    allocate(transportObjectState :: state)

  end subroutine allocateState

  !!
  !!
  !!
  elemental subroutine accumulateValue(self, value)
    class(randomWalker), intent(inout) :: self
    real(defReal), intent(in)          :: value

    self % accumulatedValue = self % accumulatedValue + value

  end subroutine accumulateValue

  !!
  !!
  !!
  function castRandomWalkerPtr(source, fatal) result(ptr)
    class(transportObject), intent(in)     :: source
    logical(defBool), intent(in), optional :: fatal
    logical(defBool)                       :: throwError
    type(randomWalker), pointer            :: ptr
    character(*), parameter                :: here = 'castRandomWalkerPtr (randomWalker_class.f90)'

    select type(temp => source)
      type is(randomWalker)
        ptr => temp

      class default
        ptr => null()

    end select

    ! Throw error if requested.
    throwError = .false.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) call fatalError(here, "Transport object is not of type 'randomWalker'.")

  end function castRandomWalkerPtr

  !!
  !!
  !!
  elemental function getAccumulatedValue(self) result(accumulatedValue)
    class(randomWalker), intent(in) :: self
    real(defReal)                   :: accumulatedValue

    accumulatedValue = self % accumulatedValue

  end function getAccumulatedValue

  !!
  !!
  !!
  elemental function getType(self) result(type)
    class(randomWalker), intent(in) :: self
    integer(shortInt)               :: type

    type = P_RANDOM_WALKER

  end function getType

  !!
  !!
  !!
  function newRandomWalker() result(new)
    type(randomWalker) :: new

  end function newRandomWalker

  !!
  !!
  !!
  elemental subroutine setAccumulatedValue(self, value)
    class(randomWalker), intent(inout) :: self
    real(defReal), intent(in)          :: value

    self % accumulatedValue = value

  end subroutine setAccumulatedValue

end module randomWalker_class