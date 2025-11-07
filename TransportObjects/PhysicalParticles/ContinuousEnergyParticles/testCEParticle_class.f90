module testCEParticle_class

  use CEParticle_inter,           only : CEParticle
  use CEParticleState_class,      only : CEParticleState
  use errors_mod,                 only : fatalError
  use numPrecision
  use transportObject_inter,      only : transportObject
  use transportObjectState_class, only : transportObjectState
  use universalVariables,         only : P_TEST_TRANSPORT_OBJECT

  implicit none
  private

  ! Public procedures.
  public :: castTestCEParticlePtr

  !!
  !!
  !!
  type, public, extends(CEParticle) :: testCEParticle
    private
  contains
    procedure :: allocateState
    procedure :: getSpeed
    procedure :: getType
  end type testCEParticle

contains
  !!
  !!
  !!
  subroutine allocateState(self, state)
    class(testCEParticle), intent(inout)                  :: self
    class(transportObjectState), allocatable, intent(out) :: state

    allocate(CEParticleState :: state)

  end subroutine allocateState

  !!
  !!
  !!
  function castTestCEParticlePtr(source, fatal) result(ptr)
    class(transportObject), intent(in)     :: source
    logical(defBool), intent(in), optional :: fatal
    logical(defBool)                       :: throwError
    type(testCEParticle), pointer          :: ptr
    character(*), parameter                :: HERE = 'castTestCEParticlePtr (testCEParticle_class.f90)'

    select type(temp => source)
      type is(testCEParticle)
        ptr => temp

      class default
        ptr => null()

    end select

    ! Throw error if requested.
    throwError = .false.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) &
    call fatalError(HERE, "Transport object is not of type 'testCEParticle'.")

  end function castTestCEParticlePtr

  !!
  !!
  !!
  function getSpeed(self) result(speed)
    class(testCEParticle), intent(in) :: self
    real(defReal)                     :: speed

    speed = ZERO

  end function getSpeed

  !!
  !!
  !!
  elemental function getType(self) result(type)
    class(testCEParticle), intent(in) :: self
    integer(shortInt)                 :: type

    type = P_TEST_TRANSPORT_OBJECT

  end function getType

end module testCEParticle_class