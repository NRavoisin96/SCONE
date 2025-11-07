module testTransportObject_class

  use numPrecision
  use transportObject_inter,      only : transportObject
  use transportObjectState_class, only : transportObjectState
  use universalVariables,         only : P_TEST_TRANSPORT_OBJECT

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(transportObject) :: testTransportObject
    private
  contains
    procedure :: allocateState
    procedure :: getType
  end type testTransportObject

contains
  !!
  !!
  !!
  subroutine allocateState(self, state)
    class(testTransportObject), intent(inout)             :: self
    class(transportObjectState), allocatable, intent(out) :: state

    allocate(transportObjectState :: state)

  end subroutine allocateState

  !!
  !!
  !!
  elemental function getType(self) result(type)
    class(testTransportObject), intent(in) :: self
    integer(shortInt)                      :: type

    type = P_TEST_TRANSPORT_OBJECT

  end function getType

end module testTransportObject_class