module testResponse_test

  use dictionary_class,          only : dictionary
  use funit
  use numPrecision
  use testResponse_class,        only : testResponse
  use testTransportObject_class, only : testTransportObject

  implicit none

@testCase
  type, extends(TestCase) :: test_testResponse
    private
    type(testResponse)        :: response
    type(testTransportObject) :: testObject
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_testResponse


contains

  !!
  !! Sets up test_testResponse object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_testResponse), intent(inout) :: this
    type(dictionary)                        :: tempDict

    call tempDict % init(1)
    call tempDict % store('value', 1.3_defReal)
    call this % response % init(tempDict)

    call this % testObject % init()

  end subroutine setUp

  !!
  !! Kills test_testResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_testResponse), intent(inout) :: this

    call this % testObject % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the filter
  !!
@Test
  subroutine testResponseing(this)
    class(test_testResponse), intent(inout) :: this
    real(defReal)                           :: result
    real(defReal), parameter                :: TOL = 1.0e-9_defReal

    call this % response % get(this % testObject, result)
    @assertEqual(1.3_defReal, result, TOL)

  end subroutine testResponseing

end module testResponse_test