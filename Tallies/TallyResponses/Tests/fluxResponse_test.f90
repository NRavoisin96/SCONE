module fluxResponse_test

  use dictionary_class,          only : dictionary
  use fluxResponse_class,        only : fluxResponse
  use funit
  use numPrecision
  use testTransportObject_class, only : testTransportObject

  implicit none

@testCase
  type, extends(TestCase) :: test_fluxResponse
    private
    type(fluxResponse)        :: response
    type(testTransportObject) :: testObject
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_fluxResponse

contains
  !!
  !! Sets up test_fluxResponse object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_fluxResponse), intent(inout) :: this

    call this % testObject % init()

  end subroutine setUp

  !!
  !! Kills test_fluxResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_fluxResponse), intent(inout) :: this

    call this % testObject % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the response
  !!
@Test
  subroutine fluxResponseing(this)
    class(test_fluxResponse), intent(inout) :: this
    real(defReal)                           :: result
    real(defReal), parameter                :: TOL = 1.0e-9_defReal

    call this % response % get(this % testObject, result)
    @assertEqual(ONE, result, TOL)

  end subroutine fluxResponseing

end module fluxResponse_test