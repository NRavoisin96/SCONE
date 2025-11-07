module testFilter_test

  use dictionary_class,           only : dictionary
  use funit
  use numPrecision
  use testFilter_class,           only : testFilter
  use transportObjectState_class, only : transportObjectState

  implicit none

@testCase
  type, extends(TestCase) :: test_testFilter
    private
    type(testFilter) :: filter
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_testFilter


contains
@Before
  !!
  !! Sets up test_testFilter object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_testFilter), intent(inout) :: this
    type(dictionary)                      :: tempDict

    call tempDict % init(2)
    call tempDict % store('minIdx', 4)
    call tempDict % store('maxIdx', 6)

    call this % filter % init(tempDict)

  end subroutine setUp

@After
  !!
  !! Kills test_testFilter object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_testFilter), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the filter
  !!
@Test
  subroutine testFiltering(this)
    class(test_testFilter), intent(inout) :: this
    type(transportObjectState)            :: state
    logical(defBool)                      :: filterRes

    call state % setMaterialIdx(1)
    filterRes = this % filter % isPass(state)
    @assertFalse(filterRes)

    call state % setMaterialIdx(4)
    filterRes = this % filter % isPass(state)
    @assertTrue(filterRes)

    call state % setMaterialIdx(6)
    filterRes = this % filter % isPass(state)
    @assertTrue(filterRes)

    call state % setMaterialIdx(7)
    filterRes = this % filter % isPass(state)
    @assertFalse(filterRes)

  end subroutine testFiltering

end module testFilter_test