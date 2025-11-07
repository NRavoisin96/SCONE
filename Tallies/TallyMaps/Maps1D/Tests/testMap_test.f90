module testMap_test
  
  use dictionary_class,           only : dictionary
  use funit
  use numPrecision
  use testMap_class,              only : testMap
  use transportObjectState_class, only : transportObjectState

  implicit none

@testCase
  type, extends(TestCase) :: test_testMap
    private
    type(testMap) :: map
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_testMap

contains
@Before
  !!
  !! Sets up test_intMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_testMap), intent(inout) :: this
    type(dictionary) :: tempDict

    call tempDict % init(1)
    call tempDict % store('maxIdx',4)
    call this % map % init(tempDict)

  end subroutine setUp

@After
  !!
  !! Kills test_intMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_testMap), intent(inout) :: this

    call this % map % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test that maps performs as expected
  !!
@Test
  subroutine testMapping(this)
    class(test_testMap), intent(inout) :: this
    type(transportObjectState)         :: state

    call state % setMaterialIdx(0)
    @assertEqual(0, this % map % map(state), 'Invalid idx case: ')

    call state % setMaterialIdx(1)
    @assertEqual(1, this % map % map(state), 'Normal idx case: ')

    call state % setMaterialIdx(2)
    @assertEqual(2, this % map % map(state), 'Normal idx case: ')

    call state % setMaterialIdx(4)
    @assertEqual(4, this % map % map(state), 'Normal idx case: ')

    call state % setMaterialIdx(5)
    @assertEqual(0, this % map % map(state), 'Invalid idx case: ')

    call state % setMaterialIdx(-1)
    @assertEqual(0, this % map % map(state), 'Invalid idx case: ')

  end subroutine testMapping

  !!
  !! Test getting number of bins and dimensions
  !!
@Test
  subroutine testNumBin(this)
    class(test_testMap), intent(inout) :: this

    ! Test number of bins
    @assertEqual(4, this % map % bins(0))
    @assertEqual(4, this % map % bins(1))
    @assertEqual(0, this % map % bins(2))
    @assertEqual(0, this % map % bins(-2))

    ! Test dimension
    @assertEqual(1, this % map % dimensions())

  end subroutine testNumBin

end module testMap_test