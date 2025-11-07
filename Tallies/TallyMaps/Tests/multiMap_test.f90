module multiMap_test
  
  use dictionary_class,           only : dictionary
  use dictParser_func,            only : charToDict
  use funit
  use multiMap_class,             only : multiMap
  use numPrecision
  use outputFile_class,           only : outputFile
  use transportObjectState_class, only : transportObjectState

  implicit none

@testCase
  type, extends(TestCase) :: test_multiMap
    type(multiMap) :: map
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_multiMap

contains
@Before
  !!
  !! Sets up test_intMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_multiMap), intent(inout) :: this
    type(dictionary)                    :: tempDict
    character(*), parameter :: def = "                                        &
      &type multiMap;                                                        &
      &maps (map1 map2 map3);                                                &
      &map1 {type spaceMap; axis x; grid unstruct; bins (0.0 1.0 2.0); }     &
      &map2 {type spaceMap; axis y; grid unstruct; bins (0.0 2.0 4.0 6.0); } &
      &map3 {type spaceMap; axis z; grid unstruct; bins (0.0 3.0 6.0); }     "

    call charToDict(tempDict, def)
    call this % map % init(tempDict)

    call tempDict % kill()

  end subroutine setUp

@After
  !!
  !! Kills test_intMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_multiMap), intent(inout) :: this

    call this % map % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test that map correctly returns number of bins and dimensions
  !!
@Test
  subroutine testBinAndDimension(this)
    class(test_multiMap), intent(inout) :: this

    ! Test approperiate dimension
    @assertEqual(3, this % map % dimensions())

    ! Test Bin number
    @assertEqual(12, this % map % bins(0))
    @assertEqual(2, this % map % bins(1))
    @assertEqual(3, this % map % bins(2))
    @assertEqual(2, this % map % bins(3))

    @assertEqual(0, this % map % bins(-1))
    @assertEqual(0, this % map % bins(4))

  end subroutine testBinAndDimension

  !!
  !! Test mapping
  !!
@Test
  subroutine testMapping(this)
    class(test_multiMap), intent(inout) :: this
    real(defReal), dimension(3)         :: r
    type(transportObjectState)          :: state

    ! Map to bin (1 1 1) -> idx == 1
    r = [0.1_defReal, 0.1_defReal, 0.1_defReal]
    call state % setGlobalPosition(r)
    @assertEqual(1, this % map % map(state))

    ! Map to bin (2 3 2) -> idx == 12
    r = [1.1_defReal, 5.1_defReal, 5.1_defReal]
    call state % setGlobalPosition(r)
    @assertEqual(12, this % map % map(state))

    ! Map to bin (2, 2, 1) -> idx == 4
    r = [1.1_defReal, 3.1_defReal, 2.1_defReal]
    call state % setGlobalPosition(r)
    @assertEqual(4, this % map % map(state))

    ! Map outside division -> idx == 0
    r = [-1.1_defReal, 5.1_defReal, 5.1_defReal]
    call state % setGlobalPosition(r)
    @assertEqual(0, this % map % map(state))

    r = [1.1_defReal, 50.1_defReal, 5.1_defReal]
    call state % setGlobalPosition(r)
    @assertEqual(0, this % map % map(state))

    r = [1.1_defReal, 5.1_defReal, -5.1_defReal]
    call state % setGlobalPosition(r)
    @assertEqual(0, this % map % map(state))

  end subroutine testMapping

  !!
  !! Test correctness of calls when printing
  !!
@Test
  subroutine testPrint(this)
    class(test_multiMap), intent(inout) :: this
    type(outputFile)                     :: out

    call out % init('dummyPrinter', fatalErrors = .false.)

    call this % map % print(out)
    @assertTrue(out % isValid(), 'Incorrect printing sequence: ')
    call out % reset()

  end subroutine testPrint

end module multiMap_test