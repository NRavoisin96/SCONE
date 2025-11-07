module weightWindowsField_iTest
  
  use dictionary_class,         only : dictionary
  use dictParser_func,          only : charToDict
  use funit
  use numPrecision
  use testCEParticle_class,     only : testCEParticle
  use weightWindowsField_class, only : weightWindowsField

  implicit none

@testCase
  type, extends(TestCase)    :: test_weightWindows
    private
    type(testCEParticle)     :: p
    type(weightWindowsField) :: wwField
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_weightWindows

  !!
  !! Weight Windows Definition
  !!
  character(*), parameter :: DICT_DEF = &
  & "file ./IntegrationTestFiles/testWW ;"

contains
@Before
  !!
  !! Sets up test_weightWindows object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_weightWindows), intent(inout) :: this
    type(dictionary)                         :: dict

    call charToDict(dict, DICT_DEF)
    call this % p % init()
    call this % wwField % init(dict)

  end subroutine setUp

@After
  !!
  !! Kills test_weightWindows object
  !!
  subroutine tearDown(this)
    class(test_weightWindows), intent(inout) :: this

    call this % p % kill()
    call this % wwField % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test retrieving the weight window values
  !!
@Test
  subroutine testGetValue(this)
    class(test_weightWindows), intent(inout) :: this
    real(defReal), dimension(3)              :: bins, EXPECTED_BINS
    real(defReal), parameter                 :: TOL = 1.0e-6_defReal

    call this % p % setGlobalPosition([0.5_defReal, 7.0_defReal, ZERO])
    call this % p % setEnergy(10.0_defReal)

    bins = this % wwField % at(this % p)
    EXPECTED_BINS = [0.4_defReal, 1.5_defReal, 0.8_defReal]

    @assertEqual(EXPECTED_BINS, bins, tolerance = TOL)

    call this % p % setGlobalPosition([-0.5_defReal, 7.0_defReal, ZERO])
    call this % p % setEnergy(10.0_defReal)

    bins = this % wwField % at(this % p)
    EXPECTED_BINS = ZERO

    @assertEqual(EXPECTED_BINS, bins, tolerance = TOL)

  end subroutine testGetValue

end module weightWindowsField_iTest
