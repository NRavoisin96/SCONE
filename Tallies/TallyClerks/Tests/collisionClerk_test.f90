module collisionClerk_test

  use ceNeutronCache_mod,         only : cache_init => init, trackingCache
  use collisionClerk_class,       only : collisionClerk
  use dictionary_class,           only : dictionary
  use errors_mod,                 only : fatalError
  use funit
  use genericProcedures,          only : numToChar
  use numPrecision
  use outputFile_class,           only : outputFile
  use scoreMemory_class,          only : scoreMemory
  use testNeutronDatabase_class,  only : testNeutronDatabase
  use testPhysicalParticle_class, only : testPhysicalParticle
  use transportObjectState_class, only : transportObjectState

  implicit none

  ! Parameters.
  real(defReal), parameter :: SCORE_1 = 0.7_defReal / 0.3_defReal, SCORE_2 = 1.3_defReal / 0.3_defReal

@testParameter(constructor = new_testNumber)
  type, extends(AbstractTestParameter) :: testNumber
    integer(shortInt) :: i = 0
  contains
    procedure :: toString
  end type testNumber

@testCase(constructor = newTest)
  type, extends(ParameterizedTestCase) :: test_collisionClerk
    private
    integer(longInt), dimension(:), allocatable :: bins
    logical(defBool)                            :: hasFilter = .false., hasMap = .false., has2Res = .false.
    real(defReal), dimension(:), allocatable    :: results
    type(collisionClerk)                        :: clerk
    type(dictionary)                            :: clerkDict, filterDict, mapDict, res1Dict, res2Dict
    type(outputFile)                            :: outF
    type(scoreMemory)                           :: mem
    type(testNeutronDatabase)                   :: nucData
    type(testPhysicalParticle)                  :: p
  contains
    procedure :: initTest
    procedure :: setUp
    procedure :: tearDown
    procedure :: verifyResults
  end type test_collisionClerk

contains
  !!
  !! Build new test parameter form integer
  !!
  function new_testNumber(i) result (tstNum)
    integer(shortInt) :: i
    type(testNumber)  :: tstNum

    tstNum % i = i

  end function new_testNumber

  !!
  !! Write test parameter to string
  !!
  function toString(this) result(string)
    class(testNumber), intent(in) :: this
    character(:), allocatable :: string
    character(nameLen)        :: str

    write (str, *) this % i
    string = str

  end function toString

  !!
  !! Construct test case
  !!
  function newTest(testParam) result(tst)
    type(testNumber), intent(in) :: testParam
    type(test_collisionClerk)    :: tst
    integer(shortInt)            :: i, nBins, testNum

    ! Set test parameters
    testNum = testParam % i - 1
    tst % hasFilter = mod(testNum, 2) == 1
    tst % hasMap = mod(testNum / 2, 2) == 1
    tst % has2Res = mod(testNum / 4, 2) == 1

    nBins = 1
    if (tst % has2Res) nBins = 2 * nBins
    if (tst % hasMap) nBins = 7 * nBins

    allocate(tst % bins(nBins), tst % results(nBins))
    tst % bins = [(int(i, longInt), i = 1, nBins)]
    tst % results = ZERO

    select case(testParam % i)
      case(1, 5)
        tst % results(1) = SCORE_1 + SCORE_2
        if (testParam % i == 5) tst % results(2) = 1.3_defReal * (SCORE_1 + SCORE_2)
          
      case(2, 4)
        tst % results(1) = SCORE_1
      
      case(3)
        tst % results(1) = SCORE_1
        tst % results(6) = SCORE_2
      
      case(6, 7, 8)
        tst % results(1) = SCORE_1
        tst % results(2) = 1.3_defReal * SCORE_1
        if (testParam % i == 7) then
          tst % results(11) = SCORE_2
          tst % results(12) = 1.3_defReal * SCORE_2

        end if

    end select

  end function newTest

@Before
  !!
  !!
  !!
  subroutine setUp(this)
    class(test_collisionClerk), intent(inout) :: this

    ! Initialise components.
    call this % clerkDict % init(7)
    call this % filterDict % init(3)
    call this % mapDict % init(2)
    call this % res1Dict % init(1)
    call this % res2Dict % init(2)
    call this % nucData % build(0.3_defReal)
    call this % p % init()
    call this % outF % init('dummyPrinter', fatalErrors = .false.)

    ! Build all common dictionary and data objects.
    call this % filterDict % store('type', 'testFilter')
    call this % filterDict % store('minIdx', 0)
    call this % filterDict % store('maxIdx', 5)

    call this % mapDict % store('type', 'testMap')
    call this % mapDict % store('maxIdx', 7)

    call this % res1Dict % store('type', 'fluxResponse')

    call this % res2Dict % store('type', 'testResponse')
    call this % res2Dict % store('value', 1.3_defReal)

  end subroutine setUp

@After
  !!
  !!
  !!
  subroutine tearDown(this)
    class(test_collisionClerk), intent(inout) :: this

    if (allocated(this % bins)) deallocate(this % bins)
    this % hasFilter = .false.
    this % hasMap = .false.
    this % has2Res = .false.
    if (allocated(this % results)) deallocate(this % results)
    call this % p % kill()
    call this % nucData % kill()
    call this % clerkDict % kill()
    call this % filterDict % kill()
    call this % mapDict % kill()
    call this % res1Dict % kill()
    call this % res2Dict % kill()
    call this % mem % kill()
    call this % clerk % kill()
    call this % outF % reset()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!!  PRIVATE HELPER ROUTINES
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !!
  !!
  !!
  subroutine initTest(this, handleVirtual, caseDescription)
    class(test_collisionClerk), intent(inout) :: this
    logical(defBool), intent(in)              :: handleVirtual
    character(:), allocatable, intent(out)    :: caseDescription
    character(nameLen), parameter             :: NAME = 'myClerk', RES_1_NAME = 'flux', RES_2_NAME = 'testResponse'

    ! Build configuration for the specific test run.
    caseDescription = 'Vanilla case with:'
    if (this % hasFilter) caseDescription = caseDescription//' Filter'
    if (this % hasMap) caseDescription = caseDescription//' Map'
    if (this % has2Res) caseDescription = caseDescription//' 2nd Response'

    call this % clerkDict % store('type', 'collisionClerk')
    call this % clerkDict % store('handleVirtual', merge(1, 0, handleVirtual))
    call this % clerkDict % store(RES_1_NAME, this % res1Dict)
    call this % clerkDict % store(RES_2_NAME, this % res2Dict)

    if (this % hasFilter) call this % clerkDict % store('filter', this % filterDict)
    if (this % hasMap) call this % clerkDict % store('map', this % mapDict)

    if (this % has2Res) then
      call this % clerkDict % store('response', [RES_1_NAME, RES_2_NAME])

    else
      call this % clerkDict % store('response', [RES_1_NAME])

    end if

    call this % clerk % init(this % clerkDict, NAME)
    call this % clerk % setMemAddress(1_longInt)

    call this % mem % init(int(this % clerk % getSize(), longInt), 1)

  end subroutine initTest

  !!
  !!
  !!
  subroutine verifyResults(this, caseDescription)
    class(test_collisionClerk), intent(inout) :: this
    character(*), intent(in)                  :: caseDescription
    integer(shortInt)                         :: i
    real(defReal)                             :: res
    real(defReal), parameter                  :: TOL = 1.0e-9_defReal

    ! Close cycle.
    call this % mem % closeCycle(ONE)

    ! Verify results of scoring
    do i = 1, size(this % bins)
      call this % mem % getResult(res, this % bins(i))
      @assertEqual(this % results(i), res, TOL, caseDescription//' BIN: '//numToChar(i))

    end do

    ! Verify that size of memory returned is correct
    @assertEqual(size(this % bins), this % clerk % getSize(), caseDescription//' Memory size test: ')

    ! Verify that output calls are correct
    call this % clerk % print (this % outF, this % mem)
    @assertTrue(this % outF % isValid(), caseDescription)

  end subroutine verifyResults

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !!
  !! Scoring test
  !!
@Test(cases = [1, 2, 3, 4, 5, 6, 7, 8])
  subroutine testScoring(this)
    class(test_collisionClerk), intent(inout) :: this
    character(:), allocatable                 :: caseDescription

    ! Initialise test.
    call this % initTest(.false., caseDescription)

    ! Perform scoring operations.
    call this % p % setMaterialIdx(1)
    call this % p % setWeight(0.7_defReal)
    call this % clerk % reportInColl(this % p, .false., this % nucData, this % mem)
    
    call this % p % setMaterialIdx(6)
    call this % p % setWeight(1.3_defReal)
    call this % clerk % reportInColl(this % p, .false., this % nucData, this % mem)

    ! Virtual scoring should not contribute to score in this case.
    call this % p % setWeight(1000.3_defReal)
    call this % clerk % reportInColl(this % p, .true., this % nucData, this % mem)

    ! Verify results.
    call this % verifyResults(caseDescription)

  end subroutine testScoring

@Test(cases = [1, 2, 3, 4, 5, 6, 7, 8])
  subroutine testScoringVirtual(this)
    class(test_collisionClerk), intent(inout) :: this
    character(:), allocatable                 :: caseDescription

    ! Initialise test.
    call this % initTest(.true., caseDescription)

    ! Perform scoring operations.
    call this % p % setMaterialIdx(1)
    call this % p % setWeight(0.7_defReal)
    call this % clerk % reportInColl(this % p, .true., this % nucData, this % mem)
    
    call this % p % setMaterialIdx(6)
    call this % p % setWeight(1.3_defReal)
    call this % clerk % reportInColl(this % p, .false., this % nucData, this % mem)

    ! Verify results.
    call this % verifyResults(caseDescription)

  end subroutine testScoringVirtual

end module collisionClerk_test