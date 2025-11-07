module simpleFMClerk_test

  use dictionary_class,            only : dictionary
  use funit
  use numPrecision
  use outputFile_class,            only : outputFile
  use particleDungeon_class,       only : particleDungeon
  use physicalParticleState_class, only : castPhysicalParticleStatePtr, physicalParticleState
  use scoreMemory_class,           only : scoreMemory
  use simpleFMClerk_class,         only : simpleFMClerk, FMResult
  use tallyResult_class,           only : tallyResult
  use testNeutronDatabase_class,   only : testNeutronDatabase
  use testPhysicalParticle_class,  only : testPhysicalParticle

  implicit none

@testCase
  type, extends(TestCase) :: test_simpleFMClerk
    private
    type(simpleFMClerk)        :: clerk
    type(testPhysicalParticle) :: testParticle
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_simpleFMClerk

contains
@Before
  !!
  !! Sets up test_simpleFMClerk object we can use in a number of tests
  !!
  !! Simple 3x3 fission matrix divided with test map
  !!
  subroutine setUp(this)
    class(test_simpleFMClerk), intent(inout) :: this
    character(nameLen)                       :: name
    type(dictionary)                         :: dict, mapDict

    call mapDict % init(2)
    call mapDict % store('type','testMap')
    call mapDict % store('maxIdx',3)

    ! Build intput dictionary
    call dict % init(2)
    call dict % store('type','simpleFMClerk')
    call dict % store('map', mapDict)

    name = 'testClerk'
    call this % clerk % init(dict,name)
    call this % testParticle % init()

    call mapDict % kill()
    call dict % kill()

  end subroutine setUp

@After
  !!
  !! Kills test_simpleFMClerk object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_simpleFMClerk), intent(inout) :: this

    call this % clerk % kill()
    call this % testParticle % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correctness in a simple use case
  !!
@Test
  subroutine testSimpleUseCase(this)
    class(test_simpleFMClerk), intent(inout) :: this
    class(tallyResult), allocatable          :: res
    real(defReal)                            :: val
    type(particleDungeon)                    :: pop
    type(physicalParticleState)              :: phase
    type(physicalParticleState), pointer     :: preHistoryStatePtr
    type(scoreMemory)                        :: mem
    type(testNeutronDatabase)                :: xsData
    real(defReal), parameter                 :: TOL = 1.0e-7_defReal

    ! Create score memory
    call mem % init(int(this % clerk % getSize(), longInt), 1, batchSize = 1)
    call this % clerk % setMemAddress(1_longInt)

    ! Create test transport Nuclear Data
    call xsData % build(1.1_defReal, fissionXS = 1.1_defReal, nuFissionXS = 2.0_defReal)

    ! Crate dungeon of original events
    ! One particle born in matIdx 1 and other in 2
    call pop % init(3)

    call phase % setMaterialIdx(2)
    call pop % detain(phase)

    call phase % setMaterialIdx(1)
    call pop % detain(phase)

    call this % clerk % reportCycleStart(pop, mem)

    ! Score some events
    call this % testParticle % setMaterialIdx(2)
    call this % testParticle % setWeight(0.7_defReal)
    preHistoryStatePtr => castPhysicalParticleStatePtr(this % testParticle % getPreHistoryStatePtr(), .true.)
    call preHistoryStatePtr % setMaterialIdx(2)
    call this % clerk % reportInColl(this % testParticle, .false., xsData, mem)

    call this % testParticle % setMaterialIdx(1)
    call this % testParticle % setWeight(1.1_defReal)
    call preHistoryStatePtr % setMaterialIdx(2)
    call this % clerk % reportInColl(this % testParticle, .false., xsData, mem)


    call this % testParticle % setMaterialIdx(1)
    call this % testParticle % setWeight(ONE)
    call preHistoryStatePtr % setMaterialIdx(1)
    call this % clerk % reportInColl(this % testParticle, .false., xsData, mem)

    call this % clerk % reportCycleEnd(pop, mem)

    ! Close cycle
    call mem % closeCycle(ONE)

    ! Verify results

    ! Fission matrix
    ! 1 -> 1 Transition
    call mem % getResult(val, 1_longInt)
    @assertEqual(1.818181818181_defReal, val, TOL)

    ! 1 -> 2 Transition
    call mem % getResult(val, 2_longInt)
    @assertEqual(ZERO, val, TOL)

    ! 1 -> 3 Transition
    call mem % getResult(val, 3_longInt)
    @assertEqual(ZERO, val, TOL)

    ! 2 -> 1 Transition
    call mem % getResult(val, 4_longInt)
    @assertEqual(2.0_defReal, val, TOL)

    ! 2 -> 2 Transition
    call mem % getResult(val, 5_longInt)
    @assertEqual(1.27272727272727_defReal, val, TOL)

    ! Verify run-time result
    call this % clerk % getResult(res, mem)

    select type(res)
      class is (FMresult)
        @assertEqual(3, res % N)

        ! 1 -> 1 Transition
        @assertEqual(1.818181818181_defReal, res % FM(1, 1, 1) , TOL)

        ! 1 -> 2 Transition
        @assertEqual(ZERO, res % FM(2, 1, 1), TOL)

        ! 1 -> 3 Transition
        @assertEqual(ZERO, res % FM(3, 1, 1), TOL)

        ! 2 -> 1 Transition
        @assertEqual(2.0_defReal, res % FM(1, 2, 1), TOL)

        ! 2 -> 2 Transition
        @assertEqual(1.27272727272727_defReal, res % FM(2, 2, 1), TOL)

        ! Clean all entries
        res % FM = ZERO

      class default
        @assertEqual(1,2)

    end select

    ! Get result again -> verify correcness of reallocation logic by code coverage
    call this % clerk % getResult(res, mem)
    select type(res)
      class is (FMresult)
        ! 1 -> 1 Transition
        @assertEqual(1.818181818181_defReal, res % FM(1, 1, 1), TOL)

        ! Change size of matrix
        res % N = 2
        deallocate(res % FM)
        allocate(res % FM(2, 2, 1))

    end select
    ! Get result yet again. This time with wrong size
    call this % clerk % getResult(res, mem)

    select type(res)
      class is (FMresult)
        @assertEqual(3, res % N)
        @assertEqual([3, 3, 2], shape(res % FM))
    end select

    ! Clean
    call xsData % kill()
    call pop % kill()

  end subroutine testSimpleUseCase

  !!
  !! Test correctness of the printing calls
  !!
@Test
  subroutine testPrintingCorrectness(this)
    class(test_simpleFMClerk), intent(inout) :: this
    type(outputFile)                         :: outF
    type(scoreMemory)                        :: mem

    ! Create score memory
    call mem % init(int(this % clerk % getSize(), longInt), 1)
    call this % clerk % setMemAddress(1_longInt)

    ! Verify that output calls are correct
    call outF % init('dummyPrinter', fatalErrors = .false.)
    call this % clerk % print(outF, mem)

    @assertTrue(outF % isValid())

  end subroutine testPrintingCorrectness

end module simpleFMClerk_test