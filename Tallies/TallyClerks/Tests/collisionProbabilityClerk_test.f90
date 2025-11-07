module collisionProbabilityClerk_test

  use collisionProbabilityClerk_class, only : collisionProbabilityClerk, CPMResult
  use dictionary_class,                only : dictionary
  use funit
  use numPrecision
  use outputFile_class,                only : outputFile
  use particleDungeon_class,           only : particleDungeon
  use physicalParticleState_class,     only : castPhysicalParticleStatePtr, physicalParticleState
  use scoreMemory_class,               only : scoreMemory
  use tallyResult_class,               only : tallyResult
  use testNeutronDatabase_class,       only : testNeutronDatabase
  use testPhysicalParticle_class,      only : testPhysicalParticle

  implicit none

@testCase
  type, extends(TestCase)           :: test_collisionProbabilityClerk
    private
    type(collisionProbabilityClerk) :: clerk
    type(testPhysicalParticle)      :: testParticle
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_collisionProbabilityClerk

contains
@Before
  !!
  !! Sets up test_collisionProbabilityClerk object we can use in a number of tests
  !!
  !! Simple 2x2 collision probability score (with outside region) divided with test map
  !!
  subroutine setUp(this)
    class(test_collisionProbabilityClerk), intent(inout) :: this
    character(nameLen)                                   :: name
    type(dictionary)                                     :: dict, mapDict

    call mapDict % init(2)
    call mapDict % store('type', 'testMap')
    call mapDict % store('maxIdx', 2)

    ! Build intput dictionary
    call dict % init(2)
    call dict % store('type', 'collisionProbabilityClerk')
    call dict % store('map', mapDict)

    name = 'testClerk'
    call this % clerk % init(dict, name)

    call this % testParticle % init()

    call mapDict % kill()
    call dict % kill()

  end subroutine setUp

@After
  !!
  !! Kills test_collisionProbabilityClerk object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_collisionProbabilityClerk), intent(inout) :: this

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
    class(test_collisionProbabilityClerk), intent(inout) :: this
    class(tallyResult), allocatable                      :: res
    real(defReal)                                        :: val
    type(particleDungeon)                                :: pop
    type(physicalParticleState), pointer                 :: preCollisionStatePtr
    type(scoreMemory)                                    :: mem
    type(testNeutronDatabase)                            :: xsData
    real(defReal), parameter                             :: TOL = 1.0E-7

    ! Create score memory
    call mem % init(int(this % clerk % getSize(), longInt) , 1, batchSize = 1)
    call this % clerk % setMemAddress(1_longInt)

    ! Create test transport Nuclear Data
    call xsData % build(1.1_defReal, fissionXS = 1.1_defReal, nuFissionXS = 2.0_defReal)

    ! Create one particle that can be made to collide
    ! repeatedly in several materials

    ! Score some events
    ! Particle starts in material 2 and collides in material 2
    call this % testParticle % setMaterialIdx(2)
    call this % testParticle % setWeight(0.7_defReal)
    preCollisionStatePtr => castPhysicalParticleStatePtr(this % testParticle % getPreCollisionStatePtr(), .true.)
    call preCollisionStatePtr % setMaterialIdx(2)
    call this % clerk % reportInColl(this % testParticle, .false., xsData, mem)

    ! Particle starts in material 1 and collides in material 2
    call this % testParticle % setMaterialIdx(2)
    call this % testParticle % setWeight(1.1_defReal)
    call preCollisionStatePtr % setMaterialIdx(1)
    call this % clerk % reportInColl(this % testParticle, .false., xsData, mem)

    ! Particle starts in material 1 and collides in material 1
    call this % testParticle % setMaterialIdx(1)
    call this % testParticle % setWeight(ONE)
    call preCollisionStatePtr % setMaterialIdx(1)
    call this % clerk % reportInColl(this % testParticle, .false., xsData, mem)

    ! Particle starts in material 2 and collides in material 1
    call this % testParticle % setMaterialIdx(1)
    call this % testParticle % setWeight(1.4_defReal)
    call preCollisionStatePtr % setMaterialIdx(2)
    call this % clerk % reportInColl(this % testParticle, .false., xsData, mem)

    ! Particle starts in material 2 and collides in another, unknown material
    call this % testParticle % setMaterialIdx(7)
    call this % testParticle % setWeight(ONE)
    call preCollisionStatePtr % setMaterialIdx(2)
    call this % clerk % reportInColl(this % testParticle, .false., xsData, mem)

    ! Particle starts in an unknown material and collides in material 1
    call this % testParticle % setMaterialIdx(1)
    call this % testParticle % setWeight(0.9_defReal)
    call preCollisionStatePtr % setMaterialIdx(88)
    call this % clerk % reportInColl(this % testParticle, .false., xsData, mem)
    call this % clerk % reportCycleEnd(pop, mem)

    ! Close cycle
    call mem % closeCycle(ONE)

    ! Verify results

    ! Collision probability matrix

    ! outside -> outside Transition
    call mem % getResult(val, 1_longInt)
    @assertEqual(ZERO, val, TOL)

    ! outside -> 1 Transition
    call mem % getResult(val, 2_longInt)
    @assertEqual(ONE, val, TOL)

    ! outside -> 2 Transition
    call mem % getResult(val, 3_longInt)
    @assertEqual(ZERO, val, TOL)

    ! 1 -> outside Transition
    call mem % getResult(val, 4_longInt)
    @assertEqual(ZERO, val, TOL)

    ! 1 -> 1 Transition
    call mem % getResult(val, 5_longInt)
    @assertEqual(0.47619047619_defReal, val, TOL)

    ! 1 -> 2 Transition
    call mem % getResult(val, 6_longInt)
    @assertEqual(0.52380952381_defReal, val, TOL)

    ! 2 -> outside Transition
    call mem % getResult(val, 7_longInt)
    @assertEqual(0.32258064516_defReal, val, TOL)

    ! 2 -> 1 Transition
    call mem % getResult(val, 8_longInt)
    @assertEqual(0.45161290322_defReal ,val, TOL)

    ! 2 -> 2 Transition
    call mem % getResult(val, 9_longInt)
    @assertEqual(0.22580645161_defReal, val, TOL)

    ! Verify run-time result
    call this % clerk % getResult(res, mem)

    select type(res)
      class is (CPMresult)
        @assertEqual(3, res % N)

        ! outside -> outside Transition
        @assertEqual(ZERO, res % CPM(1, 1, 1), TOL)

        ! outside -> 1 Transition
        @assertEqual(ONE, res % CPM(2, 1, 1), TOL)

        ! outside -> 2 Transition
        @assertEqual(ZERO, res % CPM(3, 1, 1), TOL)

        ! 1 -> outside Transition
        @assertEqual(ZERO, res % CPM(1, 2, 1), TOL)

        ! 1 -> 1 Transition
        @assertEqual(0.47619047619_defReal, res % CPM(2, 2, 1), TOL)

        ! 1 -> 2 Transition
        @assertEqual(0.52380952381_defReal, res % CPM(3, 2, 1), TOL)

        ! 2 -> outside Transition
        @assertEqual(0.32258064516_defReal, res % CPM(1, 3, 1), TOL)

        ! 2 -> 1 Transition
        @assertEqual(0.45161290322_defReal, res % CPM(2, 3, 1), TOL)

        ! 2 -> 2 Transition
        @assertEqual(0.22580645161_defReal, res % CPM(3, 3, 1), TOL)

        ! Clean all entries
        res % CPM = ZERO

      class default
        @assertEqual(1, 2)

    end select

    ! Get result again -> verify correcness of reallocation logic by code coverage
    call this % clerk % getResult(res, mem)
    select type(res)
      class is (CPMresult)
        ! 1 -> 1 Transition
        @assertEqual(0.47619047619_defReal, res % CPM(2, 2, 1), TOL)

        ! Change size of matrix
        res % N = 2
        deallocate(res % CPM)
        allocate(res % CPM(2, 2, 1))

    end select

    ! Get result yet again to ensure the size was not incorrectly modified
    call this % clerk % getResult(res, mem)

    select type(res)
      class is (CPMresult)
        @assertEqual(3, res % N)
        @assertEqual([3, 3, 2], shape(res % CPM))
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
    class(test_collisionProbabilityClerk), intent(inout) :: this
    type(outputFile)                                     :: outF
    type(scoreMemory)                                    :: mem

    ! Create score memory
    call mem % init(int(this % clerk % getSize(), longInt), 1)
    call this % clerk % setMemAddress(1_longInt)

    ! Verify that output calls are correct
    call outF % init('dummyPrinter', fatalErrors = .false.)
    call this % clerk % print (outF, mem)

    @assertTrue(outF % isValid())

  end subroutine testPrintingCorrectness

end module collisionProbabilityClerk_test