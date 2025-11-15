module keffImplicitClerk_test

  use dictionary_class,            only : dictionary
  use endfConstants
  use funit
  use keffAnalogClerk_class,       only : keffResult
  use keffImplicitClerk_class,     only : keffImplicitClerk
  use numPrecision
  use outputFile_class,            only : outputFile
  use particleDungeon_class,       only : particleDungeon
  use physicalParticle_inter,      only : physicalParticle
  use physicalParticleState_class, only : castPhysicalParticleStatePtr, physicalParticleState
  use scoreMemory_class,           only : scoreMemory
  use tallyCodes
  use tallyResult_class,           only : tallyResult
  use testNeutronDatabase_class,   only : testNeutronDatabase
  use testPhysicalParticle_class,  only : testPhysicalParticle

  implicit none

@testCase
  type, extends(TestCase) :: test_keffImplicitClerk
    private
    class(physicalParticle), allocatable :: testParticle
    type(keffImplicitClerk)              :: clerk
    type(testNeutronDatabase)            :: nucData
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_keffImplicitClerk

contains
@Before
  !!
  !! Sets up test_keffImplicitClerk object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_keffImplicitClerk), intent(inout) :: this
    character(nameLen), parameter                :: name = 'dummyClerk'
    type(dictionary)                             :: dict

    call dict % init(2)
    call this % clerk % init(dict, name)
    call dict % kill()

    call this % nucData % build(ONE, captureXS = 2.0_defReal, fissionXS = ONE, nuFissionXS = 3.0_defReal)
    allocate(testPhysicalParticle :: this % testParticle)
    call this % testParticle % init()

  end subroutine setUp

@After
  !!
  !! Kills test_keffImplicitClerk object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_keffImplicitClerk), intent(inout) :: this

    call this % clerk % kill()
    call this % nucData % kill()
    if (allocated(this % testParticle)) then
      call this % testParticle % kill()
      deallocate(this % testParticle)

    end if

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test for 1 cycle batches
  !!
@Test
  subroutine test1CycleBatch(this)
    class(test_keffImplicitClerk), intent(inout) :: this
    class(tallyResult), allocatable              :: res
    type(particleDungeon)                        :: pit
    type(physicalParticleState), pointer         :: preCollisionStatePtr
    type(scoreMemory)                            :: mem
    real(defReal), parameter                     :: TOL = 1.0E-9_defReal

    ! Configure memory
    call mem % init(10_longInt, 1)
    call this % clerk % setMemAddress(1_longInt)
    call pit % init(4)

    ! Configure particle
    call this % testParticle % setFate(LEAK_FATE)

    !*** Start cycle 1
    ! Score implicit reaction rates
    call this % testParticle % setWeight(0.7_defReal)
    call this % clerk % reportInColl(this % testParticle, .false., this % nucData, mem)

    ! Score analog production
    preCollisionStatePtr => castPhysicalParticleStatePtr(this % testParticle % getPreCollisionStatePtr())
    call preCollisionStatePtr % setWeight(0.1_defReal)
    call this % clerk % reportOutColl(this % testParticle, N_2N, 0.5_defReal, this % nucData, mem)

    ! Score leakage
    call this % testParticle % setWeight(0.3_defReal)
    call this % clerk % reportHist(this % testParticle, this % nucData, mem)

    ! End cycle
    call pit % detain(this % testParticle)
    call this % clerk % reportCycleEnd(pit, mem)
    call pit % release(this % testParticle)
    call mem % closeCycle(ONE)

    call this % testParticle % setFate(LEAK_FATE)

    !*** Start cycle 2
    ! Score implicit reaction rates
    call this % testParticle % setWeight(0.6_defReal)
    call this % clerk % reportInColl(this % testParticle, .false., this % nucData, mem)

    ! Score analog production
    preCollisionStatePtr => castPhysicalParticleStatePtr(this % testParticle % getPreCollisionStatePtr())
    call preCollisionStatePtr % setWeight(0.1_defReal)
    call this % clerk % reportOutColl(this % testParticle, N_2N, 0.5_defReal, this % nucData, mem)

    ! Score leakage
    call this % testParticle % setWeight(0.3_defReal)
    call this % clerk % reportHist(this % testParticle, this % nucData, mem)

    ! End cycle
    call pit % detain(this % testParticle)
    call this % clerk % reportCycleEnd(pit, mem)
    call pit % release(this % testParticle)
    call mem % closeCycle(ONE)

    ! Verify result
    call this % clerk % getResult(res, mem)
    select type(res)
      type is(keffResult)
        @assertEqual(0.906521739130435_defReal, res % keff(1), TOL, '1 Cycle Batch, keff from result: ')
        @assertEqual(0.006521739130435_defReal, res % keff(2), TOL, '1 Cycle Batch, keff STD from result: ')

      class default
        @assertTrue(.false., 'Result is not a keffResult.')

    end select

    ! Clean.
    call mem % kill()
    call pit % kill()

  end subroutine test1CycleBatch

  !!
  !! Test getSize() and print
  !!
@Test
  subroutine testMisc(this)
    class(test_keffImplicitClerk), intent(inout) :: this
    type(outputFile)                             :: out
    type(scoreMemory)                            :: mem

    ! Configure memory
    call mem % init(10_longInt, 1)
    call this % clerk % setMemAddress(1_longInt)
    call out % init('dummyPrinter', fatalErrors = .false.)

    ! Test getting size
    @assertEqual(5, this % clerk % getSize(), 'Test getSize(): ')

    ! Test correctness of output calls
    call this % clerk % print(out, mem)
    @assertTrue(out % isValid(), 'Test print(): ')

  end subroutine testMisc

end module keffImplicitClerk_test