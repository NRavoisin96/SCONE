module particleDungeon_test
  
  use CEParticleState_class,       only : CEParticleState
  use funit
  use MGParticleState_class,       only : castMGParticleStatePtr, MGParticleState
  use numPrecision
  use particleDungeon_class,       only : particleDungeon
  use physicalParticle_inter,      only : physicalParticle
  use physicalParticleState_class, only : castPhysicalParticleStatePtr, physicalParticleState
  use RNG_class,                   only : RNG
  use testCEParticle_class,        only : castTestCEParticlePtr, testCEParticle
  use testMGParticle_class,        only : testMGParticle
  use testPhysicalParticle_class,  only : testPhysicalParticle

  implicit none

contains

  !!
  !! Test stack like access. Test is a dummy use case
  !!
@Test
  subroutine testStackInterface()
    class(physicalParticle), allocatable :: p
    integer(shortInt)                    :: i
    type(CEParticleState)                :: phase
    type(particleDungeon)                :: dungeon
    type(testCEParticle), pointer        :: testCEParticlePtr

    ! Is empty result for uninitialised dungeon
    @assertTrue(dungeon % isEmpty())

    ! Initialise
    call dungeon % init(10)

    ! Is empty initialised. No particles stored
    @assertTrue(dungeon % isEmpty())
    @assertEqual(0, dungeon % popSize())

    ! Store some particles with energy
    allocate(testCEParticle :: p)
    testCEParticlePtr => castTestCEParticlePtr(p, .true.)
    call testCEParticlePtr % init()
    call testCEParticlePtr % setEnergy(8.6_defReal)
    do i = 1, 5
      call dungeon % detain(testCEParticlePtr)

    end do

    ! Store some phase coordinates with energy
    call phase % setEnergy(3.4_defReal)
    do i = 1, 4
      call dungeon % detain(phase)

    end do

    ! Verify size
    @assertFalse(dungeon % isEmpty())
    @assertEqual(9, dungeon % popSize())

    ! Remove particles
    do i = 1, 4
      call dungeon % release(p)
      testCEParticlePtr => castTestCEParticlePtr(p, .true.)
      @assertEqual(3.4_defReal, testCEParticlePtr % getEnergy())

    end do

     do i = 1, 5
      call dungeon % release(p)
      testCEParticlePtr => castTestCEParticlePtr(p, .true.)
      @assertEqual(8.6_defReal, testCEParticlePtr % getEnergy())

    end do

    ! Verify size
    @assertTrue(dungeon % isEmpty())
    @assertEqual(0, dungeon % popSize())

    ! Clean.
    call p % kill()
    deallocate(p)
    call dungeon % kill()

  end subroutine testStackInterface

  !!
  !! Test array like interface
  !!
@Test
  subroutine testArrayInterface
    class(physicalParticle), allocatable :: p
    integer(shortInt)                    :: i
    type(CEParticleState)                :: phase
    type(particleDungeon)                :: dungeon
    type(testCEParticle), pointer        :: testCEParticlePtr

    ! Is empty result for uninitialised dungeon
    @assertTrue(dungeon % isEmpty())

    ! Initialise to fixed size
    call dungeon % setSize(2)

    ! Is filled with random particles
    @assertFalse(dungeon % isEmpty())
    @assertEqual(2, dungeon % popSize())

    ! Extend size
    call dungeon % setSize(9)
    @assertFalse(dungeon % isEmpty())
    @assertEqual(9, dungeon % popSize())

    ! Store some particles with energy
    allocate(testCEParticle :: p)
    testCEParticlePtr => castTestCEParticlePtr(p, .true.)
    call testCEParticlePtr % init()
    call testCEParticlePtr % setEnergy(8.6_defReal)
    do i = 1, 5
      call dungeon % replace(testCEParticlePtr, i)

    end do

    ! Store some phase coordinates with energy
    call phase % setEnergy(3.4_defReal)
    do i = 1, 4
      call dungeon % replace(phase, 5 + i)

    end do

    ! Verify size
    @assertFalse(dungeon % isEmpty())
    @assertEqual(9, dungeon % popSize())

    ! Raplace by particle and phaseCoords
    call testCEParticlePtr % setEnergy(13.0_defReal)
    call phase % setEnergy(17.0_defReal)
    call dungeon % replace(testCEParticlePtr, 9)
    call dungeon % replace(phase, 1)

    ! Verify particles by copies
    p = dungeon % copy(9)
    testCEParticlePtr => castTestCEParticlePtr(p, .true.)
    @assertEqual(13.0_defReal, testCEParticlePtr % getEnergy())

    do i = 8, 6, -1
      p = dungeon % copy(i)
      testCEParticlePtr => castTestCEParticlePtr(p, .true.)
      @assertEqual(3.4_defReal, testCEParticlePtr % getEnergy())

    end do

    do i = 5, 2, -1
      p = dungeon % copy(i)
      testCEParticlePtr => castTestCEParticlePtr(p, .true.)
      @assertEqual(8.6_defReal, testCEParticlePtr % getEnergy())

    end do
    p = dungeon % copy(1)
    testCEParticlePtr => castTestCEParticlePtr(p, .true.)
    @assertEqual(17.0_defReal, testCEParticlePtr % getEnergy())

    ! Verify that population has not changed
    @assertFalse(dungeon % isEmpty())
    @assertEqual(9, dungeon % popSize())

    ! Shrink size
    call dungeon % setSize(2)
    @assertFalse(dungeon % isEmpty())
    @assertEqual(2, dungeon % popSize())

    ! Clean population
    call dungeon % cleanPop()
    @assertTrue(dungeon % isEmpty())
    @assertEqual(0, dungeon % popSize())

    ! Clean memory
    call dungeon % kill()

  end subroutine testArrayInterface

  !!
  !! Test weight normalisation and inquiry
  !!
@Test
  subroutine testWeightNorm()
    class(physicalParticle), allocatable :: testParticle
    integer(shortInt)                    :: i
    type(particleDungeon)                :: dungeon
    real(defReal), parameter             :: TOL = 1.0e-9_defReal

    ! Initialise
    call dungeon % init(10)
    @assertEqual(ZERO, dungeon % popWeight(), TOL)

    ! Store some particles with non-uniform weight
    allocate(testPhysicalParticle :: testParticle)
    call testParticle % init()
    do i = 1, 5
      call testParticle % setWeight(0.5_defReal + 0.1_defReal * i)
      call dungeon % detain(testParticle)

    end do

    ! Verify total weight
    @assertEqual(4.0_defReal, dungeon % popWeight(), TOL)

    ! Normalise weight
    call dungeon % normWeight(2.0_defReal)
    @assertEqual(2.0_defReal, dungeon % popWeight(), TOL)

    ! Get particles and compare weight
    do i = 5, 1, -1
      call dungeon % release(testParticle)
      @assertEqual(0.25_defReal + 0.05_defReal * i, testParticle % getWeight(), TOL)

    end do

    ! Verify weight
    @assertEqual(0.0_defReal, dungeon % popWeight(), TOL)

    ! Clean
    call dungeon % kill()

  end subroutine testWeightNorm

  !!
  !! Test normalisation of population to smaller number
  !! Particles with non-uniform weight
  !!  NOTE: Weight preservation is disabled for now
  !!
@Test
  subroutine testNormPopDown()
    class(physicalParticle), allocatable :: testParticle
    integer(shortInt)                    :: i
    type(particleDungeon)                :: dungeon
    type(RNG)                            :: pRNG
    real(defReal), parameter             :: TOL = 1.0e-9_defReal


    ! Initialise
    call dungeon % init(10)
    call pRNG % init(7865856_longInt)

    ! Store some particles with non-uniform weight
    allocate(testPhysicalParticle :: testParticle)
    call testParticle % init()
    do i = 1, 10
      call testParticle % setWeight(0.5_defReal + 0.1_defReal * i)
      call testParticle % setBroodId(1) ! Avoid triggering error on sort by broodID
      call dungeon % detain(testParticle)

    end do

    ! Normalise population
    call dungeon % normSize(5, pRNG)

    ! Verify size
    @assertEqual(5, dungeon % popSize())

    ! Verify weight *** DISABLED
    !@assertEqual(6.05_defReal, dungeon % popWeight(), TOL)

    ! Clean memory
    call dungeon % kill()

  end subroutine testNormPopDown

  !!
  !! Test normalisation of population to smaller number
  !! Particles with non-uniform weight
  !!  NOTE: Weight preservation is disabled for now
  !!
@Test
  subroutine testNormPopUp()
    class(physicalParticle), allocatable :: testParticle
    integer(shortInt)                    :: i
    type(particleDungeon)                :: dungeon
    type(RNG)                            :: pRNG
    real(defReal), parameter             :: TOL = 1.0e-9_defReal

    ! Initialise
    call dungeon % init(20)
    call pRNG % init(435468_longInt)

    ! Store some particles with non-uniform weight
    allocate(testPhysicalParticle :: testParticle)
    call testParticle % init()
    do i = 1, 10
      call testParticle % setWeight(0.5_defReal + 0.1_defReal * i)
      call testParticle % setBroodId(1) ! Avoid triggering error on sort by broodID
      call dungeon % detain(testParticle)

    end do

    ! Normalise population
    call dungeon % normSize(15, pRNG)

    ! Verify size
    @assertEqual(15, dungeon % popSize())

    ! Verify weight *** DISABLED
    !@assertEqual(18.15_defReal, dungeon % popWeight(), TOL)

    ! Clean memory
    call dungeon % kill()

  end subroutine testNormPopUp

  !!
  !! Test sorting of the population by brood ID without duplicates
  !!
@Test
  subroutine testSortingByBroodID()
    class(physicalParticle), allocatable :: testParticle
    integer(shortInt)                    :: i
    type(particleDungeon)                :: dungeon
    type(physicalParticleState), pointer :: statePtr
    real(defReal), parameter             :: TOL = 1.0e-9_defReal

    ! Initialise
    call dungeon % init(10)

    ! Store some particles with brood ID in reverse order
    allocate(testPhysicalParticle :: testParticle)
    call testParticle % init()
    do i = 1, 10
      call testParticle % setBroodId(10 - i + 1)
      call dungeon % detain(testParticle)

    end do

    ! Sort by brood ID
    call dungeon % sortByBroodID(10)

    ! Verify order
    do i = 1, 10
      statePtr => castPhysicalParticleStatePtr(dungeon % get(i))
      @assertEqual(i, statePtr % getBroodId())

    end do

    call testParticle % kill()
    deallocate(testParticle)
    call dungeon % kill()

  end subroutine testSortingByBroodID


  !!
  !! Test sorting of the population by brood ID with duplicates
  !!
  @Test
  subroutine testSortingByBroodID_withDuplicates()
    class(physicalParticleState), allocatable :: testParticleState
    integer(shortInt)                         :: i, j
    type(MGParticleState), pointer            :: MGParticleStatePtr
    type(particleDungeon)                     :: dungeon
    integer(shortInt), parameter              :: N_duplicates = 7
    real(defReal), parameter                  :: TOL = 1.0e-9_defReal

    ! Initialise
    call dungeon % init(10 * N_duplicates)

    ! Store some particles with brood ID in reverse order
    ! Use the group number to distinguish duplicates and make sure
    ! that the insertion order is preserved (for particles with the same brood ID)
    allocate(MGParticleState :: testParticleState)
    MGParticleStatePtr => castMGParticleStatePtr(testParticleState, .true.)
    do j = 1, N_duplicates
      do i = 1, 10
        call MGParticleStatePtr % setBroodId(10 - i + 1)
        call MGParticleStatePtr % setEnergyGroup(j)
        call dungeon % detain(MGParticleStatePtr)

      end do

    end do

    ! Sort by brood ID
    call dungeon % sortByBroodID(10)

    ! Verify order
    do i = 1, 10
      do j = 1, N_duplicates
        MGParticleStatePtr => castMGParticleStatePtr(dungeon % get(j + (i - 1) * N_duplicates), .true.)
        @assertEqual(i, MGParticleStatePtr % getBroodId())
        @assertEqual(j, MGParticleStatePtr % getEnergyGroup())

      end do

    end do

    call dungeon % kill()

  end subroutine testSortingByBroodID_withDuplicates

end module particleDungeon_test