module mgXsClerk_test

  use CEParticleState_class,      only : castCEParticleStatePtr, CEParticleState
  use dictionary_class,           only : dictionary
  use endfConstants
  use funit
  use genericProcedures,          only : numToChar
  use mgXsClerk_class,            only : mgXsClerk
  use numPrecision
  use outputFile_class,           only : outputFile
  use scoreMemory_class,          only : scoreMemory
  use testCEParticle_class,       only : testCEParticle
  use testNeutronDatabase_class,  only : testNeutronDatabase
  use transportObjectState_class, only : transportObjectState

  implicit none

  @testCase
    type, extends(TestCase) :: test_mgXsClerk
      private
      type(mgXsClerk)           :: testClerk1, testClerk2
      type(testNeutronDatabase) :: nucData
      type(testCEParticle)      :: testParticle
    contains
      procedure :: setUp
      procedure :: tearDown
  end type test_mgXsClerk

contains
@Before
  !!
  !! Sets up test_mgXsClerk object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_mgXsClerk), intent(inout) :: this
    character(nameLen)                   :: temp
    type(dictionary)                     :: tempDict, energyDict, spaceDict

    ! Build energy map
    call energyDict % init(3)
    call energyDict % store('type', 'energyMap')
    call energyDict % store('grid', 'unstruct')
    call energyDict % store('bins', [1.0E-03_defReal, ONE, 10.0_defReal])

    ! Build material map
    call spaceDict % init(2)
    call spaceDict % store('type', 'testMap')
    call spaceDict % store('maxIdx', 2)

    ! Define first clerk, with high order scattering and spatial map
    call tempDict % init(2)
    call tempDict % store('energyMap', energyDict)
    call tempDict % store('spaceMap', spaceDict)

    temp = 'MGxs1'
    call this % testClerk1 % init(tempDict, temp)
    call spaceDict % kill()
    call energyDict % kill()
    call tempDict % kill()

    ! Build energy map
    call energyDict % init(3)
    call energyDict % store('type', 'energyMap')
    call energyDict % store('grid', 'unstruct')
    call energyDict % store('bins', [1.0E-11_defReal, 0.6_defReal, 1.2_defReal, 20.0_defReal])

    ! Define second clerk, without spatial map and high order scattering
    call tempDict % init(2)
    call tempDict % store('energyMap', energyDict)
    call tempDict % store('PN', 0)

    temp = 'MGxs2'
    call this % testClerk2 % init(tempDict, temp)
    call energyDict % kill()
    call tempDict % kill()

    ! Build test neutronDatabase
    call this % nucData % build(ONE, captureXS = 2.0_defReal, fissionXS = 1.5_defReal, nuFissionXS = 3.0_defReal)

    ! Build test particle.
    call this % testParticle % init()

  end subroutine setUp

@After
  !!
  !! Kills test case object
  !!
  subroutine tearDown(this)
    class(test_mgXsClerk), intent(inout) :: this

    call this % testClerk1 % kill()
    call this % testClerk2 % kill()
    call this % nucData % kill()
    call this % testParticle % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Scoring test for clerk 1
  !!
@Test
  subroutine testScoring_clerk1(this)
    class(test_mgXsClerk), intent(inout)        :: this
    type(CEParticleState), pointer              :: currentStatePtr, preCollisionStatePtr
    real(defReal), dimension(:, :), allocatable :: capt, chi, fiss, nu, prod, P0, P1, P2, P3, &
                                                   P4, P5, P6, P7, transFL, transOS
    type(scoreMemory)                           :: mem
    type(outputFile)                            :: out
    real(defReal), parameter                    :: TOL = 1.0E-9

    ! Configure memory
    call mem % init(1000_longInt, 1)
    call this % testClerk1 % setMemAddress(1_longInt)

    ! Report fission particles
    call this % testParticle % setWeight(0.5_defReal)
    call this % testParticle % setEnergy(0.3_defReal)
    call this % testParticle % setMaterialIdx(2)

    ! Scoring
    currentStatePtr => castCEParticleStatePtr(this % testParticle % updateAndGetCurrentStatePtr(), .true.)
    call this % testClerk1 % reportSpawn(N_FISSION, this % testParticle, currentStatePtr, this % nucData, mem)

    call this % testParticle % setEnergy(3.0_defReal)

    ! Scoring
    call this % testClerk1 % reportSpawn(N_FISSION, this % testParticle, currentStatePtr, this % nucData, mem)

    ! Scoring
    call this % testClerk1 % reportInColl(this % testParticle, .false., this % nucData, mem)

    preCollisionStatePtr => castCEParticleStatePtr(this % testParticle % getPreCollisionStatePtr(), .true.)
    call preCollisionStatePtr % setWeight(0.2_defReal)
    call preCollisionStatePtr % setEnergy(3.0_defReal)
    call preCollisionStatePtr % setMaterialIdx(2)
    call this % testParticle % setEnergy(0.1_defReal)

    call this % testClerk1 % reportOutColl(this % testParticle, N_2N, 0.75_defReal, this % nucData, mem)
    call mem % closeCycle(ONE)

    ! Process and get results
    call this % testClerk1 % processRes(mem, capt, fiss, transFL, transOS, nu, chi, P0, P1, prod)
    call this % testClerk1 % processPN(mem, P2, P3, P4, P5, P6, P7)

    ! Verify results of scoring
    @assertEqual([ZERO, ZERO, TWO, ZERO], capt(1, :), TOL, 'Capture XS')
    @assertEqual([ZERO, ZERO, 1.5_defReal, ZERO], fiss(1, :), TOL, 'Fission XS')
    @assertEqual([ZERO, ZERO, TWO, ZERO], nu(1, :), TOL, 'NuFission XS')
    @assertEqual([ZERO, ZERO, HALF, HALF], chi(1, :), TOL, 'Chi')
    @assertEqual([ZERO, ZERO, 4.0_defReal, ZERO], transOS(1, :), TOL, 'Transport XS O.S.')
    @assertEqual([ZERO, ZERO, 5.5_defReal, ZERO], transFL(1, :), TOL, 'Transport XS F.L.')
    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, TWO, ZERO, ZERO], P0(1, :), TOL, 'P0')
    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, 1.5_defReal, ZERO, ZERO], P1(1, :), TOL, 'P1')
    @assertEqual([ONE, ONE, ONE, ONE, ONE, TWO, ONE, ONE], prod(1, :), TOL, 'prod')

    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, 0.6875_defReal, ZERO, ZERO], P2(1, :), TOL, 'P2')
    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, -0.140625_defReal, ZERO, ZERO], P3(1, :), TOL, 'P3')
    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, -0.7001953125_defReal, ZERO, ZERO], P4(1, :), TOL, 'P4')
    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, -0.8327636719_defReal, ZERO, ZERO], P5(1, :), TOL, 'P5')
    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, -0.5615539551_defReal, ZERO, ZERO], P6(1, :), TOL, 'P6')
    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, -0.0683670044_defReal, ZERO, ZERO], P7(1, :), TOL, 'P7')

    ! Test getting size
    @assertEqual(100, this % testClerk1 % getSize(), 'Test getSize(): ')

    ! Test correctness of output calls
    call out % init('dummyPrinter', fatalErrors = .false.)
    call this % testClerk1 % print(out, mem)
    @assertTrue(out % isValid(), 'Test print(): ')

  end subroutine testScoring_clerk1

  !!
  !! Scoring test for clerk 2
  !!
@Test
  subroutine testScoring_clerk2(this)
    class(test_mgXsClerk), intent(inout)        :: this
    real(defReal), dimension(:, :), allocatable :: capt, chi, fiss, nu, prod, P0, P1, transFL, transOS
    type(CEParticleState), pointer              :: currentStatePtr, preCollisionStatePtr
    type(scoreMemory)                           :: mem
    type(outputFile)                            :: out
    real(defReal), parameter                    :: TOL = 1.0E-9

    ! Configure memory
    call mem % init(1000_longInt, 1)
    call this % testClerk2 % setMemAddress(1_longInt)

    ! Report fission particles
    call this % testParticle % setWeight(0.5_defReal)
    call this % testParticle % setEnergy(3.0_defReal)
    call this % testParticle % setMaterialIdx(2)

    ! Scoring
    currentStatePtr => castCEParticleStatePtr(this % testParticle % updateAndGetCurrentStatePtr(), .true.)
    call this % testClerk2 % reportSpawn(N_FISSION, this % testParticle, currentStatePtr, this % nucData, mem)

    call this % testParticle % setEnergy(0.3_defReal)

    ! Scoring
    call this % testClerk2 % reportSpawn(N_FISSION, this % testParticle, currentStatePtr, this % nucData, mem)

    ! Scoring
    call this % testClerk2 % reportInColl(this % testParticle, .false., this % nucData, mem)

    preCollisionStatePtr => castCEParticleStatePtr(this % testParticle % getPreCollisionStatePtr(), .true.)
    call preCollisionStatePtr % setWeight(0.2_defReal)
    call preCollisionStatePtr % setEnergy(0.3_defReal)
    call this % testParticle % setEnergy(1.1_defReal)

    call this % testClerk2 % reportOutColl(this % testParticle, N_2N, 0.75_defReal, this % nucData, mem)
    call mem % closeCycle(ONE)

    ! Process and get results
    call this % testClerk2 % processRes(mem, capt, fiss, transFL, transOS, nu, chi, P0, P1, prod)

    ! Verify results of scoring
    @assertEqual([ZERO, ZERO, TWO], capt(1, :), TOL, 'Capture XS')
    @assertEqual([ZERO, ZERO, 1.5_defReal], fiss(1, :), TOL, 'Fission XS')
    @assertEqual([ZERO, ZERO, TWO], nu(1, :), TOL, 'NuFission XS')
    @assertEqual([HALF, ZERO, HALF], chi(1, :), TOL, 'Chi')
    @assertEqual([ZERO, ZERO, 4.0_defReal], transOS(1, :), TOL, 'Transport XS O.S.')
    @assertEqual([ZERO, ZERO, 5.5_defReal], transFL(1, :), TOL, 'Transport XS F.L.')
    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, TWO, ZERO], P0(1, :), TOL, 'P0')
    @assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, 1.5_defReal, ZERO], P1(1, :), TOL, 'P1')
    @assertEqual([ONE, ONE, ONE, ONE, ONE, ONE, ONE, TWO, ONE], prod(1, :), TOL, 'prod')

    ! Test getting size
    @assertEqual(48, this % testClerk2 % getSize(), 'Test getSize(): ')

    ! Test correctness of output calls
    call out % init('dummyPrinter', fatalErrors = .false.)
    call this % testClerk2 % print(out, mem)
    @assertTrue(out % isValid(), 'Test print(): ')

  end subroutine testScoring_clerk2

end module mgXsClerk_test