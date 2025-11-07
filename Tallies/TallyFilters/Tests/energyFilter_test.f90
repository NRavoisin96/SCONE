module energyFilter_test

  use CEParticleState_class, only : CEParticleState
  use dictionary_class,      only : dictionary
  use energyFilter_class,    only : energyFilter
  use funit
  use MGParticleState_class, only : MGParticleState
  use numPrecision

  implicit none

@testCase
  type, extends(TestCase) :: test_energyFilter
    private
    type(energyFilter) :: cemgFilter, filter, mgFilter
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_energyFilter

  !! Parameters
  integer(shortInt), parameter :: G_LOW = 7, G_TOP = 3
  real(defReal), parameter     :: E_Delta = 0.999_defReal, E_MAX = 1.34_defReal, E_MIN = 1.27E-6_defReal ! E_Delta Must be < 1 for correct tests !
contains
@Before
  !!
  !! Sets up test_energyFilter object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_energyFilter), intent(inout) :: this
    type(dictionary)                        :: tempDict

    ! Build CE filter
    call tempDict % init(4)
    call tempDict % store('Emin', E_MIN)
    call tempDict % store('Emax', E_MAX)
    call this % filter % init(tempDict)

    ! Build MG-CE filter
    call tempDict % store('Gtop', G_TOP)
    call tempDict % store('Glow', G_LOW)
    call this % cemgFilter % init(tempDict)

    ! Build MG filter
    call tempDict % kill()
    call tempDict % init(2)
    call tempDict % store('Gtop', G_TOP)
    call tempDict % store('Glow', G_LOW)
    call this % mgFilter % init(tempDict)

  end subroutine setUp

@After
  !!
  !! Kills test_energyFilter object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_energyFilter), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the CE filter
  !!
@Test
  subroutine testCEFilter(this)
    class(test_energyFilter), intent(inout) :: this
    logical(defBool)                        :: filterRes
    real(defReal)                           :: testE
    type(CEParticleState)                   :: state

    ! Below Emin -> FALSE
    testE = E_MIN * E_Delta
    call state % setEnergy(testE)
    filterRes = this % filter % isPass(state)
    @assertFalse(filterRes)

    ! Emin -> TRUE
    testE = E_MIN
    call state % setEnergy(testE)
    filterRes = this % filter % isPass(state)
    @assertTrue(filterRes)

    ! Below E max -> TRUE
    testE = E_MAX * E_Delta
    call state % setEnergy(testE)
    filterRes = this % filter % isPass(state)
    @assertTrue(filterRes)

    ! E max -> TRUE
    testE = E_MAX
    call state % setEnergy(testE)
    filterRes = this % filter % isPass(state)
    @assertTrue(filterRes)

    ! Above Emax -> FALSE
    testE = E_MAX / E_Delta
    call state % setEnergy(testE)
    filterRes = this % filter % isPass(state)
    @assertFalse(filterRes)

  end subroutine testCEFilter

  !!
  !! Test correct behaviour of the CEMG filter
  !!
@Test
  subroutine testMGFilter(this)
    class(test_energyFilter), intent(inout) :: this
    integer(shortInt)                       :: testG
    logical(defBool)                        :: filterRes
    type(MGParticleState)                   :: state

    ! Below G_TOP -> FALE
    testG = G_TOP - 1
    call state % setEnergyGroup(testG)
    filterRes = this % cemgFilter % isPass(state)
    @assertFalse(filterRes)

    ! Exactly G_TOP -> TRUE
    testG = G_TOP
    call state % setEnergyGroup(testG)
    filterRes = this % cemgFilter % isPass(state)
    @assertTrue(filterRes)

    ! Above G_TOP -> TRUE
    testG = G_TOP + 1
    call state % setEnergyGroup(testG)
    filterRes = this % cemgFilter % isPass(state)
    @assertTrue(filterRes)

    ! Above G_LOW -> FALSE
    testG = G_LOW + 1
    call state % setEnergyGroup(testG)
    filterRes = this % cemgFilter % isPass(state)
    @assertFalse(filterRes)

    ! Exactly G_LOW -> TRUE
    testG = G_LOW
    call state % setEnergyGroup(testG)
    filterRes = this % cemgFilter % isPass(state)
    @assertTRUE(filterRes)

  end subroutine testMGFilter

  !!
  !! Test correct behaviour of the CEMG filter
  !!
@Test
  subroutine testCEMGFilter(this)
    class(test_energyFilter), intent(inout) :: this
    integer(shortInt)                       :: testG
    logical(defBool)                        :: filterRes
    real(defReal)                           :: testE
    type(CEParticleState)                   :: CEState
    type(MGParticleState)                   :: MGState

    ! Below Emin -> FALSE
    testE = E_MIN * E_Delta
    call CEState % setEnergy(testE)
    filterRes = this % cemgFilter % isPass(CEState)
    @assertFalse(filterRes)

    ! Emin -> TRUE
    testE = E_MIN
    call CEState % setEnergy(testE)
    filterRes = this % cemgFilter % isPass(CEState)
    @assertTrue(filterRes)

    ! Below E max -> TRUE
    testE = E_MAX * E_Delta
    call CEState % setEnergy(testE)
    filterRes = this % cemgFilter % isPass(CEState)
    @assertTrue(filterRes)

    ! E max -> TRUE
    testE = E_MAX
    call CEState % setEnergy(testE)
    filterRes = this % cemgFilter % isPass(CEState)
    @assertTrue(filterRes)

    ! Above Emax -> FALSE
    testE = E_MAX / E_Delta
    call CEState % setEnergy(testE)
    filterRes = this % cemgFilter % isPass(CEState)
    @assertFalse(filterRes)

    ! Below G_TOP -> FALE
    testG = G_TOP - 1
    call MGState % setEnergyGroup(testG)
    filterRes = this % cemgFilter % isPass(MGState)
    @assertFalse(filterRes)

    ! Exactly G_TOP -> TRUE
    testG = G_TOP
    call MGState % setEnergyGroup(testG)
    filterRes = this % cemgFilter % isPass(MGState)
    @assertTrue(filterRes)

    ! Above G_TOP -> TRUE
    testG = G_TOP + 1
    call MGState % setEnergyGroup(testG)
    filterRes = this % cemgFilter % isPass(MGState)
    @assertTrue(filterRes)

    ! Above G_LOW -> FALSE
    testG = G_LOW + 1
    call MGState % setEnergyGroup(testG)
    filterRes = this % cemgFilter % isPass(MGState)
    @assertFalse(filterRes)

    ! Exactly G_LOW -> TRUE
    testG = G_LOW
    call MGState % setEnergyGroup(testG)
    filterRes = this % cemgFilter % isPass(MGState)
    @assertTrue(filterRes)

  end subroutine testCEMGFilter

end module energyFilter_test