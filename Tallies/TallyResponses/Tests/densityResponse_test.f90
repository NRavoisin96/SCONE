module densityResponse_test

  use CENeutron_class,       only : CENeutron
  use CEPhoton_class,        only : CEPhoton
  use densityResponse_class, only : densityResponse
  use dictionary_class,      only : dictionary
  use funit
  use numPrecision
  use universalVariables,    only : neutronMass, lightSpeed

  implicit none

@testCase
  type, extends(TestCase) :: test_densityResponse
    private
    type(densityResponse) :: response
    type(CENeutron)       :: testNeutron
    type(CEPhoton)        :: testPhoton
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_densityResponse


contains

  !!
  !! Sets up test_densityResponse object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_densityResponse), intent(inout) :: this
    type(dictionary)                           :: tempDict

    call this % testNeutron % init()
    call this % testPhoton % init()

  end subroutine setUp

  !!
  !! Kills test_densityResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_densityResponse), intent(inout) :: this

    call this % testNeutron % kill()
    call this % testPhoton % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the response
  !!
@Test
  subroutine densityResponseing(this)
    class(test_densityResponse), intent(inout) :: this
    real(defReal)                              :: E, ref, result
    real(defReal), parameter                   :: tol = 1.0e-9_defReal

    ! Test neutron density with different particle energies
    E = ONE
    call this % testNeutron % setEnergy(E)
    ref = ONE / lightSpeed / sqrt(TWO * E / neutronMass)
    call this % response % get(this % testNeutron, result)
    @assertEqual(ref, result, ref * tol)

    E = 1.6e-06_defReal
    call this % testNeutron % setEnergy(E)
    ref = ONE / lightSpeed / sqrt(TWO * E / neutronMass)
    call this % response % get(this % testNeutron, result)
    @assertEqual(ref, result, ref * tol)

    ! Test photon density
    ref = ONE / lightSpeed
    call this % response % get(this % testPhoton, result)
    @assertEqual(ref, result, ref * tol)

  end subroutine densityResponseing

end module densityResponse_test
