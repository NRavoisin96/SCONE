module densityResponse_test

  use densityResponse_class, only : densityResponse
  use dictionary_class,      only : dictionary
  use funit
  use numPrecision
  use particle_class,        only : particle, P_NEUTRON, P_PHOTON
  use universalVariables,    only : neutronMass, lightSpeed

  implicit none

@testCase
  type, extends(TestCase) :: test_densityResponse
    private
    type(densityResponse) :: response
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

  end subroutine setUp

  !!
  !! Kills test_densityResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_densityResponse), intent(inout) :: this

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
    type(particle)                             :: p
    real(defReal)                              :: ref, result
    real(defReal), parameter                   :: tol = 1.0e-9_defReal

    ! Test neutron density with different particle energies
    p % type = P_NEUTRON
    p % isMG = .false.
    p % E = ONE
    ref = ONE / lightSpeed / sqrt(TWO * p % E / neutronMass)
    call this % response % get(p, result)
    @assertEqual(ref, result, ref * tol)

    p % E = 1.6e-06_defReal
    ref = ONE / lightSpeed / sqrt(TWO * p % E / neutronMass)
    call this % response % get(p, result)
    @assertEqual(ref, result, ref * tol)

    ! Test photon density
    p % type = P_PHOTON
    ref = ONE / lightSpeed
    call this % response % get(p, result)
    @assertEqual(ref, result, ref * tol)

  end subroutine densityResponseing

end module densityResponse_test
