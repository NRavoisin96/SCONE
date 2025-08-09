module macroResponse_test

  use dictionary_class,          only : dictionary
  use endfConstants
  use funit
  use macroResponse_class,       only : macroResponse
  use numPrecision
  use particle_class,            only : particle, P_NEUTRON
  use testNeutronDatabase_class, only : testNeutronDatabase

  implicit none

@testCase
  type, extends(TestCase)     :: test_macroResponse
    private
    type(macroResponse)       :: response_absorption
    type(macroResponse)       :: response_capture
    type(macroResponse)       :: response_fission
    type(macroResponse)       :: response_heating
    type(macroResponse)       :: response_nuFission
    type(macroResponse)       :: response_total
    type(testNeutronDatabase) :: xsData
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_macroResponse


contains

  !!
  !! Sets up test_macroResponse object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_macroResponse), intent(inout) :: this
    type(dictionary)                         :: tempDict

    ! Allocate and initialise test nuclearData
    ! Cross-sections:          Total        eScatering   IeScatter Capture      Fission      nuFission    Heating
    call this % xsData % build(6.0_defReal, 3.0_defReal, ZERO,     2.0_defReal, 1.0_defReal, 1.5_defReal, 9.0_defReal)

    ! Set up responses
    ! Total
    call tempDict % init(2)
    call tempDict % store('type', 'macroResponse')
    call tempDict % store('MT', macroTotal)
    call this % response_total % init(tempDict)
    call tempDict % kill()

    ! Capture
    call tempDict % init(2)
    call tempDict % store('type', 'macroResponse')
    call tempDict % store('MT', macroCapture)
    call this % response_capture % init(tempDict)
    call tempDict % kill()

    ! Fission
    call tempDict % init(2)
    call tempDict % store('type', 'macroResponse')
    call tempDict % store('MT', macroFission)
    call this % response_fission % init(tempDict)
    call tempDict % kill()

    ! nuFission
    call tempDict % init(2)
    call tempDict % store('type', 'macroResponse')
    call tempDict % store('MT', macroNuFission)
    call this % response_nuFission % init(tempDict)
    call tempDict % kill()

    ! Absorbtion
    call tempDict % init(2)
    call tempDict % store('type', 'macroResponse')
    call tempDict % store('MT', macroAbsorption)
    call this % response_absorption % init(tempDict)
    call tempDict % kill()

    ! Heating
    call tempDict % init(2)
    call tempDict % store('type', 'macroResponse')
    call tempDict % store('MT', macroHeating)
    call this % response_heating % init(tempDict)
    call tempDict % kill()

  end subroutine setUp

  !!
  !! Kills test_macroResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_macroResponse), intent(inout) :: this

    ! Kill and deallocate testTransportNuclearData
    call this % xsData % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the filter
  !!
@Test
  subroutine testGettingResponse(this)
    class(test_macroResponse), intent(inout) :: this
    type(particle)                           :: p
    real(defReal)                            :: result
    real(defReal), parameter                 :: tol = 1.0E-9

    p % type = P_NEUTRON

    ! Test response values
    call this % response_total % get(p, result, this % xsData)
    @assertEqual(6.0_defReal, result, tol)

    call this % response_capture % get(p, result, this % xsData)
    @assertEqual(2.0_defReal, result, tol)

    call this % response_fission % get(p, result, this % xsData)
    @assertEqual(1.0_defReal, result, tol)

    call this % response_nuFission % get(p, result, this % xsData)
    @assertEqual(1.5_defReal, result, tol)

    call this % response_absorption % get(p, result, this % xsData)
    @assertEqual(3.0_defReal, result, tol)

    call this % response_heating % get(p, result, this % xsData)
    @assertEqual(9.0_defReal, result, tol)

  end subroutine testGettingResponse

end module macroResponse_test
