module microResponse_test

  use CENeutron_class,           only : CENeutron
  use dictionary_class,          only : dictionary
  use endfConstants
  use funit
  use materialMenu_mod,          only : init, kill
  use microResponse_class,       only : microResponse
  use numPrecision
  use testNeutronDatabase_class, only : testNeutronDatabase

  implicit none

@testCase
  type, extends(TestCase)     :: test_microResponse
    private
    type(CENeutron)           :: testNeutron
    type(microResponse)       :: response_absorption, response_capture, response_eScatter, &
                                 response_fission, response_heating, response_total
    type(testNeutronDatabase) :: xsData
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_microResponse


contains

  !!
  !! Sets up test_microResponse object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_microResponse), intent(inout) :: this
    type(dictionary)                         :: tempDict, dictMat1, dictMat2, dictMat3

    ! Allocate and initialise test nuclearData
    ! Cross-sections:          Total        eScatering   IeScatter Capture     Fission       nuFission    Heating
    call this % xsData % build(6.0_defReal, 3.0_defReal, ZERO,     2.0_defReal, 1.0_defReal, 1.5_defReal, 9.0_defReal)

    ! Set dictionaries to initialise material
    call dictMat1 % init(1)
    call dictMat2 % init(2)
    call dictMat3 % init(1)

    call dictMat3 % store('54135.03', 2.0_defReal)

    call dictMat2 % store('temp', 300.0_defReal)
    call dictMat2 % store('composition', dictMat3)

    call dictMat1 % store('Xenon', dictMat2)

    ! Initialise material
    call init(dictMat1)

    ! Set up responses
    ! Total
    call tempDict % init(3)
    call tempDict % store('type', 'microResponse')
    call tempDict % store('MT', N_TOTAL)
    call tempDict % store('material', 'Xenon')
    call this % response_total % init(tempDict)
    call tempDict % kill()

    ! Capture
    call tempDict % init(3)
    call tempDict % store('type', 'microResponse')
    call tempDict % store('MT', N_GAMMA)
    call tempDict % store('material', 'Xenon')
    call this % response_capture % init(tempDict)
    call tempDict % kill()

    ! Fission
    call tempDict % init(3)
    call tempDict % store('type', 'microResponse')
    call tempDict % store('MT', N_FISSION)
    call tempDict % store('material', 'Xenon')
    call this % response_fission % init(tempDict)
    call tempDict % kill()

    ! nuFission
    call tempDict % init(3)
    call tempDict % store('type', 'microResponse')
    call tempDict % store('MT', N_N_ELASTIC)
    call tempDict % store('material', 'Xenon')
    call this % response_eScatter % init(tempDict)
    call tempDict % kill()

    ! Absorbtion
    call tempDict % init(3)
    call tempDict % store('type', 'microResponse')
    call tempDict % store('MT', N_ABSORPTION)
    call tempDict % store('material', 'Xenon')
    call this % response_absorption % init(tempDict)
    call tempDict % kill()

    ! Heating.
    call tempDict % init(3)
    call tempDict % store('type', 'microResponse')
    call tempDict % store('MT', N_heating)
    call tempDict % store('material', 'Xenon')
    call this % response_heating % init(tempDict)
    call tempDict % kill()

    call this % testNeutron % init()

  end subroutine setUp

  !!
  !! Kills test_microResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_microResponse), intent(inout) :: this

    ! Kill and deallocate testTransportNuclearData
    call this % xsData % kill()
    call this % testNeutron % kill()
    call kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the filter
  !!
@Test
  subroutine testGettingResponse(this)
    class(test_microResponse), intent(inout) :: this
    real(defReal)                            :: result
    real(defReal), parameter                 :: TOL = 1.0E-9

    ! Test response values
    call this % response_total % get(this % testNeutron, result, this % xsData)
    @assertEqual(3.0_defReal, result, TOL)

    call this % response_capture % get(this % testNeutron, result, this % xsData)
    @assertEqual(1.0_defReal, result, TOL)

    call this % response_fission % get(this % testNeutron, result, this % xsData)
    @assertEqual(0.5_defReal, result, TOL)

    call this % response_eScatter % get(this % testNeutron, result, this % xsData)
    @assertEqual(1.5_defReal, result, TOL)

    call this % response_absorption % get(this % testNeutron, result, this % xsData)
    @assertEqual(1.5_defReal, result, TOL)

    call this % response_heating % get(this % testNeutron, result, this % xsData)
    @assertEqual(4.5_defReal, result, TOL)

  end subroutine testGettingResponse

end module microResponse_test