module weightResponse_test

  use dictionary_class,           only : dictionary
  use endfConstants
  use funit
  use numPrecision
  use testNeutronDatabase_class,  only : testNeutronDatabase
  use testPhysicalParticle_class, only : testPhysicalParticle
  use weightResponse_class,       only : weightResponse

  implicit none

@testCase
  type, extends(TestCase) :: test_weightResponse
    private
    type(testNeutronDatabase)  :: xsData
    type(testPhysicalParticle) :: testParticle
    type(weightResponse)       :: response_weight_m0, response_weight_m2
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_weightResponse


contains

  !!
  !! Sets up test_macroResponse object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_weightResponse), intent(inout) :: this
    type(dictionary)                          :: tempDict

    ! Cross-sections:         Total
    call this % xsData % build(4.0_defReal)

    ! Set up weight response
    call tempDict % init(1)
    call tempDict % store('moment', 0)
    call this % response_weight_m0 % init(tempDict)
    call tempDict % kill()

    call tempDict % init(1)
    call tempDict % store('moment', 2)
    call this % response_weight_m2 % init(tempDict)
    call tempDict % kill()

    call this % testParticle % init()

  end subroutine setUp

  !!
  !! Kills test_weightResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_weightResponse), intent(inout) :: this

    ! Kill and deallocate testTransportNuclearData
    call this % xsData % kill()
    call this % testParticle % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the filter
  !!
@Test
  subroutine testGettingResponse(this)
    class(test_weightResponse), intent(inout) :: this
    real(defReal)                             :: result
    real(defReal), parameter                  :: TOL = 1.0E-9

    call this % testParticle % setWeight(TWO)

    ! Test response values
    call this % response_weight_m0 % get(this % testParticle, result, this % xsData)
    @assertEqual(TWO, result, TOL)

    call this % response_weight_m2 % get(this % testParticle, result, this % xsData)
    @assertEqual(8.0_defReal, result, TOL)

  end subroutine testGettingResponse

end module weightResponse_test