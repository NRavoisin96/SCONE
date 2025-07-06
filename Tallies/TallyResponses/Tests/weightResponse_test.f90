module weightResponse_test

  use dictionary_class,          only : dictionary
  use endfConstants
  use funit
  use numPrecision
  use particle_class,            only : particle, P_NEUTRON
  use testNeutronDatabase_class, only : testNeutronDatabase
  use weightResponse_class,      only : weightResponse

  implicit none

@testCase
  type, extends(TestCase) :: test_weightResponse
    private
    type(weightResponse)      :: response_weight_m0
    type(weightResponse)      :: response_weight_m2
    type(testNeutronDatabase) :: xsData
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

  end subroutine setUp

  !!
  !! Kills test_weightResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_weightResponse), intent(inout) :: this

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
    class(test_weightResponse), intent(inout) :: this
    type(particle)                            :: p
    real(defReal)                             :: result
    real(defReal), parameter                  :: tol = 1.0E-9

    p % type = P_NEUTRON
    p % w = 2.0_defReal

    ! Test response values
    call this % response_weight_m0 % get(p, result, this % xsData)
    @assertEqual(2.0_defReal, result, tol)

    call this % response_weight_m2 % get(p, result, this % xsData)
    @assertEqual(8.0_defReal, result, tol)

  end subroutine testGettingResponse

end module weightResponse_test
