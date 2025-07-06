module uniFissSitesField_test
  
  use dictionary_class,         only : dictionary
  use dictParser_func,          only : charToDict
  use funit
  use geometry_inter,           only : geometry
  use numPrecision
  use particle_class,           only : particle, particleState
  use uniFissSitesField_class,  only : uniFissSitesField
  use universalVariables,       only : ONE, ZERO

  implicit none


@testCase
  type, extends(TestCase) :: test_uniFissSitesField
    private
    type(uniFissSitesField) :: ufsField
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_uniFissSitesField

  !!
  !! Map definition
  !!
  character(*), parameter :: DICT_DEF = &
  " type spaceMap;  axis z;  grid unstruct; &
    &bins (0.0 20.0 40.0 60.0 80.0); "

contains

  !!
  !! Sets up test_weightWindows object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_uniFissSitesField), intent(inout) :: this
    type(dictionary)                             :: dict, dictMap
    class(geometry), pointer                     :: geom
    integer(shortInt)                            :: type

    call charToDict(dictMap, DICT_DEF)

    ! Initialise dictionaries
    call dict % init(3)

    ! Build material map definition
    call dict % store('type', 'uniFissSitesField')
    call dict % store('uniformVolMap', 1)
    call dict % store('map', dictMap)

    call this % ufsField % init(dict)
    call this % ufsField % estimateVol(geom, type)

  end subroutine setUp

  !!
  !! Kills test_weightWindows object
  !!
  subroutine tearDown(this)
    class(test_uniFissSitesField), intent(inout) :: this

    call this % ufsField % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test retrieving the ufs values
  !!
@Test
  subroutine testGetValue(this)
    class(test_uniFissSitesField), intent(inout) :: this
    type(particle)                               :: p
    type(particleState)                          :: state
    real(defReal), dimension(3)                  :: bins, EXPECTED_BINS

    ! Test case in the map
    call p % coords % setPosition([0.5_defReal, 7.0_defReal, 50.0_defReal], 1)

    bins = this % ufsField % at(p)
    EXPECTED_BINS = [0.25_defReal, 0.25_defReal, ZERO]
    @assertEqual(EXPECTED_BINS, bins, tolerance=1.0e-6_defReal)

    ! Test case outside the map
    call p % coords % setPosition([0.5_defReal, 7.0_defReal, 100.0_defReal], 1)

    bins = this % ufsField % at(p)
    EXPECTED_BINS = [ONE, ONE, ONE]
    @assertEqual(EXPECTED_BINS,bins)

    ! Modify the map by storing fission sites
    state % r   = [0.5_defReal, 7.0_defReal, 12.0_defReal]
    state % wgt = 0.2_defReal

    call this % ufsField % storeFS(state)

    state % r   = [0.5_defReal, 7.0_defReal, 23.2_defReal]
    state % wgt = 0.8_defReal
    call this % ufsField % storeFS(state)

    call this % ufsField % updateMap()

    ! Test case in the updated map
    call p % coords % setPosition([0.5_defReal, 7.0_defReal, 18.1_defReal], 1)

    bins = this % ufsField % at(p)
    EXPECTED_BINS = [0.25_defReal, 0.06666666667_defReal, ZERO]
    @assertEqual(EXPECTED_BINS, bins, tolerance=1.0e-6_defReal)

  end subroutine testGetValue

end module uniFissSitesField_test