module uniFissSitesField_test
  
  use dictionary_class,           only : dictionary
  use dictParser_func,            only : charToDict
  use funit
  use geometry_inter,             only : geometry
  use numPrecision
  use testTransportObject_class,  only : testTransportObject
  use transportObjectState_class, only : transportObjectState
  use uniFissSitesField_class,    only : uniFissSitesField
  use universalVariables,         only : ONE, ZERO

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

  ! Variables.
  type(testTransportObject) :: testObject

contains
@Before
  !!
  !! Sets up test_weightWindows object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_uniFissSitesField), intent(inout) :: this
    class(geometry), pointer                     :: geom
    integer(shortInt)                            :: type
    type(dictionary)                             :: dict, dictMap

    call charToDict(dictMap, DICT_DEF)

    ! Initialise dictionaries
    call dict % init(3)

    ! Build material map definition
    call dict % store('type', 'uniFissSitesField')
    call dict % store('uniformVolMap', 1)
    call dict % store('map', dictMap)
    call this % ufsField % init(dict)
    call this % ufsField % estimateVol(geom, type)

    ! Initialise test transport object.
    call testObject % init()

  end subroutine setUp

@After
  !!
  !! Kills test_weightWindows object
  !!
  subroutine tearDown(this)
    class(test_uniFissSitesField), intent(inout) :: this

    call this % ufsField % kill()
    call testObject % kill()

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
    real(defReal), dimension(3)                  :: bins, EXPECTED_BINS
    type(transportObjectState)                   :: testState
    real(defReal), parameter                     :: TOL = 1.0e-6_defReal

    ! Test case in the map
    call testObject % setGlobalPosition([0.5_defReal, 7.0_defReal, 50.0_defReal])
    bins = this % ufsField % at(testObject)
    EXPECTED_BINS = [0.25_defReal, 0.25_defReal, ZERO]
    @assertEqual(EXPECTED_BINS, bins, tolerance = TOL)

    ! Test case outside the map
    call testObject % setGlobalPosition([0.5_defReal, 7.0_defReal, 100.0_defReal])
    bins = this % ufsField % at(testObject)
    EXPECTED_BINS = [ONE, ONE, ONE]
    @assertEqual(EXPECTED_BINS, bins, tolerance = TOL)

    ! Modify the map by storing fission sites
    call testState % setGlobalPosition([0.5_defReal, 7.0_defReal, 12.0_defReal])
    call testState % setWeight(0.2_defReal)
    call this % ufsField % storeFS(testState)

    call testState % setGlobalPosition([0.5_defReal, 7.0_defReal, 23.2_defReal])
    call testState % setWeight(0.8_defReal)
    call this % ufsField % storeFS(testState)

    call this % ufsField % updateMap()

    ! Test case in the updated map
    call testObject % setGlobalPosition([0.5_defReal, 7.0_defReal, 18.1_defReal])
    bins = this % ufsField % at(testObject)
    EXPECTED_BINS = [0.25_defReal, 0.06666666667_defReal, ZERO]
    @assertEqual(EXPECTED_BINS, bins, tolerance = TOL)

  end subroutine testGetValue

end module uniFissSitesField_test