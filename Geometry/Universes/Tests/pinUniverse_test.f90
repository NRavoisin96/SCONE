module pinUniverse_test

  use cellShelf_class,    only : cellShelf
  use charMap_class,      only : charMap
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use funit
  use meshShelf_class,    only : meshShelf
  use numPrecision
  use pinUniverse_class,  only : pinUniverse, MOVING_IN, MOVING_OUT
  use publicObjects,      only : coordData, newCoordData
  use surfaceShelf_class, only : surfaceShelf
  use universalVariables, only : HALF, INF, ONE, SURF_TOL, ZERO
  
  implicit none

  ! Parameters
  character(*), parameter :: UNI_DEF = &
  "id 7; type pinUniverse; origin (0.0 0.0 0.0); rotation (0.0 0.0 0.0); &
  &radii (2.5 1.5 0.0); fills (u<7> u<14> void);"

  ! Variables
  type(surfaceShelf) :: surfs
  type(cellShelf)    :: cells
  type(meshShelf)    :: meshes
  type(charMap)      :: mats
  type(pinUniverse)  :: uni


contains

  !!
  !! Set-up test environment
  !!
@Before
  subroutine setup()
    character(nameLen)                           :: name
    integer(shortInt), dimension(:), allocatable :: fills
    type(dictionary)                             :: dict

    ! Load void material
    name = 'void'
    call mats % add(name, 13)

    ! Build universe
    call charToDict(dict, UNI_DEF)
    call uni % init(dict, mats, fills, cells, surfs, meshes)

    ! Set index
    call uni % setIdx(3)

    ! Verify fill array
    @assertEqual([-14, -7, 13], fills)


  end subroutine setup

  !!
  !! Clean after test
  !!
@After
  subroutine clean()

    call surfs % kill()
    call cells % kill()
    call mats % kill()
    call uni % kill()

  end subroutine clean

  !!
  !! Test miscellaneous functionality
  !!
@Test
  subroutine test_misc()

    ! Get id
    @assertEqual(7, uni % id())

    ! Set ID
    call uni % setId(7)
    @assertEqual(7, uni % id())

  end subroutine test_misc

  !!
  !! Test entering a universe
  !!
@Test
  subroutine test_enter()
    type(coordData)             :: data
    real(defReal), dimension(3) :: r_ref, u_ref, r, u
    real(defReal), parameter    :: TOL = 1.0E-7_defReal

    ! ** Enter into local cell 1
    r = [ZERO, ONE, ZERO]
    u = [ZERO, ZERO, ONE]
    data = newCoordData(r, u)
    call uni % enter(data)

    ! Verify location
    r_ref = r
    u_ref = u
    @assertEqual(r_ref, data % r, TOL)
    @assertEqual(u_ref, data % u, TOL)
    @assertEqual(3, data % universeIdx)
    @assertEqual(1, data % localId)
    @assertEqual(0, data % cellIdx)

    ! ** Enter into local cell 2
    r = [2.3_defReal, ZERO, -980.0_defReal]
    u = [ZERO, ZERO, ONE]
    data = newCoordData(r, u)
    call uni % enter(data)

    ! Verify location
    r_ref = r
    u_ref = u
    @assertEqual(r_ref, data % r, TOL)
    @assertEqual(u_ref, data % u, TOL)
    @assertEqual(3, data % universeIdx)
    @assertEqual(2, data % localId)
    @assertEqual(0, data % cellIdx)

    ! ** Enter into local cell 3
    r = [2.6_defReal, ZERO, -980.0_defReal]
    u = [ZERO, ZERO, ONE]
    data = newCoordData(r, u)
    call uni % enter(data)

    ! Verify location
    r_ref = r
    u_ref = u
    @assertEqual(r_ref, data % r, TOL)
    @assertEqual(u_ref, data % u, TOL)
    @assertEqual(3, data % universeIdx)
    @assertEqual(3, data % localId)
    @assertEqual(0, data % cellIdx)

    ! VERIFY THAT ROTATION IS NOT SET (all angles were 0.0)
    @assertFalse(data % isRotated)

  end subroutine test_enter

  !!
  !! Test distance calculation
  !!
@Test
  subroutine test_distance()
    type(coordData)          :: data
    real(defReal)            :: ref
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! ** In local cell 1 distance to boundary
    data = newCoordData([ONE, ZERO, ZERO], [ONE, ZERO, ZERO], localId = 1, universeIdx = 3)
    call uni % distance(data)
    ref = HALF
    @assertEqual(ref, data % d, ref * tol)
    @assertEqual(MOVING_OUT, data % surfaceIdx)

    ! ** In outermost cell moving away
    data = newCoordData([2.0_defReal, 1.6_defReal, ZERO], [ONE, ZERO, ZERO], localId = 3, universeIdx = 3)
    call uni % distance(data)
    @assertEqual(INF, data % d)
    ! Surface momento is undefined -> No crossing

    ! In ordinary cell in-between
    data = newCoordData([ZERO, 1.6_defReal, ZERO], [ZERO, -ONE, ZERO], localId = 2, universeIdx = 3)
    call uni % distance(data)
    ref = 0.1_defReal
    @assertEqual(ref, data % d, ref * tol)
    @assertEqual(MOVING_IN, data % surfaceIdx)

  end subroutine test_distance

  !!
  !! Test cell-to cell crossing
  !!
@Test
  subroutine test_cross()
    type(coordData) :: data
    real(defReal)   :: eps

    ! Cross from cell 1 to cell 2
    eps = HALF * SURF_TOL
    data = newCoordData([ZERO, 1.5_defReal - eps, ZERO], [ZERO, ONE, ZERO], localId = 1, surfaceIdx = MOVING_OUT, universeIdx = 8)
    call uni % cross(data)
    @assertEqual(2, data % localId)

    ! Cross form cell 2 to cell 1
    data = newCoordData([ZERO, 1.5_defReal + eps, ZERO], [ZERO, -ONE, ZERO], localId = 2, surfaceIdx = MOVING_IN, universeIdx = 8)
    call uni % cross(data)
    @assertEqual(1, data % localId)

  end subroutine test_cross

  !!
  !! Test surface transitions
  !!
  !! Check that there is no problem with distance calculations
  !! if particle is placed very close to an annulus surface (within SURF_TOL)
  !!
@Test
  subroutine test_edgeCases()
    type(coordData)             :: data
    real(defReal)               :: eps
    real(defReal), parameter    :: TOL = 1.0E-7_defReal

    ! At boundary between cell 1 and 2
    eps = HALF * SURF_TOL
    data = newCoordData([ZERO, 1.5_defReal - eps, ZERO], [ONE, -0.00001_defReal, ZERO], universeIdx = 8)

    ! Should find particle in cell 1
    ! And return very small distance -> MOVING OUT
    call uni % findCell(data)
    @assertEqual(1, data % localId)
    call uni % distance(data)

    @assertEqual(ZERO, data % d, 1.0E-3_defReal)
    @assertEqual(MOVING_OUT, data % surfaceIdx)

  end subroutine test_edgeCases

end module pinUniverse_test