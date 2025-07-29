module cellUniverse_test

  use cellShelf_class,    only : cellShelf
  use cellUniverse_class, only : cellUniverse
  use charMap_class,      only : charMap
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use funit
  use genericProcedures
  use meshShelf_class,    only : meshShelf
  use numPrecision
  use publicObjects,      only : coordData, newCoordData
  use universalVariables, only : ONE, UNDEF_MAT, ZERO
  use surfaceShelf_class, only : surfaceShelf

  implicit none

  ! Parameters
  character(*), parameter :: SURF_DEF = &
  " surf1 { id 1; type sphere; origin (0.0 0.0 0.0); radius 2;}&
  & surf2 { id 2; type sphere; origin (4.0 0.0 0.0); radius 1;}"

  character(*), parameter :: CELL_DEF = &
  " cell1 {id 1; type simpleCell; surfaces (-1); filltype uni; universe 3;} &
  & cell2 {id 2; type simpleCell; surfaces (1 2); filltype uni; universe 4;}"

  !
  ! Note that rotation is such that following axis transformation applies:
  !   x -> z
  !   y -> -y
  !   z -> x
  !
  character(*), parameter :: UNI_DEF = &
  "id 1; type cellUniverse; origin (2.0 0.0 0.0); rotation (90.0 90.0 90.0); cells (1 2);"

  ! Variables
  type(surfaceShelf) :: surfs
  type(meshShelf)    :: meshes
  type(cellShelf)    :: cells
  type(charMap)      :: mats
  type(cellUniverse) :: uni

contains

  !!
  !! Setup environment
  !!
@Before
  subroutine setUp()
    integer(shortInt), dimension(:), allocatable :: fills
    type(dictionary) :: dict

    ! Build surfaces.
    call charToDict(dict, SURF_DEF)
    call surfs % init(dict)
    call dict % kill()

    ! Build cells.
    call charToDict(dict, CELL_DEF)
    call cells % init(dict, surfs, mats)
    call dict % kill()

    ! Build universe.
    call charToDict(dict, UNI_DEF)
    call uni % init(dict, mats, fills, cells, surfs, meshes)
    call dict % kill()

    ! Set index.
    call uni % setIdx(8)

    ! Verify fill.
    @assertEqual([-3, -4, UNDEF_MAT], fills)

  end subroutine setUp

  !!
  !! Clean environment.
  !!
@After
  subroutine clean()

    call surfs % kill()
    call cells % kill()
    call mats % kill()
    call uni % kill()

  end subroutine clean

  !!
  !! Test miscellaneous functionality (of generic universe)
  !!
@Test
  subroutine test_misc()

    ! Get id
    @assertEqual(1, uni % id())

    ! Set ID
    call uni % setId(7)
    @assertEqual(7, uni % id())

  end subroutine test_misc

  !!
  !! Test entering a universe
  !!
@Test
  subroutine test_enter()
    type(coordData)                :: data
    real(defReal), dimension(3)    :: r_ref, u_ref, r, u
    real(defReal), parameter       :: TOL = 1.0E-7_defReal

    ! ** Enter into local cell 1
    data = newCoordData([ZERO, ZERO, 3.0_defReal], [ZERO, ZERO, ONE])
    call uni % enter(data)

    ! Verify location
    r_ref = [ONE, ZERO, ZERO]
    u_ref = [ONE, ZERO, ZERO]
    @assertEqual(r_ref, data % r, TOL)
    @assertEqual(u_ref, data % u, TOL)
    @assertEqual(8, data % universeIdx)
    @assertEqual(1, data % localId)
    @assertEqual(cells % getIdx(1), data % cellIdx)

    ! ** Enter into local cell 2
    data = newCoordData([2.0_defReal, ZERO, ONE], [ZERO, ONE, ZERO])
    call uni % enter(data)

    ! Verify location
    r_ref = [-ONE, ZERO, 2.0_defReal]
    u_ref = [ZERO, -ONE, ZERO]
    @assertEqual(r_ref, data % r, TOL)
    @assertEqual(u_ref, data % u, TOL)
    @assertEqual(8, data % universeIdx)
    @assertEqual(2, data % localId)
    @assertEqual(cells % getIdx(2), data % cellIdx)

    ! ** Enter into the UNDEFINED cell
    data = newCoordData([ZERO, ZERO, 6.5_defReal], [ONE, ZERO, ZERO])
    call uni % enter(data)

    ! Verify location
    r_ref = [4.5_defReal, ZERO, ZERO]
    u_ref = [ZERO, ZERO, ONE]
    @assertEqual(r_ref, data % r, TOL)
    @assertEqual(u_ref, data % u, TOL)
    @assertEqual(8, data % universeIdx)
    @assertEqual(3, data % localId)
    @assertEqual(0, data % cellIdx)

    ! Verify rotation settings in coord
    ! * Do it only once
    @assertTrue(data % isRotated)
    @assertEqual([ZERO, ZERO,  ONE], data % rotationMatrix(1, :), TOL)
    @assertEqual([ZERO, -ONE, ZERO], data % rotationMatrix(2, :), TOL)
    @assertEqual([ONE , ZERO, ZERO], data % rotationMatrix(3, :), TOL)


  end subroutine test_enter

  !!
  !! Test distance calculation
  !!
@Test
  subroutine test_distance()
    real(defReal)            :: d, ref
    type(coordData)          :: data
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! ** In local cell 1 distance to boundary.
    data = newCoordData([-ONE, ZERO, ZERO], [ONE, ZERO, ZERO], &
                        cellIdx = cells % getIdx(1), localId = 1, universeIdx = 8)
    call uni % distance(data)
    ref = 3.0_defReal
    @assertEqual(ref, data % d, TOL * ref)
    @assertEqual(surfs % getIdx(1), data % surfaceIdx)

    ! ** In local cell 2 distance to surface 2.
    data = newCoordData([7.0_defReal, ZERO, ZERO], [-ONE, ZERO, ZERO], &
                        cellIdx = cells % getIdx(2), localId = 2, universeIdx = 8)
    call uni % distance(data)
    ref = 2.0_defReal
    @assertEqual(ref, data % d, TOL * ref)
    @assertEqual(surfs % getIdx(2), data % surfaceIdx)

    ! ** In local cell 2 distance to infinity.
    ! surfIdx must be set to 0.
    data = newCoordData([7.0_defReal, ZERO, ZERO], [ONE, ZERO, ZERO], &
                        cellIdx = cells % getIdx(2), localId = 2, universeIdx = 8)
    call uni % distance(data)
    @assertEqual(INF, data % d)
    @assertEqual(0, data % surfaceIdx)

  end subroutine test_distance

  !!
  !! Test cell-to cell crossing
  !!
@Test
  subroutine test_cross()
    type(coordData)   :: data

    ! Cross from cell 1 to cell 2.
    data = newCoordData([ZERO, 2.0_defReal, ZERO], [ZERO, ONE, ZERO], cells % getIdx(1), 1, surfs % getIdx(1), 8)
    call uni % cross(data)
    @assertEqual(2, data % localId)
    @assertEqual(cells % getIdx(2), data % cellIdx)

  end subroutine test_cross

end module cellUniverse_test