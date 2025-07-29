module rootUniverse_test

  use cellShelf_class,    only : cellShelf
  use charMap_class,      only : charMap
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use funit
  use meshShelf_class,    only : meshShelf
  use numPrecision
  use publicObjects,      only : coordData, newCoordData
  use rootUniverse_class, only : rootUniverse
  use surfaceShelf_class, only : surfaceShelf
  use universalVariables, only : ONE, OUTSIDE_MAT, ZERO

  implicit none

  ! Parameters
  character(*), parameter :: SURF_DEF = &
  "surf1 { id 4; type sphere; origin (0.0 5.0 0.0); radius 0.5;}&
  &surf2 { id 1; type sphere; origin (0.0 0.0 0.0); radius 2;}"

  character(*), parameter :: UNI_DEF = &
  "id 1; type rootUniverse; border 1; fill u<17>;"

  ! Variables
  type(surfaceShelf) :: surfs
  type(cellShelf)    :: cells
  type(meshShelf)    :: meshes
  type(charMap)      :: mats
  type(rootUniverse) :: uni

contains

  !!
  !! Setup environment
  !!
@Before
  subroutine setUp()
    integer(shortInt), dimension(:), allocatable :: fills
    type(dictionary) :: dict

    ! Build surfaces and MATS

    call charToDict(dict, SURF_DEF)
    call surfs % init(dict)
    call dict % kill()

    ! Build universe
    call charToDict(dict, UNI_DEF)
    call uni % init(dict, mats, fills, cells, surfs, meshes)
    call dict % kill()

    ! Set index
    call uni % setIdx(8)

    ! Verify fill
    @assertEqual([-17, OUTSIDE_MAT], fills)

  end subroutine setUp

  !!
  !! Clean environment
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

    ! Test boundary surface
    @assertEqual(surfs % getIdx(1), uni % border())

  end subroutine test_misc

  !!
  !! Test entering a universe
  !!
@Test
  subroutine test_enter()
    type(coordData)             :: data
    real(defReal), dimension(3) :: r, r_ref, u, u_ref
    real(defReal), parameter    :: TOL = 1.0E-7_defReal

    ! Enter inside
    r = [ONE, -ONE, ONE]
    u = [ONE, ZERO, ZERO]
    data = newCoordData(r, u)
    call uni % enter(data)

    r_ref = r
    u_ref = u
    @assertEqual(r_ref, data % r, TOL)
    @assertEqual(u_ref, data % u, TOL)
    @assertEqual(8, data % universeIdx)
    @assertEqual(0, data % cellIdx)
    @assertEqual(1, data % localId)

    ! Enter outside
    r = [2.0_defReal, -2.0_defReal, ONE]
    u = [ONE, ZERO, ZERO]
    data = newCoordData(r, u)
    call uni % enter(data)

    r_ref = r
    u_ref = u
    @assertEqual(r_ref, data % r, TOL)
    @assertEqual(u_ref, data % u, TOL)
    @assertEqual(8, data % universeIdx)
    @assertEqual(0, data % cellIdx)
    @assertEqual(2, data % localId)

  end subroutine test_enter

  !!
  !! Test distance calculation
  !!
@Test
  subroutine test_distance()
    type(coordData)          :: data
    real(defReal)            :: ref
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Distance from inside -> only relevant
    data = newCoordData([ONE, ZERO, ZERO], [ONE, ZERO, ZERO], localId = 1, universeIdx = 8)
    call uni % distance(data)
    ref = ONE
    @assertEqual(ref, data % d, ref * TOL)
    @assertEqual(surfs % getIdx(1), data % surfaceIdx)

  end subroutine test_distance

  !!
  !! Test cell-to cell crossing
  !!
@Test
  subroutine test_cross()
    type(coordData) :: data

    ! Cross into outside
    data = newCoordData([2.0_defReal, ZERO, ZERO], [ONE, ZERO, ZERO], localId = 1, surfaceIdx = surfs % getIdx(1), universeIdx = 8)
    call uni % cross(data)
    @assertEqual(2, data % localId)

  end subroutine test_cross

end module rootUniverse_test