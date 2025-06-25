module rootUniverse_test

  use numPrecision
  use universalVariables, only : ONE, OUTSIDE_MAT, ZERO
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use coord_class,        only : coord
  use surfaceShelf_class, only : surfaceShelf
  use cellShelf_class,    only : cellShelf
  use meshShelf_class,    only : meshShelf
  use rootUniverse_class, only : rootUniverse
  use funit

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
    type(coord) :: new
    real(defReal), dimension(3) :: r, u
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Enter inside
    r = [ONE, -ONE, ONE]
    u = [ONE, ZERO, ZERO]
    call uni % enter(r, u, new)

    @assertEqual(r, new % getPosition(), TOL)
    @assertEqual(u, new % getDirection(), TOL)
    @assertEqual(8, new % getUniIdx())
    @assertEqual(0, new % getCellIdx())
    @assertEqual(1, new % getLocalId())

    ! Enter outside
    r = [2.0_defReal, -2.0_defReal, ONE]
    u = [ONE, ZERO, ZERO]
    call uni % enter(r, u, new)

    @assertEqual(r, new % getPosition(), TOL)
    @assertEqual(u, new % getDirection(), TOL)
    @assertEqual(8, new % getUniIdx())
    @assertEqual(0, new % getCellIdx())
    @assertEqual(2, new % getLocalId())

  end subroutine test_enter

  !!
  !! Test distance calculation
  !!
@Test
  subroutine test_distance()
    real(defReal)            :: d, ref
    integer(shortInt)        :: surfIdx
    type(coord)              :: pos
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Distance from inside -> only relevant
    call pos % setPosition([ONE, ZERO, ZERO])
    call pos % setDirection([ONE, ZERO, ZERO])
    call pos % setUniIdx(8)
    call pos % setLocalId(1)

    call uni % distance(pos, d, surfIdx)

    ref = ONE
    @assertEqual(ref, d, ref * TOL)
    @assertEqual(surfs % getIdx(1), surfIdx)

  end subroutine test_distance

  !!
  !! Test cell-to cell crossing
  !!
@Test
  subroutine test_cross()
    type(coord)       :: pos
    integer(shortInt) :: idx

    ! Cross into outside
    call pos % setPosition([2.0_defReal, ZERO, ZERO])
    call pos % setDirection([ONE, ZERO, ZERO])
    call pos % setUniIdx(8)
    call pos % setLocalId(1)

    idx = surfs % getIdx(1)
    call uni % cross(pos, idx)

    @assertEqual(2, pos % getLocalId())

  end subroutine test_cross

  !!
  !! Test cell offset
  !!
@Test
  subroutine test_cellOffset()
    type(coord)       :: pos

    ! Inside
    call pos % setPosition([1.5_defReal, ZERO, ZERO])
    call pos % setDirection([ONE, ZERO, ZERO])
    call pos % setUniIdx(8)
    call pos % setLocalId(1)

    @assertEqual([ZERO, ZERO, ZERO], uni % cellOffset(pos))

    ! Outside
    call pos % setPosition([2.5_defReal, ZERO, ZERO])
    call pos % setDirection([ONE, ZERO, ZERO])
    call pos % setUniIdx(8)
    call pos % setLocalId(2)

    @assertEqual([ZERO, ZERO, ZERO], uni % cellOffset(pos))

  end subroutine test_cellOffset

end module rootUniverse_test