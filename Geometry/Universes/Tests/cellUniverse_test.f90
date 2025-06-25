module cellUniverse_test

  use numPrecision
  use genericProcedures
  use universalVariables, only : ONE, UNDEF_MAT, ZERO
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use coord_class,        only : coord
  use surfaceShelf_class, only : surfaceShelf
  use meshShelf_class,    only : meshShelf
  use cellShelf_class,    only : cellShelf
  use cellUniverse_class, only : cellUniverse
  use funit

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
    type(coord)                    :: new
    real(defReal), dimension(3)    :: r_ref, u_ref, r, u
    real(defReal), parameter       :: TOL = 1.0E-7_defReal
    real(defReal), dimension(3, 3) :: rotationMatrix

    ! ** Enter into local cell 1
    r = [ZERO, ZERO, 3.0_defReal]
    u = [ZERO, ZERO, ONE]
    call uni % enter(r, u, new)

    ! Verify location
    r_ref = [ONE, ZERO, ZERO]
    u_ref = [ONE, ZERO, ZERO]
    @assertEqual(r_ref, new % getPosition(), TOL)
    @assertEqual(u_ref, new % getDirection(), TOL)
    @assertEqual(8, new % getUniIdx())
    @assertEqual(1, new % getLocalId())
    @assertEqual(cells % getIdx(1), new % getCellIdx())

    ! ** Enter into local cell 2
    r = [2.0_defReal, ZERO, ONE]
    u = [ZERO, ONE, ZERO]
    call uni % enter(r, u, new)

    ! Verify location
    r_ref = [-ONE, ZERO, 2.0_defReal]
    u_ref = [ZERO, -ONE, ZERO]
    @assertEqual(r_ref, new % getPosition(), TOL)
    @assertEqual(u_ref, new % getDirection(), TOL)
    @assertEqual(8, new % getUniIdx())
    @assertEqual(2, new % getLocalId())
    @assertEqual(cells % getIdx(2), new % getCellIdx())

    ! ** Enter into the UNDEFINED cell
    r = [ZERO, ZERO, 6.5_defReal]
    u = [ONE, ZERO, ZERO]
    call uni % enter(r, u, new)

    ! Verify location
    r_ref = [4.5_defReal, ZERO, ZERO]
    u_ref = [ZERO, ZERO, ONE]
    @assertEqual(r_ref, new % getPosition(), TOL)
    @assertEqual(u_ref, new % getDirection(), TOL)
    @assertEqual(8, new % getUniIdx())
    @assertEqual(3, new % getLocalId())
    @assertEqual(0, new % getCellIdx())

    ! Verify rotation settings in coord
    ! * Do it only once
    @assertTrue(new % getIsRotated())
    rotationMatrix = new % getRotationMatrix()
    @assertEqual([ZERO, ZERO,  ONE], rotationMatrix(1, :), TOL)
    @assertEqual([ZERO, -ONE, ZERO], rotationMatrix(2, :), TOL)
    @assertEqual([ONE , ZERO, ZERO], rotationMatrix(3, :), TOL)


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

    ! ** In local cell 1 distance to boundary.
    call pos % setPosition([-ONE, ZERO, ZERO])
    call pos % setDirection([ONE, ZERO, ZERO])
    call pos % setUniIdx(8)
    call pos % setCellIdx(cells % getIdx(1))
    call pos % setLocalId(1)

    call uni % distance(pos, d, surfIdx)

    ref = 3.0_defReal
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(surfs % getIdx(1), surfIdx)


    ! ** In local cell 2 distance to surface 2.
    call pos % setPosition([7.0_defReal, ZERO, ZERO])
    call pos % setDirection([-ONE, ZERO, ZERO])
    call pos % setCellIdx(cells % getIdx(2))
    call pos % setLocalId(2)

    call uni % distance(pos, d, surfIdx)

    ref = 2.0_defReal
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(surfs % getIdx(2), surfIdx)

    ! ** In local cell 2 distance to infinity.
    ! surfIdx must be set to 0.
    call pos % setDirection([ONE, ZERO, ZERO])
    call uni % distance(pos, d, surfIdx)

    @assertEqual(INF, d)
    @assertEqual(0, surfIdx)

  end subroutine test_distance

  !!
  !! Test cell-to cell crossing
  !!
@Test
  subroutine test_cross()
    type(coord)       :: pos
    integer(shortInt) :: idx

    ! Cross from cell 1 to cell 2.
    call pos % setPosition([ZERO, 2.0_defReal, ZERO])
    call pos % setDirection([ZERO, ONE, ZERO])
    call pos % setUniIdx(8)
    call pos % setCellIdx(cells % getIdx(1))
    call pos % setLocalId(1)

    idx = surfs % getIdx(1)
    call uni % cross(pos, idx)

    @assertEqual(2, pos % getLocalId())
    @assertEqual(cells % getIdx(2), pos % getCellIdx())

  end subroutine test_cross

  !!
  !! Test cell offset
  !!
@Test
  subroutine test_cellOffset()
    type(coord) :: pos

    ! Cell 1.
    call pos % setPosition([ZERO, ONE, ZERO])
    call pos % setDirection([ZERO, ONE, ZERO])
    call pos % setUniIdx(8)
    call pos % setCellIdx(cells % getIdx(1))
    call pos % setLocalId(1)

    @assertEqual([ZERO, ZERO, ZERO], uni % cellOffset(pos))

    ! Cell 2.
    call pos % setPosition([-7.0_defReal, 2.0_defReal, ZERO])
    call pos % setDirection([ZERO, ONE, ZERO])
    call pos % setUniIdx(8)
    call pos % setCellIdx(cells % getIdx(2))
    call pos % setLocalId(2)

    @assertEqual([ZERO, ZERO, ZERO], uni % cellOffset(pos))

  end subroutine test_cellOffset

end module cellUniverse_test