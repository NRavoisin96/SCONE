module latUniverse_test

  use numPrecision
  use genericProcedures
  use universalVariables, only : HALF, ONE, UNDEF_MAT, ZERO
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use coord_class,        only : coord
  use surfaceShelf_class, only : surfaceShelf
  use cellShelf_class,    only : cellShelf
  use meshShelf_class,    only : meshShelf
  use latUniverse_class,  only : latUniverse
  use funit

  implicit none

  ! Parameters
  character(*), parameter :: UNI1_DEF = &
  "id 1; type latUniverse; origin (0.0 0.0 0.0); rotation (0.0 0.0 0.0); &
  &pitch (1.0 2.0 3.0); shape (3 2 2); padMat void; &
  &map ( 3 4 5 &
  &      7 4 8 &
  &            &
  &      1 2 3 &
  &      4 5 6); "

  character(*), parameter :: UNI2_DEF = &
  "id 2; type latUniverse; pitch (1.0 2.0 0.0); shape (2 1 0); padMat u<1>; &
  &map (1 2); "

  ! Variables
  type(surfaceShelf) :: surfs
  type(cellShelf)    :: cells
  type(meshShelf)    :: meshes
  type(charMap)      :: mats
  type(latUniverse)  :: uni1
  type(latUniverse)  :: uni2

contains

  !!
  !! Setup environment
  !!
@Before
  subroutine setUp()
    integer(shortInt), dimension(:), allocatable :: fills
    type(dictionary)   :: dict
    character(nameLen) :: name
    integer(shortInt), dimension(:), allocatable :: ref

    ! Add materials.
    name = 'void'
    call mats % add(name, 3)

    ! Build universe 1.
    call charToDict(dict, UNI1_DEF)
    call uni1 % init(dict, mats, fills, cells, surfs, meshes)
    call dict % kill()
    call uni1 % setIdx(8)

    ! Verify fill vector.
    ref = [-4, -5, -6, -1, -2, -3, -7, -4, -8, -3, -4, -5, 3]
    @assertEqual(ref, fills)

    ! Build universe 2.
    call charToDict(dict, UNI2_DEF)
    call uni2 % init(dict, mats, fills, cells, surfs, meshes)
    call dict % kill()
    call uni2 % setIdx(3)

    ref = [-1, -2, -1]
    @assertEqual(ref, fills)

  end subroutine setUp

  !!
  !! Clean environment
  !!
@After
  subroutine clean()

    call surfs % kill()
    call cells % kill()
    call mats % kill()
    call uni1 % kill()
    call uni2 % kill()

  end subroutine clean

  !!
  !! Test miscellaneous functionality (of generic universe)
  !!
@Test
  subroutine test_misc()

    ! * Single universe is fine here
    ! Get id
    @assertEqual(1, uni1 % id())

    ! Set ID
    call uni1 % setId(7)
    @assertEqual(7, uni1 % id())

  end subroutine test_misc

  !!
  !! Test entering a universe
  !!
@Test
  subroutine test_enter()
    type(coord) :: new
    real(defReal), dimension(3) :: r_ref, u_ref, r, u
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! ** 3D universe
    ! Enter inside -> Away from surface
    r = [ONE, ONE, HALF]
    u = [ZERO, ZERO, ONE]

    call uni1 % enter(r, u, new)

    r_ref = r
    u_ref = u
    @assertEqual(r_ref, new % getPosition(), TOL)
    @assertEqual(u_ref, new % getDirection(), TOL)
    @assertEqual(8, new % getUniIdx())
    @assertEqual(12, new % getLocalId())
    @assertEqual(0, new % getCellIdx())

    ! Enter outside
    r = [1.6_defReal, HALF, HALF]
    u = [ZERO, ZERO, ONE]

    call uni1 % enter(r, u, new)

    r_ref = r
    u_ref = u
    @assertEqual(r_ref, new % getPosition(), TOL)
    @assertEqual(u_ref, new % getDirection(), TOL)
    @assertEqual(8, new % getUniIdx())
    @assertEqual(13, new % getLocalId())
    @assertEqual(0, new % getCellIdx())

    ! Enter in a corner.
    r = [-HALF, ZERO, ZERO]
    u = [-ONE, ONE, -ONE]
    u = u / norm2(u)

    call uni1 % enter(r, u, new)

    r_ref = r
    u_ref = u
    @assertEqual(r_ref, new % getPosition(), TOL)
    @assertEqual(u_ref, new % getDirection(), TOL)
    @assertEqual(8, new % getUniIdx())
    @assertEqual(4, new % getLocalId())
    @assertEqual(0, new % getCellIdx())

    ! ** 2D Universe
    ! Enter inside -> Away from surface
    r = [HALF, HALF, 13.5_defReal]
    u = [ZERO, ZERO, ONE]

    call uni2 % enter(r, u, new)

    r_ref = r
    u_ref = u
    @assertEqual(r_ref, new % getPosition(), TOL)
    @assertEqual(u_ref, new % getDirection(), TOL)
    @assertEqual(3, new % getUniIdx())
    @assertEqual(2, new % getLocalId())
    @assertEqual(0, new % getCellIdx())

    ! Enter outside
    r = [1.6_defReal, HALF, HALF]
    u = [ZERO, ZERO, ONE]

    call uni2 % enter(r, u, new)

    r_ref = r
    u_ref = u
    @assertEqual(r_ref, new % getPosition(), TOL)
    @assertEqual(u_ref, new % getDirection(), TOL)
    @assertEqual(3, new % getUniIdx())
    @assertEqual(3, new % getLocalId())
    @assertEqual(0, new % getCellIdx())

    ! Enter on a face
    r = [ZERO, ZERO, ZERO]
    u = [-ONE, ONE, -ONE]
    u = u / norm2(u)

    call uni2 % enter(r, u, new)

    r_ref = r
    u_ref = u
    @assertEqual(r_ref, new % getPosition(), TOL )
    @assertEqual(u_ref, new % getDirection(), TOL)
    @assertEqual(3, new % getUniIdx())
    @assertEqual(1, new % getLocalId())
    @assertEqual(0, new % getCellIdx())


  end subroutine test_enter

  !!
  !! Test distance calculation
  !!
@Test
  subroutine test_distance()
    real(defReal)            :: d, ref, eps
    integer(shortInt)        :: surfIdx
    type(coord)              :: pos
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! ** 3D universe
    ! Well inside a cell
    call pos % setPosition([ZERO, 0.1_defReal, HALF])
    call pos % setDirection([ZERO, -ZERO, ONE])
    call pos % setUniIdx(8)
    call pos % setCellIdx(0)
    call pos % setLocalId(11)

    call uni1 % distance(pos, d, surfIdx)

    ref = 2.5_defReal
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(-6, surfIdx)

    ! From outside -> miss
    call pos % setPosition([-4.0_defReal, 0.1_defReal, HALF])
    call pos % setDirection([ONE, ONE, ZERO])
    call pos % setDirection(pos % getDirection() / norm2(pos % getDirection()))
    call pos % setUniIdx(8)
    call pos % setCellIdx(0)
    call pos % setLocalId(13)

    call uni1 % distance(pos, d, surfIdx)

    @assertEqual(INF, d)
    @assertEqual(-7, surfIdx)

    ! After a surface undershoot
    eps = HALF * SURF_TOL
    call pos % setPosition([-ONE, ZERO - eps, -HALF])
    call pos % setDirection([ONE, ONE, ZERO])
    call pos % setDirection(pos % getDirection() / norm2(pos % getDirection()))
    call pos % setUniIdx(8)
    call pos % setCellIdx(0)
    call pos % setLocalId(4)

    call uni1 % distance(pos, d, surfIdx)

    ref = SQRT2 * HALF
    @assertEqual(ref, d, ref * TOL)
    @assertEqual(-2, surfIdx)

    ! After overshoot via a corner
    call pos % setPosition([-HALF + eps, ZERO + eps, -HALF])
    call pos % setDirection([ONE, ONE, ZERO])
    call pos % setDirection(pos % getDirection() / norm2(pos % getDirection()))
    call pos % setUniIdx(8)
    call pos % setCellIdx(0)
    call pos % setLocalId(4)

    call uni1 % distance(pos, d, surfIdx)
    @assertEqual(ZERO, d,  TOL)
    @assertEqual(-2, surfIdx)

    !** 2D universe
    ! Well inside a cell -> Vertical
    call pos % setPosition([HALF, 0.6_defReal, HALF])
    call pos % setDirection([ZERO, ZERO, ONE])
    call pos % setUniIdx(3)
    call pos % setCellIdx(0)
    call pos % setLocalId(2)

    call uni2 % distance(pos, d, surfIdx)

    @assertEqual(INF, d)

    ! Well inside a cell -> Shallow hit
    call pos % setPosition([HALF, 0.6_defReal, HALF])
    call pos % setDirection([ZERO, 0.01_defReal, ONE])
    call pos % setDirection(pos % getDirection() / norm2(pos % getDirection()))

    call uni2 % distance(pos, d, surfIdx)

    ref = sqrt(40.0_defReal**2 + 0.4_defReal**2)
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(-4, surfIdx)

    ! From outside -> Hit
    call pos % setPosition([-1.5_defReal, 0.6_defReal, HALF])
    call pos % setDirection([ONE, ZERO, ONE])
    call pos % setDirection(pos % getDirection() / norm2(pos % getDirection()))
    call pos % setUniIdx(3)
    call pos % setCellIdx(0)
    call pos % setLocalId(3)

    call uni2 % distance(pos, d, surfIdx)

    ref = HALF * SQRT2
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(-7, surfIdx)

  end subroutine test_distance

  !!
  !! Test cell-to cell crossing
  !!
@Test
  subroutine test_cross()
    type(coord)       :: pos

    ! *** 3D Lattice
    ! Cross inside
    call pos % setPosition([-ONE, ZERO, -HALF])
    call pos % setDirection([-ONE, ONE, -ONE])
    call pos % setDirection(pos % getDirection() / norm2(pos % getDirection()))
    call pos % setUniIdx(8)
    call pos % setCellIdx(0)
    call pos % setLocalId(1)

    call uni1 % cross(pos, -4)

    @assertEqual(4, pos % getLocalId())

    ! Cross from outside
    call pos % setPosition([ONE, 2.0_defReal, -HALF])
    call pos % setDirection([ONE, -ONE, -ONE])
    call pos % setDirection(pos % getDirection() / norm2(pos % getDirection()))
    call pos % setLocalId(13)

    call uni1 % cross(pos, -7)

    @assertEqual(6, pos % getLocalId())

    ! Cross to outside
    call pos % setPosition([1.5_defReal, ONE, -ONE])
    call pos % setDirection([ONE, ZERO, ZERO])
    call pos % setLocalId(6)

    call uni1 % cross(pos, -2)

    @assertEqual(13, pos % getLocalId())

    ! *** 2D Lattice
    call pos % setPosition([ZERO, ZERO, 16.5_defReal])
    call pos % setDirection([ONE, ONE, -ONE])
    call pos % setDirection(pos % getDirection() / norm2(pos % getDirection()))
    call pos % setUniIdx(3)
    call pos % setCellIdx(0)
    call pos % setLocalId(1)

    call uni2 % cross(pos, -2)

    @assertEqual(2, pos % getLocalId())

    ! Cross from outside
    call pos % setPosition([-ONE, -HALF, -78.5_defReal])
    call pos % setDirection([ONE, ONE, ZERO])
    call pos % setDirection(pos % getDirection() / norm2(pos % getDirection()))
    call pos % setLocalId(3)

    call uni2 % cross(pos, -7)

    @assertEqual(1, pos % getLocalId())

  end subroutine test_cross

  !!
  !! Test cell offset
  !!
@Test
  subroutine test_cellOffset()
    type(coord)                 :: pos
    real(defReal), dimension(3) :: ref
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! ** 3D lattice
    ! Inside
    call pos % setPosition([ZERO, ZERO, HALF])
    call pos % setDirection([-ONE, ONE, -ONE])
    call pos % setDirection(pos % getDirection() / norm2(pos % getDirection()))
    call pos % setUniIdx(8)
    call pos % setCellIdx(0)
    call pos % setLocalId(11)

    ref = [0.0_defReal, 1.0_defReal, 1.5_defReal]
    @assertEqual(ref, uni1 % cellOffset(pos), TOL)

    ! Outside
    call pos % setPosition([-7.0_defReal, ZERO, HALF])
    call pos % setDirection([-ONE, ONE, -ONE])
    call pos % setDirection(pos % getDirection() / norm2(pos % getDirection()))
    call pos % setLocalId(13)

    ref = ZERO
    @assertEqual(ref, uni1 % cellOffset(pos), TOL)

    ! ** 2D Lattice
    call pos % setPosition([HALF, ZERO, HALF])
    call pos % setDirection([-ONE, ONE, -ONE])
    call pos % setDirection(pos % getDirection() / norm2(pos % getDirection()))
    call pos % setUniIdx(3)
    call pos % setCellIdx(0)
    call pos % setLocalId(2)

    ref = [HALF, ZERO, ZERO]
    @assertEqual(ref, uni2 % cellOffset(pos), TOL)

    ! Outside
    call pos % setPosition([-7.0_defReal, ZERO, HALF])
    call pos % setDirection([-ONE, ONE, -ONE])
    call pos % setDirection(pos % getDirection() / norm2(pos % getDirection()))
    call pos % setLocalId(3)

    ref = ZERO
    @assertEqual(ref, uni2 % cellOffset(pos), TOL)

  end subroutine test_cellOffset


end module latUniverse_test
