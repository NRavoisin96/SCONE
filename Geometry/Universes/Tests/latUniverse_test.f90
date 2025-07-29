module latUniverse_test

  use cellShelf_class,    only : cellShelf
  use charMap_class,      only : charMap
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use funit
  use genericProcedures
  use latUniverse_class,  only : latUniverse
  use meshShelf_class,    only : meshShelf
  use numPrecision
  use publicObjects,      only : coordData, newCoordData
  use surfaceShelf_class, only : surfaceShelf
  use universalVariables, only : HALF, ONE, UNDEF_MAT, ZERO

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
    type(coordData)             :: data
    real(defReal), dimension(3) :: r_ref, u_ref, r, u
    real(defReal), parameter    :: TOL = 1.0E-7_defReal

    ! ** 3D universe
    ! Enter inside -> Away from surface
    r = [ONE, ONE, HALF]
    u = [ZERO, ZERO, ONE]
    data = newCoordData(r, u)
    call uni1 % enter(data)
    r_ref = r
    u_ref = u
    @assertEqual(r_ref, data % r, TOL)
    @assertEqual(u_ref, data % u, TOL)
    @assertEqual(8, data % universeIdx)
    @assertEqual(12, data % localId)
    @assertEqual(0, data % cellIdx)

    ! Enter outside
    r = [1.6_defReal, HALF, HALF]
    u = [ZERO, ZERO, ONE]
    data = newCoordData(r, u)
    call uni1 % enter(data)
    r_ref = r
    u_ref = u
    @assertEqual(r_ref, data % r, TOL)
    @assertEqual(u_ref, data % u, TOL)
    @assertEqual(8, data % universeIdx)
    @assertEqual(13, data % localId)
    @assertEqual(0, data % cellIdx)

    ! Enter in a corner.
    r = [-HALF, ZERO, ZERO]
    u = [-ONE, ONE, -ONE]
    data = newCoordData(r, u)
    call uni1 % enter(data)
    r_ref = r
    u_ref = u / norm2(u)
    @assertEqual(r_ref, data % r, TOL)
    @assertEqual(u_ref, data % u, TOL)
    @assertEqual(8, data % universeIdx)
    @assertEqual(4, data % localId)
    @assertEqual(0, data % cellIdx)

    ! ** 2D Universe
    ! Enter inside -> Away from surface
    r = [HALF, HALF, 13.5_defReal]
    u = [ZERO, ZERO, ONE]
    data = newCoordData(r, u)
    call uni2 % enter(data)
    r_ref = r
    u_ref = u
    @assertEqual(r_ref, data % r, TOL)
    @assertEqual(u_ref, data % u, TOL)
    @assertEqual(3, data % universeIdx)
    @assertEqual(2, data % localId)
    @assertEqual(0, data % cellIdx)

    ! Enter outside
    r = [1.6_defReal, HALF, HALF]
    u = [ZERO, ZERO, ONE]
    data = newCoordData(r, u)
    call uni2 % enter(data)
    r_ref = r
    u_ref = u
    @assertEqual(r_ref, data % r, TOL)
    @assertEqual(u_ref, data % u, TOL)
    @assertEqual(3, data % universeIdx)
    @assertEqual(3, data % localId)
    @assertEqual(0, data % cellIdx)

    ! Enter on a face
    r = [ZERO, ZERO, ZERO]
    u = [-ONE, ONE, -ONE]
    data = newCoordData(r, u)
    call uni2 % enter(data)
    r_ref = r
    u_ref = u / norm2(u)
    @assertEqual(r_ref, data % r, TOL )
    @assertEqual(u_ref, data % u, TOL)
    @assertEqual(3, data % universeIdx)
    @assertEqual(1, data % localId)
    @assertEqual(0, data % cellIdx)


  end subroutine test_enter

  !!
  !! Test distance calculation
  !!
@Test
  subroutine test_distance()
    type(coordData)             :: data
    real(defReal)               :: ref, eps
    real(defReal), parameter    :: TOL = 1.0E-7_defReal

    ! ** 3D universe
    ! Well inside a cell
    data = newCoordData([ZERO, 0.1_defReal, HALF], [ZERO, -ZERO, ONE], &
                        cellIdx = 0, localId = 11, universeIdx = 8)
    call uni1 % distance(data)
    ref = 2.5_defReal
    @assertEqual(ref, data % d, TOL * ref)
    @assertEqual(-6, data % surfaceIdx)

    ! From outside -> miss
    data = newCoordData([-4.0_defReal, 0.1_defReal, HALF], [ONE, ONE, ZERO], &
                        cellIdx = 0, localId = 13, universeIdx = 8)
    call uni1 % distance(data)
    @assertEqual(INF, data % d)
    @assertEqual(-7, data % surfaceIdx)

    ! After a surface undershoot
    eps = HALF * SURF_TOL
    data = newCoordData([-ONE, ZERO - eps, -HALF], [ONE, ONE, ZERO], &
                        cellIdx = 0, localId = 4, universeIdx = 8)
    call uni1 % distance(data)
    ref = SQRT2 * HALF
    @assertEqual(ref, data % d, ref * TOL)
    @assertEqual(-2, data % surfaceIdx)

    ! After overshoot via a corner
    data = newCoordData([-HALF + eps, ZERO + eps, -HALF], [ONE, ONE, ZERO], &
                        cellIdx = 0, localId = 4, universeIdx = 8)
    call uni1 % distance(data)
    @assertEqual(ZERO, data % d,  TOL)
    @assertEqual(-2, data % surfaceIdx)

    !** 2D universe
    ! Well inside a cell -> Vertical
    data = newCoordData([HALF, 0.6_defReal, HALF], [ZERO, ZERO, ONE], &
                        cellIdx = 0, localId = 2, universeIdx = 3)
    call uni2 % distance(data)
    @assertEqual(INF, data % d)

    ! Well inside a cell -> Shallow hit
    data = newCoordData([HALF, 0.6_defReal, HALF], [ZERO, 0.01_defReal, ONE], &
                        cellIdx = 0, localId = 2, universeIdx = 3)
    call uni2 % distance(data)
    ref = sqrt(40.0_defReal ** 2 + 0.4_defReal ** 2)
    @assertEqual(ref, data % d, TOL * ref)
    @assertEqual(-4, data % surfaceIdx)

    ! From outside -> Hit
    data = newCoordData([-1.5_defReal, 0.6_defReal, HALF], [ONE, ZERO, ONE], &
                        cellIdx = 0, localId = 3, universeIdx = 3)
    call uni2 % distance(data)
    ref = HALF * SQRT2
    @assertEqual(ref, data % d, TOL * ref)
    @assertEqual(-7, data % surfaceIdx)

  end subroutine test_distance

  !!
  !! Test cell-to cell crossing
  !!
@Test
  subroutine test_cross()
    type(coordData) :: data

    ! *** 3D Lattice
    ! Cross inside
    data = newCoordData([-ONE, ZERO, -HALF], [-ONE, ONE, -ONE], 0, 1, -4, 8)
    call uni1 % cross(data)
    @assertEqual(4, data % localId)

    ! Cross from outside
    data = newCoordData([ONE, 2.0_defReal, -HALF], [ONE, -ONE, -ONE], 0, 13, -7, 8)
    call uni1 % cross(data)
    @assertEqual(6, data % localId)

    ! Cross to outside
    data = newCoordData([1.5_defReal, ONE, -ONE], [ONE, ZERO, ZERO], 0, 6, -2, 8)
    call uni1 % cross(data)
    @assertEqual(13, data % localId)

    ! *** 2D Lattice
    data = newCoordData([ZERO, ZERO, 16.5_defReal], [ONE, ONE, -ONE], 0, 1, -2, 3)
    call uni2 % cross(data)
    @assertEqual(2, data % localId)

    ! Cross from outside
    data = newCoordData([-ONE, -HALF, -78.5_defReal], [ONE, ONE, ZERO], 0, 3, -7, 3)
    call uni2 % cross(data)
    @assertEqual(1, data % localId)

  end subroutine test_cross

  !!
  !! Test cell offset
  !!
@Test
  subroutine test_cellOffset()
    type(coordData)             :: data
    real(defReal), dimension(3) :: ref
    real(defReal), parameter    :: TOL = 1.0E-7_defReal

    ! ** 3D lattice
    ! Inside
    data = newCoordData([ZERO, ZERO, HALF], [-ONE, ONE, -ONE], 0, 11, 8)
    ref = [0.0_defReal, 1.0_defReal, 1.5_defReal]
    @assertEqual(ref, uni1 % cellOffset(data % localId), TOL)

    ! Outside
    data = newCoordData([-7.0_defReal, ZERO, HALF], [-ONE, ONE, -ONE], 0, 13, 8)
    ref = ZERO
    @assertEqual(ref, uni1 % cellOffset(data % localId), TOL)

    ! ** 2D Lattice
    data = newCoordData([HALF, ZERO, HALF], [-ONE, ONE, -ONE], 0, 2, 3)
    ref = [HALF, ZERO, ZERO]
    @assertEqual(ref, uni2 % cellOffset(data % localId), TOL)

    ! Outside
    data = newCoordData([-7.0_defReal, ZERO, HALF], [-ONE, ONE, -ONE], 0, 3, 3)
    ref = ZERO
    @assertEqual(ref, uni2 % cellOffset(data % localId), TOL)

  end subroutine test_cellOffset

end module latUniverse_test