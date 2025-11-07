module coord_test

  use coord_class,        only : coord
  use coordList_class,    only : coordList
  use funit
  use numPrecision
  use universalVariables, only : ONE, ZERO

  implicit none

  ! Variables
  type(coordList) :: coords

contains

  !!
  !! Set up test envoroment
  !!
  !! Coords Placed at 3 levels
  !!   Level 1 -> global
  !!   Level 2 -> Rotated universe
  !!   Level 3 -> Translated universe
  !!
@Before
  subroutine set_up()
    real(defReal), dimension(3,3) :: mat

    ! Set Nesting
    call coords % setNesting(3)
    call coords % setMaterialIdx(2)
    call coords % setUniqueId(7)

    ! Set Level 1
    call coords % setPosition([ONE, ZERO, -ONE], 1)
    call coords % setDirection([ZERO, ONE, ZERO], 1)
    call coords % setUniverseIdx(1, 1)
    call coords % setUniverseRootId(1, 1)
    call coords % setLocalId(1, 1)
    call coords % setCellIdx(1, 1)

    ! Set Level 2
    ! Rotation Y -> Z; Z -> -Y
    mat = ZERO
    mat(1,1) = ONE
    mat(3,2) = -ONE
    mat(2,3) = ONE
    call coords % setPosition([ONE, ZERO, -ONE], 2)
    call coords % setDirection([ZERO, ZERO, -ONE], 2)
    call coords % setUniverseIdx(2, 2)
    call coords % setUniverseRootId(6, 2)
    call coords % setLocalId(3, 2)
    call coords % setIsRotated(.true., 2)
    call coords % setRotationMatrix(mat, 2)
    call coords % setCellIdx(3, 2)

    ! Set Level 3
    ! Translation to origin
    call coords % setPosition([ZERO, ZERO, ZERO], 3)
    call coords % setDirection([ZERO, ZERO, -ONE], 3)
    call coords % setUniverseIdx(4, 3)
    call coords % setUniverseRootId(12, 3)
    call coords % setLocalId(2, 3)
    call coords % setCellIdx(0, 3)

  end subroutine set_up

  !!
  !! Clean test enviroment
  !!
@After
  subroutine clean_up()

    call coords % kill()

  end subroutine clean_up

  !!
  !! Test state changing procedures
  !!
@Test
  subroutine test_changing_state()

    ! Test State
    @assertTrue(coords % isPlaced())
    @assertFalse(coords % isAbove())
    @assertFalse(coords % isUninitialised())

    ! Change to above
    call coords % takeAboveGeom()
    @assertFalse(coords % isPlaced())
    @assertTrue(coords % isAbove())
    @assertFalse(coords % isUninitialised())

    ! Chenge to uninitialised
    call coords % kill()
    @assertFalse(coords % isPlaced())
    @assertFalse(coords % isAbove())
    @assertTrue(coords % isUninitialised())

  end subroutine test_changing_state

  !!
  !! Test nesting level changes & cell inquiry
  !!
@Test
  subroutine test_nesting_level()

    ! Move deeper
    call coords % addLevel()
    @assertEqual(4, coords % getNesting())
    @assertEqual(0, coords % getLowestCellIdx())

    ! Move to higher level
    call coords % decreaseLevel(2)
    @assertEqual(2, coords % getNesting())
    @assertEqual(3, coords % getLowestCellIdx())

  end subroutine test_nesting_level

  !!
  !! Test rotation
  !!
  !! Verifies only the deflection by mu !
  !!
@Test
  subroutine test_rotation()
    real(defReal)               :: mu, phi
    real(defReal), dimension(3) :: u1, u2, u3
    real(defReal), parameter    :: TOL = 1.0E-7_defReal

    mu = 0.3_defReal
    phi = 2.1_defReal

    ! Save pre-rotation direction
    u1 = coords % getDirection(1)
    u2 = coords % getDirection(2)
    u3 = coords % getDirection(3)

    call coords % rotate(mu, phi)

    ! Verify deflection
    @assertEqual(mu, dot_product(u1, coords % getDirection(1)), TOL)
    @assertEqual(mu, dot_product(u2, coords % getDirection(2)), TOL)
    @assertEqual(mu, dot_product(u3, coords % getDirection(3)), TOL)

  end subroutine test_rotation

  !!
  !! Test direction assigment
  !!
@Test
  subroutine test_direction_assigment()
    real(defReal), dimension(3) :: u1, u2, u3
    real(defReal), parameter    :: TOL = 1.0E-7_defReal

    ! Save pre-rotation direction
    u1 = coords % getDirection(1)
    u2 = coords % getDirection(2)
    u3 = coords % getDirection(3)

    ! Invert direction
    call coords % assignDirection(-u1)

    ! Verify
    @assertEqual(-u1, coords % getDirection(1))
    @assertEqual(-u2, coords % getDirection(2))
    @assertEqual(-u3, coords % getDirection(3))

  end subroutine test_direction_assigment

  !!
  !! Test Movment
  !!
@Test
  subroutine test_movement()
    real(defReal), dimension(3) :: u1, u2, u3, r1, r2, r3
    real(defReal)               :: d
    real(defReal), parameter    :: TOL = 1.0E-7_defReal

    ! Move local
    d = 0.3_defReal
    r1 = coords % getPosition(1)
    r2 = coords % getPosition(2)
    r3 = coords % getPosition(3)
    u1 = coords % getDirection(1)
    u2 = coords % getDirection(2)
    u3 = coords % getDirection(3)

    call coords % moveLocal(d, 3)

    ! Verify
    @assertEqual(r1 + d * u1, coords % getPosition(1), TOL)
    @assertEqual(r2 + d * u2, coords % getPosition(2), TOL)
    @assertEqual(r3 + d * u3, coords % getPosition(3), TOL)
    @assertTrue(coords % isPlaced())

    ! Move Global
    d = -13.0_defReal
    r1 = coords % getPosition(1)
    u1 = coords % getDirection(1)

    call coords % moveGlobal(d)

    ! Verify
    @assertEqual(r1 + d * u1, coords % getPosition(1), TOL)
    @assertTrue(coords % isAbove())

  end subroutine test_movement

  !!
  !! Test coord validation
  !!
@Test
  subroutine test_coord_valid()
    integer(shortInt) :: i, nesting

    nesting = coords % getNesting()
    do i = 1, nesting
      @assertTrue(coords % isValid(i))

    end do

    call coords % kill()

    do i = 1, nesting
      @assertFalse(coords % isValid(i))

    end do


  end subroutine test_coord_valid

end module coord_test
