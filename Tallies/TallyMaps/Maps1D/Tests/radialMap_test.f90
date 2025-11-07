module radialMap_test

  use dictionary_class,           only : dictionary
  use funit
  use numPrecision
  use outputFile_class,           only : outputFile
  use radialMap_class,            only : radialMap
  use transportObjectState_class, only : transportObjectState

  implicit none


@testCase
  type, extends(TestCase) :: test_radialMap
    private
    type(radialMap) :: map_cyl_linear, map_cyl_equivol, map_cyl_unstruct, &
                       map_sph_from_zero, map_sph_from_min, map_sph_equivol
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_radialMap

contains
@Before
  !!
  !! Sets up test_radialMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_radialMap), intent(inout) :: this
    type(dictionary)                     :: tempDict

    ! Build cylindrical map with linear bins
    call tempDict % init(4)
    call tempDict % store('axis','x')
    call tempDict % store('grid','lin')
    call tempDict % store('max', 8.0_defReal)
    call tempDict % store('N', 4)

    call this % map_cyl_linear % init(tempDict)
    call tempDict % kill()

    ! Build cylindrical map with different orientation & minimum radius
    call tempDict % init(5)
    call tempDict % store('axis','z')
    call tempDict % store('grid','equivolume')
    call tempDict % store('min', 2.0_defReal)
    call tempDict % store('max', 10.0_defReal)
    call tempDict % store('N', 5)

    call this % map_cyl_equivol % init(tempDict)
    call tempDict % kill()

    ! Build cylindrical map with different origin & unstruct bins
    call tempDict % init(4)
    call tempDict % store('axis','z')
    call tempDict % store('origin',[ONE, ONE, ZERO])
    call tempDict % store('grid','unstruct')
    call tempDict % store('bins', [1.5_defReal, 2.3_defReal, 3.8_defReal, 8.0_defReal])

    call this % map_cyl_unstruct % init(tempDict)
    call tempDict % kill()

    ! Build spherical map with default origin & minimum radius
    call tempDict % init(3)
    call tempDict % store('grid','lin')
    call tempDict % store('max', 10.0_defReal)
    call tempDict % store('N', 20)

    call this % map_sph_from_zero % init(tempDict)
    call tempDict % kill()

    ! Build spherical map with diffrent origin & minimum radius
    call tempDict % init(5)
    call tempDict % store('origin', [ONE, ONE, ONE])
    call tempDict % store('grid', 'lin')
    call tempDict % store('min', 5.0_defReal)
    call tempDict % store('max', 10.0_defReal)
    call tempDict % store('N', 5)

    call this % map_sph_from_min % init(tempDict)
    call tempDict % kill()

    ! Build spherical map with equivolume bins
    call tempDict % init(4)
    call tempDict % store('grid', 'equivolume')
    call tempDict % store('min', 2.0_defReal)
    call tempDict % store('max', 20.0_defReal)
    call tempDict % store('N', 8)

    call this % map_sph_equivol % init(tempDict)
    call tempDict % kill()

  end subroutine setUp

@After
  !!
  !! Kills test_radialMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_radialMap), intent(inout) :: this

    call this % map_cyl_linear % kill()
    call this % map_cyl_equivol % kill()
    call this % map_cyl_unstruct % kill()
    call this % map_sph_from_zero % kill()
    call this % map_sph_from_min % kill()
    call this % map_sph_equivol % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !!
  !! Test cylindrical map with different orientation & minimum radius
  !!
@Test
  subroutine testCylLinear(this)
    class(test_radialMap), intent(inout)       :: this
    integer(shortInt)                          :: i
    integer(shortInt), dimension(4)            :: idxs
    real(defReal), dimension(3, 4)             :: directions
    type(transportObjectState), dimension(4)   :: states
    integer(shortInt), dimension(4), parameter :: RES_IDXS = [1, 3, 4, 0]
    real(defReal), dimension(4), parameter     :: phi = [1.4_defReal, 3.0_defReal, 0.5_defReal, PI/2], &
                                                  r = [0.4_defReal, 5.38_defReal, 7.9_defReal, 9.1_defReal], &
                                                  z = [1.0_defReal, 39.8_defReal, 0.05_defReal, -12.2_defReal]

    ! Initialise states.
    directions = ZERO
    do i = 1, 4
      directions(:, i) = [z(i), r(i) * cos(phi(i)), r(i) * sin(phi(i))]
      call states(i) % setGlobalPosition(directions(:, i))
      idxs(i) = this % map_cyl_linear % map(states(i))

    end do
    @assertEqual(RES_IDXS, idxs)

  end subroutine testCylLinear

  !!
  !! Test cylindrical map with equivolume bins
  !!
@Test
  subroutine testCylEquivol(this)
    class(test_radialMap), intent(inout)       :: this
    integer(shortInt)                          :: i
    integer(shortInt), dimension(4)            :: idxs
    real(defReal), dimension(3, 4)             :: directions
    type(transportObjectState), dimension(4)   :: states
    integer(shortInt), dimension(4), parameter :: RES_IDXS = [0, 1, 4, 5]
    real(defReal), dimension(4), parameter     :: phi = [1.4_defReal, 3.0_defReal, 0.5_defReal, PI/2], &
                                                  r = [1.82_defReal, 4.68_defReal, 7.9_defReal, 9.01_defReal], &
                                                  z = [1.0_defReal, 39.8_defReal, 0.05_defReal, -12.2_defReal]

    ! Initialise states.
    directions = ZERO
    do i = 1, 4
      directions(:, i) = [r(i) * cos(phi(i)), r(i) * sin(phi(i)), z(i)]
      call states(i) % setGlobalPosition(directions(:, i))
      idxs(i) = this % map_cyl_equivol % map(states(i))

    end do
    @assertEqual(RES_IDXS, idxs)

  end subroutine testCylEquivol

  !!
  !! Test cylindrical map with different origin & unstruct bins
  !!
@Test
  subroutine testCylUnstruct(this)
    class(test_radialMap), intent(inout)       :: this
    integer(shortInt)                          :: i
    integer(shortInt), dimension(4)            :: idxs
    real(defReal), dimension(3, 4)             :: directions
    type(transportObjectState), dimension(4)   :: states
    integer(shortInt), dimension(4), parameter :: RES_IDXS = [1, 3, 0, 2]
    real(defReal), dimension(4), parameter     :: phi = [1.4_defReal, 3.0_defReal, 0.5_defReal, PI/2], &
                                                  r = [1.52_defReal, 5.5_defReal, 8.9_defReal, 2.88_defReal], &
                                                  z = [1.0_defReal, 39.8_defReal, 0.05_defReal, -12.2_defReal]

    ! Initialise states.
    directions = ZERO
    do i = 1, 4
      directions(:, i) = [r(i) * cos(phi(i)) + ONE, r(i) * sin(phi(i)) + ONE, z(i)]
      call states(i) % setGlobalPosition(directions(:, i))
      idxs(i) = this % map_cyl_unstruct % map(states(i))

    end do
    @assertEqual(RES_IDXS, idxs)

  end subroutine testCylUnstruct

  !!
  !! Test spherical map with default-initialised grid
  !!
@Test
  subroutine testSphFromOrigin(this)
    class(test_radialMap), intent(inout)       :: this
    integer(shortInt)                          :: i
    integer(shortInt), dimension(4)            :: idxs
    real(defReal), dimension(3, 4)             :: directions
    type(transportObjectState), dimension(4)   :: states
    integer(shortInt), dimension(4), parameter :: RES_IDXS = [1, 8, 18, 0]
    real(defReal), dimension(4), parameter     :: phi = [1.4_defReal, 3.98_defReal, HALF, PI / TWO], &
                                                  r = [0.4_defReal, 3.58_defReal, 8.9_defReal, 11.0_defReal], &
                                                  theta = [ZERO, PI / TWO, PI / 4.0_defReal, -PI / TWO]

    ! Initialise states.
    directions = ZERO
    do i = 1, 4
      directions(:, i) = [r(i) * cos(phi(i)) * sin(theta(i)), r(i) * sin(phi(i)) * sin(theta(i)), r(i) * cos(theta(i))]
      call states(i) % setGlobalPosition(directions(:, i))
      idxs(i) = this % map_sph_from_zero % map(states(i))

    end do
    @assertEqual(RES_IDXS, idxs)

  end subroutine testSphFromOrigin

  !!
  !! Test spherical map with grid with shifted origin & minimum radius
  !!
@Test
  subroutine testSphFromMin(this)
    class(test_radialMap), intent(inout)       :: this
    integer(shortInt)                          :: i
    integer(shortInt), dimension(4)            :: idxs
    real(defReal), dimension(3, 4)             :: directions
    type(transportObjectState), dimension(4)   :: states
    integer(shortInt), dimension(4), parameter :: RES_IDXS = [0, 1, 4, 0]
    real(defReal), dimension(4), parameter     :: phi = [1.4_defReal, 3.98_defReal, HALF, PI / TWO], &
                                                  r = [1.5_defReal, 5.5_defReal, 8.9_defReal, 11.0_defReal], &
                                                  theta = [ZERO, PI / TWO, PI / 4.0_defReal, -PI / TWO]

    ! Initialise states.
    directions = ZERO
    do i = 1, 4
      directions(:, i) = [r(i) * cos(phi(i)) * sin(theta(i)), r(i) * sin(phi(i)) * sin(theta(i)), r(i) * cos(theta(i))] + ONE
      call states(i) % setGlobalPosition(directions(:, i))
      idxs(i) = this % map_sph_from_min % map(states(i))

    end do
    @assertEqual(RES_IDXS, idxs)

  end subroutine testSphFromMin

  !!
  !! Test spherical map with grid with equivolume bins
  !!
@Test
  subroutine testSphEquivol(this)
    class(test_radialMap), intent(inout)       :: this
    integer(shortInt)                          :: i
    integer(shortInt), dimension(4)            :: idxs
    real(defReal), dimension(3, 4)             :: directions
    type(transportObjectState), dimension(4)   :: states
    integer(shortInt), dimension(4), parameter :: RES_IDXS = [0, 1, 7, 2]
    real(defReal), dimension(4), parameter     :: phi = [1.4_defReal, 3.98_defReal, HALF, PI / TWO], &
                                                  r = [1.5_defReal, 5.5_defReal, 18.9_defReal, 11.0_defReal], &
                                                  theta = [ZERO, PI / TWO, PI / 4.0_defReal, -PI / TWO]

    ! Initialise states.
    directions = ZERO
    do i = 1, 4
      directions(:, i) = [r(i) * cos(phi(i)) * sin(theta(i)), r(i) * sin(phi(i)) * sin(theta(i)), r(i) * cos(theta(i))]
      call states(i) % setGlobalPosition(directions(:, i))
      idxs(i) = this % map_sph_equivol % map(states(i))

    end do
    @assertEqual(RES_IDXS, idxs)

  end subroutine testSphEquivol

  !!
  !! Test bin number retrival
  !!
@Test
  subroutine testBinNumber(this)
    class(test_radialMap), intent(inout) :: this

    ! Test that map is 1D
    @assertEqual(4, this % map_cyl_linear % bins(0), 'All bins.')
    @assertEqual(4, this % map_cyl_linear % bins(1), '1st dimension.')
    @assertEqual(0, this % map_cyl_linear % bins(2), '2nd dimension.')
    @assertEqual(3, this % map_cyl_unstruct % bins(1), '1st dimension.')
    @assertEqual(5, this % map_cyl_equivol % bins(1), '1st dimension.')
    @assertEqual(20, this % map_sph_from_zero % bins(1), '1st dimension.')
    @assertEqual(20, this % map_sph_from_zero % bins(0), 'All bins.')
    @assertEqual(0,  this % map_sph_from_min % bins(2), 'Invalid dimension.')

    ! Get dimensionality
    @assertEqual(1, this % map_cyl_linear % dimensions())
    @assertEqual(1, this % map_cyl_unstruct % dimensions())
    @assertEqual(1, this % map_sph_from_min % dimensions())

  end subroutine testBinNumber

  !!
  !! Test correctness of print subroutine
  !! Does not check that values are correct, but that call sequence is without errors
  !!
@Test
  subroutine testPrint(this)
    class(test_radialMap), intent(inout) :: this
    type(outputFile)                     :: out

    call out % init('dummyPrinter', fatalErrors = .false.)

    call this % map_cyl_linear % print(out)
    @assertTrue(out % isValid(), 'Linear map case (cylindrical).')
    call out % reset()

    call this % map_cyl_equivol % print(out)
    @assertTrue(out % isValid(), 'Equivolume map case (cylindrical).')
    call out % reset()

    call this % map_cyl_unstruct % print(out)
    @assertTrue(out % isValid(), 'Unstruct map case (cylindrical).')
    call out % reset()

    call this % map_sph_from_zero % print(out)
    @assertTrue(out % isValid(), 'Linear map case from zero (spherical).')
    call out % reset()

    call this % map_sph_from_min % print(out)
    @assertTrue(out % isValid(), 'Linear map case from minimum radius (spherical).')
    call out % reset()

    call this % map_sph_equivol % print(out)
    @assertTrue(out % isValid(), 'Equivolume map case (spherical).')
    call out % reset()

  end subroutine testPrint

end module radialMap_test