module cylindricalMap_test
  
  use cylindricalMap_class,       only : cylindricalMap
  use dictionary_class,           only : dictionary
  use funit
  use numPrecision
  use outputFile_class,           only : outputFile
  use transportObjectState_class, only : transportObjectState

  implicit none


@testCase
  type, extends(TestCase) :: test_cylindricalMap
    private
    type(cylindricalMap) :: map_radial, map_unstruct, map_3d
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_cylindricalMap

contains
@Before
  !!
  !! Sets up test_cylindricalMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_cylindricalMap), intent(inout) :: this
    type(dictionary)                          :: tempDict

    ! Build radial map with different orientation & minimum radius
    call tempDict % init(5)
    call tempDict % store('orientation', 'x')
    call tempDict % store('rGrid', 'equivolume')
    call tempDict % store('Rmin', 2.0_defReal)
    call tempDict % store('Rmax', 10.0_defReal)
    call tempDict % store('rN', 5)

    call this % map_radial % init(tempDict)
    call tempDict % kill()

    ! Build map with different origin & unstruct bins
    call tempDict % init(3)
    call tempDict % store('origin', [ONE, ONE, TWO])
    call tempDict % store('rGrid', 'unstruct')
    call tempDict % store('bins', [1.5_defReal, 2.3_defReal, 3.8_defReal, 8.0_defReal])

    call this % map_unstruct % init(tempDict)
    call tempDict % kill()

    ! Build map with radial, axial and azimuthal bins
    call tempDict % init(8)
    call tempDict % store('rGrid', 'lin')
    call tempDict % store('Rmax', 20.0_defReal)
    call tempDict % store('rN', 10)
    call tempDict % store('axGrid', 'lin')
    call tempDict % store('axN', 4)
    call tempDict % store('axMin', 2.0_defReal)
    call tempDict % store('axMax', 8.0_defReal)
    call tempDict % store('azimuthalN', 4)

    call this % map_3d % init(tempDict)
    call tempDict % kill()

  end subroutine setUp

@After
  !!
  !! Kills test_cylindricalMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_cylindricalMap), intent(inout) :: this

    call this % map_radial % kill()
    call this % map_unstruct % kill()
    call this % map_3d % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !!
  !! Test radial map with different orientation & minimum radius
  !!
@Test
  subroutine testRadial(this)
    class(test_cylindricalMap), intent(inout)  :: this
    integer(shortInt)                          :: i
    integer(shortInt), dimension(4)            :: idxs
    real(defReal), dimension(3, 4)             :: positions
    type(transportObjectState), dimension(4)   :: states
    integer(shortInt), dimension(4), parameter :: RES_IDXS = [0, 2, 4, 5]
    real(defReal), dimension(4), parameter     :: phi = [1.4_defReal, 3.0_defReal, 0.5_defReal, PI / TWO], &
                                                  r = [0.4_defReal, 5.38_defReal, 8.9_defReal, 9.1_defReal], &
                                                  z = [1.0_defReal, 39.8_defReal, 0.05_defReal, -12.2_defReal]

    ! Initialise states
    positions = ZERO
    do i = 1, 4
      positions(:, i) = [z(i), r(i) * cos(phi(i)), r(i) * sin(phi(i))]
      call states(i) % setGlobalPosition(positions(:, i))
      idxs(i) = this % map_radial % map(states(i))

    end do
    @assertEqual(RES_IDXS, idxs)

  end subroutine testRadial

  !!
  !! Test map with different origin & unstruct bins
  !!
@Test
  subroutine testUnstruct(this)
    class(test_cylindricalMap), intent(inout)  :: this
    integer(shortInt)                          :: i
    integer(shortInt), dimension(4)            :: idxs
    real(defReal), dimension(3, 4)             :: positions
    type(transportObjectState), dimension(4)   :: states
    integer(shortInt), dimension(4), parameter :: RES_IDXS = [1, 3, 0, 2]
    real(defReal), dimension(4), parameter     :: phi = [1.4_defReal, 3.0_defReal, 0.5_defReal, PI / TWO], &
                                                  r = [1.52_defReal, 5.5_defReal, 8.9_defReal, 2.88_defReal], &
                                                  z = [1.0_defReal, 39.8_defReal, 0.05_defReal, -12.2_defReal]

    ! Initialise states
    positions = ZERO
    do i = 1, 4
      positions(:, i) = [r(i) * cos(phi(i)) + ONE, r(i) * sin(phi(i)) + ONE, z(i)]
      call states(i) % setGlobalPosition(positions(:, i))
      idxs(i) = this % map_unstruct % map(states(i))

    end do
    @assertEqual(RES_IDXS, idxs)

  end subroutine testUnstruct

  !!
  !! Test map with different origin & unstruct bins
  !!
@Test
  subroutine test3d(this)
    class(test_cylindricalMap), intent(inout)  :: this
    integer(shortInt)                          :: i
    integer(shortInt), dimension(4)            :: idxs
    real(defReal), dimension(3, 4)             :: positions
    type(transportObjectState), dimension(4)   :: states
    integer(shortInt), dimension(4), parameter :: RES_IDXS = [111, 123, 0, 17]
    real(defReal), dimension(4), parameter     :: phi = [1.4_defReal, 3.0_defReal, 0.5_defReal, -PI], &
                                                  r = [1.52_defReal, 5.5_defReal, 18.9_defReal, 12.88_defReal], &
                                                  z = [6.7_defReal, 2.8_defReal, 1.1_defReal, 3.9_defReal]

    ! Initialise states
    positions = ZERO
    do i = 1, 4
      positions(:, i) = [r(i) * cos(phi(i)), r(i) * sin(phi(i)), z(i)]
      call states(i) % setGlobalPosition(positions(:, i))
      idxs(i) = this % map_3d % map(states(i))

    end do
    @assertEqual(RES_IDXS, idxs)

  end subroutine test3d

  !!
  !! Test bin number retrival
  !!
@Test
  subroutine testBinNumber(this)
    class(test_cylindricalMap), intent(inout) :: this

    @assertEqual(10, this % map_3d % bins(1), '1st dimension.')
    @assertEqual(160, this % map_3d % bins(0), 'All bins.')
    @assertEqual(4, this % map_3d % bins(3), '3rd dimension.')
    @assertEqual(4, this % map_3d % bins(2), '2nd dimension.')
    @assertEqual(1, this % map_radial % bins(2), '2nd dimension.')
    @assertEqual(0, this % map_3d % bins(4), 'Invalid dimension.')

    ! Get dimensionality
    @assertEqual(3, this % map_unstruct % dimensions())

  end subroutine testBinNumber

  !!
  !! Test correctness of print subroutine
  !! Does not check that values are correct, but that call sequence is without errors
  !!
@Test
  subroutine testPrint(this)
    class(test_cylindricalMap), intent(inout) :: this
    type(outputFile)                          :: out

    call out % init('dummyPrinter', fatalErrors = .false.)

    call this % map_radial % print(out)
    @assertTrue(out % isValid(), 'Radial map case: ')
    call out % reset()

    call this % map_unstruct % print(out)
    @assertTrue(out % isValid(), 'Unstruct map case: ')
    call out % reset()

    call this % map_3d % print(out)
    @assertTrue(out % isValid(), '3d map case: ')
    call out % reset()

  end subroutine testPrint

end module cylindricalMap_test