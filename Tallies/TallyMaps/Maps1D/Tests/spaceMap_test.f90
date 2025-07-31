module spaceMap_test
  
  use dictionary_class,  only : dictionary
  use funit
  use numPrecision
  use outputFile_class,  only : outputFile
  use particle_class,    only : particleState
  use spaceMap_class,    only : spaceMap
  use universalVariables

  implicit none


@testCase
  type, extends(testCase) :: test_spaceMap
    private
    type(spaceMap), dimension(3) :: structuredMaps
    type(spaceMap), dimension(3) :: unstructuredMaps
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_spaceMap

contains
  !!
  !!
  !!
  subroutine setUp(this)
    class(test_spaceMap), intent(inout)    :: this
    type(dictionary)                       :: tempDict
    real(defReal), dimension(*), parameter :: BIN_DIV = [-10.0_defReal, &
                                                         -8.0_defReal, &
                                                         -6.0_defReal, &
                                                         -4.0_defReal, &
                                                         -2.0_defReal, &
                                                         0.0_defReal, &
                                                         2.0_defReal, &
                                                         4.0_defReal, &
                                                         6.0_defReal, &
                                                         8.0_defReal, &
                                                         10.0_defReal]
    integer(shortInt)                      :: i

    ! Create structured grids for each axis.
    do i = 1, 3
      call tempDict % init(5)
      call tempDict % store('grid','lin')
      call tempDict % store('min', -10.0_defReal)
      call tempDict % store('max', 10.0_defReal)
      call tempDict % store('N', 20)
      select case(i)
        case(1)
          call tempDict % store('axis', 'x')

        case(2)
          call tempDict % store('axis', 'y')

        case(3)
          call tempDict % store('axis', 'z')

      end select

      call this % structuredMaps(i) % init(tempDict)
      call tempDict % kill()

    end do

    ! Create unstructured grids for each axis.
    do i = 1, 3
      call tempDict % init(3)
      call tempDict % store('grid','unstruct')
      call tempDict % store('bins', BIN_DIV)
      select case(i)
        case(1)
          call tempDict % store('axis', 'x')

        case(2)
          call tempDict % store('axis', 'y')

        case(3)
          call tempDict % store('axis', 'z')

      end select
      
      call this % unstructuredMaps(i) % init(tempDict)
      call tempDict % kill()

    end do

  end subroutine setUp

  !!
  !!
  !!
  subroutine tearDown(this)
    class(test_spaceMap), intent(inout) :: this
    integer(shortInt)                   :: i

    do i = 1, size(this % structuredMaps)
      call this % structuredMaps(i) % kill()

    end do

    do i = 1, size(this % unstructuredMaps)
      call this % unstructuredMaps(i) % kill()

    end do

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test structured grid
  !!
@Test
  subroutine testStructuredGrid(this)
    class(test_spaceMap), intent(inout)        :: this
    real(defReal), dimension(2), parameter     :: positions = [0.5_defReal, -10.1_defReal]
    integer(shortInt), dimension(2), parameter :: results = [11, 0]
    integer(shortInt)                          :: i
    integer(shortInt), dimension(2)            :: idxs
    type(particleState), dimension(2)          :: states

    do i = 1, size(this % structuredMaps)
      states % r(i) = positions
      idxs = this % structuredMaps(i) % map(states)
      @assertEqual(results, idxs)

    end do

  end subroutine testStructuredGrid

@Test
  subroutine testUnstructuredGrid(this)
    class(test_spaceMap), intent(inout)        :: this
    real(defReal), dimension(2), parameter     :: positions = [0.5_defReal, -10.1_defReal]
    integer(shortInt), dimension(2), parameter :: results = [6, 0]
    integer(shortInt)                          :: i
    integer(shortInt), dimension(2)            :: idxs
    type(particleState), dimension(2)          :: states

    do i = 1, size(this % unstructuredMaps)
      states % r(i) = positions
      idxs = this % unstructuredMaps(i) % map(states)
      @assertEqual(results, idxs)

    end do

  end subroutine testUnstructuredGrid

@Test
  subroutine testBins(this)
    class(test_spaceMap), intent(inout) :: this
    integer(shortInt)                   :: i

    ! Structured grids
    do i = 1, size(this % structuredMaps)
      @assertEqual(20, this % structuredMaps(i) % bins(1), 'Normal use')
      @assertEqual(20, this % structuredMaps(i) % bins(0), 'All bins')
      @assertEqual(0, this % structuredMaps(i) % bins(-2), 'Invalid dimension')

    end do

    ! Unstructured grids
    do i = 1, size(this % unstructuredMaps)
      @assertEqual(10, this % unstructuredMaps(i) % bins(1), 'Normal use')
      @assertEqual(10, this % unstructuredMaps(i) % bins(0), 'All bins')
      @assertEqual(0, this % unstructuredMaps(i) % bins(-2), 'Invalid dimension')

    end do

  end subroutine testBins

@Test
  subroutine testPrint(this)
    class(test_spaceMap), intent(inout) :: this
    type(outputFile)                    :: out
    integer(shortInt)                   :: i

    call out % init('dummyPrinter', fatalErrors = .false.)

    do i = 1, size(this % structuredMaps)
      call this % structuredMaps(i) % print(out)
      @assertTrue(out % isValid(), 'For map with structured grid: ')
      call out % reset()

    end do

    do i = 1, size(this % unstructuredMaps)
      call this % unstructuredMaps(i) % print(out)
      @assertTrue(out % isValid(), 'For map with unstructured grid: ')
      call out % reset()

    end do

  end subroutine testPrint

end module spaceMap_test