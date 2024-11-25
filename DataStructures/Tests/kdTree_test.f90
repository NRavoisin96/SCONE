module kdTree_test
  use numPrecision
  use kdTree_class, only : kdTree
  use funit
  
  implicit none
  
  ! Variable.
  type(kdTree) :: tree
contains
  !!
  !! Setup environment.
  !!
@Before
  subroutine setUp()
    real(defReal), dimension(18, 3) :: allCoordinates
    ! Populate grid of 3-D coordinates.
    allCoordinates(1, 1) = -1.0_defReal
    allCoordinates(1, 2) = -1.0_defReal
    allCoordinates(1, 3) = -1.0_defReal
    allCoordinates(2, 1) = 0.0_defReal
    allCoordinates(2, 2) = -1.0_defReal
    allCoordinates(2, 3) = -1.0_defReal
    allCoordinates(3, 1) = -1.0_defReal
    allCoordinates(3, 2) = 0.0_defReal
    allCoordinates(3, 3) = -1.0_defReal
    allCoordinates(4, 1) = 0.0_defReal
    allCoordinates(4, 2) = 0.0_defReal
    allCoordinates(4, 3) = -1.0_defReal
    allCoordinates(5, 1) = -1.0_defReal
    allCoordinates(5, 2) = -1.0_defReal
    allCoordinates(5, 3) = 1.0_defReal
    allCoordinates(6, 1) = 0.0_defReal
    allCoordinates(6, 2) = -1.0_defReal
    allCoordinates(6, 3) = 1.0_defReal
    allCoordinates(7, 1) = -1.0_defReal
    allCoordinates(7, 2) = 0.0_defReal
    allCoordinates(7, 3) = 1.0_defReal
    allCoordinates(8, 1) = 0.0_defReal
    allCoordinates(8, 2) = 0.0_defReal
    allCoordinates(8, 3) = 1.0_defReal
    allCoordinates(9, 1) = -1.0_defReal
    allCoordinates(9, 2) = 1.0_defReal
    allCoordinates(9, 3) = -1.0_defReal
    allCoordinates(10, 1) = 0.0_defReal
    allCoordinates(10, 2) = 1.0_defReal
    allCoordinates(10, 3) = -1.0_defReal
    allCoordinates(11, 1) = -1.0_defReal
    allCoordinates(11, 2) = 1.0_defReal
    allCoordinates(11, 3) = 1.0_defReal
    allCoordinates(12, 1) = 0.0_defReal
    allCoordinates(12, 2) = 1.0_defReal
    allCoordinates(12, 3) = 1.0_defReal
    allCoordinates(13, 1) = 1.0_defReal
    allCoordinates(13, 2) = 0.0_defReal
    allCoordinates(13, 3) = -1.0_defReal
    allCoordinates(14, 1) = 1.0_defReal
    allCoordinates(14, 2) = 1.0_defReal
    allCoordinates(14, 3) = -1.0_defReal
    allCoordinates(15, 1) = 1.0_defReal
    allCoordinates(15, 2) = 0.0_defReal
    allCoordinates(15, 3) = 1.0_defReal
    allCoordinates(16, 1) = 1.0_defReal
    allCoordinates(16, 2) = 1.0_defReal
    allCoordinates(16, 3) = 1.0_defReal
    allCoordinates(17, 1) = 1.0_defReal
    allCoordinates(17, 2) = -1.0_defReal
    allCoordinates(17, 3) = -1.0_defReal
    allCoordinates(18, 1) = 1.0_defReal
    allCoordinates(18, 2) = -1.0_defReal
    allCoordinates(18, 3) = 1.0_defReal
    ! Build tree.
    call tree % init(allCoordinates, transposeData = .true.)
  end subroutine setUp
  !!
  !! Clean environment.
  !!
@After
  subroutine cleanUp()
    call tree % kill()
  end subroutine cleanUp
  
  !!
  !! Test nearest neighbour searches.
  !!
@Test
  subroutine nearestNeighbour_test()
    real(defReal), dimension(3) :: r_test
    integer(shortInt)           :: idx
    ! Create test coordinates.
    r_test = [0.0_defReal, -1.1_defReal, 1.2_defReal]
    ! Find nearest vertex in the tree.
    idx = tree % findNearestVertex(r_test)
    @assertEqual(6, idx)
    ! Create test coordinates.
    r_test = [1.05_defReal, -1.3_defReal, 1.1_defReal]
    ! Find nearest vertex in the tree.
    idx = tree % findNearestVertex(r_test)
    @assertEqual(18, idx)
  end subroutine nearestNeighbour_test
end module kdTree_test
