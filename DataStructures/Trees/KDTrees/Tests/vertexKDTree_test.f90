module vertexKDTree_test
  use numPrecision
  use vertexKDTree_class, only : vertexKDTree
  use funit
  
  implicit none
  
  ! Variable.
  type(vertexKDTree) :: tree
contains
  !!
  !! Setup environment.
  !!
@Before
  subroutine setUp()
    real(defReal), dimension(3, 18) :: allCoordinates
    ! Populate grid of 3-D coordinates.
    allCoordinates(1, 1) = -1.0_defReal
    allCoordinates(2, 1) = -1.0_defReal
    allCoordinates(3, 1) = -1.0_defReal
    allCoordinates(1, 2) = 0.0_defReal
    allCoordinates(2, 2) = -1.0_defReal
    allCoordinates(3, 2) = -1.0_defReal
    allCoordinates(1, 3) = -1.0_defReal
    allCoordinates(2, 3) = 0.0_defReal
    allCoordinates(3, 3) = -1.0_defReal
    allCoordinates(1, 4) = 0.0_defReal
    allCoordinates(2, 4) = 0.0_defReal
    allCoordinates(3, 4) = -1.0_defReal
    allCoordinates(1, 5) = -1.0_defReal
    allCoordinates(2, 5) = -1.0_defReal
    allCoordinates(3, 5) = 1.0_defReal
    allCoordinates(1, 6) = 0.0_defReal
    allCoordinates(2, 6) = -1.0_defReal
    allCoordinates(3, 6) = 1.0_defReal
    allCoordinates(1, 7) = -1.0_defReal
    allCoordinates(2, 7) = 0.0_defReal
    allCoordinates(3, 7) = 1.0_defReal
    allCoordinates(1, 8) = 0.0_defReal
    allCoordinates(2, 8) = 0.0_defReal
    allCoordinates(3, 8) = 1.0_defReal
    allCoordinates(1, 9) = -1.0_defReal
    allCoordinates(2, 9) = 1.0_defReal
    allCoordinates(3, 9) = -1.0_defReal
    allCoordinates(1, 10) = 0.0_defReal
    allCoordinates(2, 10) = 1.0_defReal
    allCoordinates(3, 10) = -1.0_defReal
    allCoordinates(1, 11) = -1.0_defReal
    allCoordinates(2, 11) = 1.0_defReal
    allCoordinates(3, 11) = 1.0_defReal
    allCoordinates(1, 12) = 0.0_defReal
    allCoordinates(2, 12) = 1.0_defReal
    allCoordinates(3, 12) = 1.0_defReal
    allCoordinates(1, 13) = 1.0_defReal
    allCoordinates(2, 13) = 0.0_defReal
    allCoordinates(3, 13) = -1.0_defReal
    allCoordinates(1, 14) = 1.0_defReal
    allCoordinates(2, 14) = 1.0_defReal
    allCoordinates(3, 14) = -1.0_defReal
    allCoordinates(1, 15) = 1.0_defReal
    allCoordinates(2, 15) = 0.0_defReal
    allCoordinates(3, 15) = 1.0_defReal
    allCoordinates(1, 16) = 1.0_defReal
    allCoordinates(2, 16) = 1.0_defReal
    allCoordinates(3, 16) = 1.0_defReal
    allCoordinates(1, 17) = 1.0_defReal
    allCoordinates(2, 17) = -1.0_defReal
    allCoordinates(3, 17) = -1.0_defReal
    allCoordinates(1, 18) = 1.0_defReal
    allCoordinates(2, 18) = -1.0_defReal
    allCoordinates(3, 18) = 1.0_defReal
    ! Build tree.
    call tree % init(allCoordinates)
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
    
    ! Point is very close to vertex 6 [0, -1, 1]
    r_test = [0.0_defReal, -1.1_defReal, 1.2_defReal]
    idx = tree % findNearestVertex(r_test)
    @assertEqual(6, idx)
    
    ! Point is very close to vertex 18 [1, -1, 1]
    r_test = [1.05_defReal, -1.3_defReal, 1.1_defReal]
    idx = tree % findNearestVertex(r_test)
    @assertEqual(18, idx)

    ! Point is very close to vertex 3 [-1, 0, -1]
    r_test = [-1.1_defReal, 0.1_defReal, -0.9_defReal]
    idx = tree % findNearestVertex(r_test)
    @assertEqual(3, idx)

    ! Point is near the center, close to vertex 8 [0, 0, 1]
    r_test = [0.1_defReal, 0.1_defReal, 0.9_defReal]
    idx = tree % findNearestVertex(r_test)
    @assertEqual(8, idx)

  end subroutine nearestNeighbour_test

end module vertexKDTree_test