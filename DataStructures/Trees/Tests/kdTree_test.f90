module kdTree_test
  
  use funit
  use kdTree_class,                 only : kdTree
  use kdTreeNode_class,             only : buildKDTreeNodePayload
  use numPrecision
  use topologicalObjectShelf_class, only : topologicalObjectShelf
  use vertex_class,                 only : buildVertexPayload
  
  implicit none
  
  ! Variable.
  type(kdTree) :: tree
contains
  !!
  !! Setup environment.
  !!
@Before
  subroutine setUp()
    type(buildVertexPayload), dimension(18) :: vertexPayloads
    type(topologicalObjectShelf), target    :: vertices
    type(buildKDTreeNodePayload)            :: payload
    integer(shortInt)                       :: i
    
    ! Populate vertex payloads.
    do i = 1, size(vertexPayloads)
      vertexPayloads(i) % idx = i

    end do
    vertexPayloads(1) % coordinates = -ONE
    vertexPayloads(2) % coordinates = [ZERO, -ONE, -ONE]
    vertexPayloads(3) % coordinates = [-ONE, ZERO, -ONE]
    vertexPayloads(4) % coordinates = [ZERO, ZERO, -ONE]
    vertexPayloads(5) % coordinates = [-ONE, -ONE, ONE]
    vertexPayloads(6) % coordinates = [ZERO, -ONE, ONE]
    vertexPayloads(7) % coordinates = [-ONE, ZERO, ONE]
    vertexPayloads(8) % coordinates = [ZERO, ZERO, ONE]
    vertexPayloads(9) % coordinates = [-ONE, ONE, -ONE]
    vertexPayloads(10) % coordinates = [ZERO, ONE, -ONE]
    vertexPayloads(11) % coordinates = [-ONE, ONE, ONE]
    vertexPayloads(12) % coordinates = [ZERO, ONE, ONE]
    vertexPayloads(13) % coordinates = [ONE, ZERO, -ONE]
    vertexPayloads(14) % coordinates = [ONE, ONE, -ONE]
    vertexPayloads(15) % coordinates = [ONE, ZERO, ONE]
    vertexPayloads(16) % coordinates = ONE
    vertexPayloads(17) % coordinates = [ONE, -ONE, -ONE]
    vertexPayloads(18) % coordinates = [ONE, -ONE, ONE]
    
    ! Build tree.
    call vertices % init(vertexPayloads)
    payload % shelf => vertices
    payload % lowerBound = 1
    payload % upperBound = 18
    payload % idxs = [(i, i = 1, 18)]
    payload % bucketSize = 4
    call tree % init(payload)

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
    idx = tree % findNearestObject(r_test)
    @assertEqual(6, idx)
    
    ! Point is very close to vertex 18 [1, -1, 1]
    r_test = [1.05_defReal, -1.3_defReal, 1.1_defReal]
    idx = tree % findNearestObject(r_test)
    @assertEqual(18, idx)

    ! Point is very close to vertex 3 [-1, 0, -1]
    r_test = [-1.1_defReal, 0.1_defReal, -0.9_defReal]
    idx = tree % findNearestObject(r_test)
    @assertEqual(3, idx)

    ! Point is near the center, close to vertex 8 [0, 0, 1]
    r_test = [0.1_defReal, 0.1_defReal, 0.9_defReal]
    idx = tree % findNearestObject(r_test)
    @assertEqual(8, idx)

  end subroutine nearestNeighbour_test

end module kdTree_test