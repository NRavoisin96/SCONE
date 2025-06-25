module polygon_class
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use numPrecision,
  use edgeShelf_class,              only : edgeShelf
  use face_inter,                   only : face, faceBox, kill_super => kill
  use genericProcedures,            only : append, areEqual, computeTriangleArea, computeTriangleCentre, &
                                           computeTriangleNormal, fatalError, findCommon, numToChar
  use triangle_class,               only : triangle
  use universalVariables,           only : INF, HALF, THIRD, SURF_TOL, ZERO
  use vertexShelf_class,            only : vertexShelf
  
  implicit none
  private
  
  !!
  !! Face of an OpenFOAM mesh. Consists of a list of vertices indices making the face up and 
  !! face-to-element connectivity information. Also contains a list of triangles into which the face
  !! is decomposed.
  !!
  type, public, extends(face) :: polygon
    private
  contains
    procedure                 :: computeComponents
    procedure                 :: createTriangle
    procedure                 :: kill
  end type polygon

contains
  
  !! Subroutine 'computeAreaAndNormal'
  !!
  !! Basic description:
  !!   Computes the area and the normal vector of the face.
  !!
  !! Detailed description:
  !!   First estimates the centroid of the face by performing a simple arithmetic average of the 
  !!   vertices' coordinates. Using this centroid, the face is then decomposed into a number of 
  !!   triangles equal to the number of edges in the face where each triangle shares the face 
  !!   centroid as a common point. The normal vector of each triangle is then computed and the 
  !!   normal vector for the face is then obtained by summing the normal vectors for each triangle.
  !!   The overall face area is obtained by normalising the overall normal vector. The face centroid 
  !!   is then recomputed by performing an area-weighted average of the triangles' centres.
  !!
  !! Arguments:
  !!   vertices -> A vertexShelf.
  !!
  !! Errors:
  !!   fatalError if area is negative.
  !!   fatalError if area is infinite.
  !!
  pure subroutine computeComponents(self, vertexIdxs, vertices, centroid, normal, area)
    class(polygon), intent(inout)               :: self
    integer(shortInt), dimension(:), intent(in) :: vertexIdxs
    type(vertexShelf), intent(in)               :: vertices
    real(defReal), dimension(3), intent(out)    :: centroid, normal
    real(defReal), intent(out)                  :: area
    integer(shortInt)                           :: i, nVertices
    real(defReal)                               :: norm, sumAreas
    real(defReal), dimension(3)                 :: C, sumNormals, sumAreasCentroid
    real(defReal), dimension(3, 3)              :: array
    
    ! Retrieve the number of vertices in the face and compute its centroid by performing 
    ! an arithmetic average of the vertices' coordinates.
    nVertices = size(vertexIdxs)
    
    ! Initialise C = ZERO and loop over all vertices.
    C = ZERO
    do i = 1, nVertices
      C = C + vertices % getVertexCoordinates(vertexIdxs(i))

    end do
    
    ! Normalise C, initialise sumAreas = ZERO and sumAreasCentroid = ZERO and loop over all vertices.
    array(3, :) = C / nVertices
    sumNormals = ZERO
    sumAreas = ZERO
    sumAreasCentroid = ZERO
    do i = 1, nVertices
      ! Set the vectors pointing to the remaining two vertices in the triangle and compute the triangle's
      ! centre and normal vector.
      array(1, :) = vertices % getVertexCoordinates(vertexIdxs(i))
      array(2, :) = vertices % getVertexCoordinates(vertexIdxs(mod(i, nVertices) + 1))
      
      ! Retrieve current triangle's normal and update sumNormals.
      normal = computeTriangleNormal(array)
      sumNormals = sumNormals + normal
      
      ! Compute Euclidian norm and update sumAreas and sumAreasCentroid.
      norm = norm2(normal)
      sumAreas = sumAreas + norm
      sumAreasCentroid = sumAreasCentroid + norm * sum(array, 1)

    end do

    ! Compute area.
    area = HALF * sumAreas
    
    ! Compute the centroid and normalised normal vector then set face area, centroid and normal vector.
    centroid = THIRD * sumAreasCentroid / sumAreas
    normal = sumNormals / norm2(sumNormals)

  end subroutine computeComponents

  !!
  !!
  !!
  pure subroutine createTriangle(self, lastNewFaceIdx, edgeIdxs, newVertices, newTriangle, vertexIdxs, boundingBox)
    class(polygon), intent(in)                     :: self
    integer(shortInt), intent(in)                  :: lastNewFaceIdx
    integer(shortInt), dimension(3), intent(in)    :: edgeIdxs
    type(vertexShelf), intent(in)                  :: newVertices
    type(faceBox), intent(inout)                   :: newTriangle
    integer(shortInt), dimension(3), intent(inout) :: vertexIdxs
    type(axisAlignedBoundingBox), intent(in)       :: boundingBox

    ! Build the new triangle.
    allocate(triangle :: newTriangle % item)
    call newTriangle % item % build(lastNewFaceIdx, self % getIdx(), self % getIsBoundary(), &
                                    vertexIdxs, newVertices, 'Triangle', boundingBox, edgeIdxs = edgeIdxs)

  end subroutine createTriangle
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(polygon), intent(inout) :: self
    
    ! Face.
    call kill_super(self)

  end subroutine kill
  
end module polygon_class