module polygon_class
  
  use numPrecision,
  use edgeShelf_class,     only : edgeShelf
  use face_inter,          only : face, faceBox, kill_super => kill
  use genericProcedures,   only : append, areEqual, computeTriangleArea, computeTriangleCentre, &
                                  computeTriangleNormal, fatalError, findCommon, numToChar
  use triangle_class,      only : triangle
  use universalVariables,  only : INF, HALF, THIRD, SURF_TOL, ZERO
  use vertexShelf_class,   only : vertexShelf
  
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
    procedure                 :: testForInclusion
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
    real(defReal), dimension(6), intent(in)        :: boundingBox

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
  
  !! Subroutine 'split'
  !!
  !! Basic description:
  !!   Splits the face into triangles.
  !!
  !! Detailed descrption:
  !!   Creates triangles by dividing the face from the vertex of smallest index (that is, all
  !!   triangles share this vertex). The remaining two vertices are then taken in a counter-
  !!   clockwise ordering, just as for regular faces in OpenFOAM. This ensures that the normal 
  !!   vectors of the resulting triangles are all pointing in the correct direction without the need
  !!   to check.
  !!
  !! Arguments:
  !!   triangles [inout]       -> A triangleShelf.
  !!   vertices [in]           -> A vertexShelf.
  !!   freeTriangleIdx [inout] -> Index of the first free item in triangleShelf.
  !!

  !!
  !!
  !!
  pure subroutine testForInclusion(self, vertices, intersectionCoords, diff, d, edgeIdx, vertexIdx)
    class(polygon), intent(in)                   :: self
    type(vertexShelf), intent(in)                :: vertices
    real(defReal), dimension(3), intent(in)      :: intersectionCoords, diff
    real(defReal), intent(inout)                 :: d
    integer(shortInt), intent(inout)             :: edgeIdx, vertexIdx
    integer(shortInt), dimension(:), allocatable :: vertexIdxs
    integer(shortInt)                            :: discardDimension, i, nVertices, nextIdx
    integer(shortInt), dimension(2)              :: dimensions
    real(defReal)                                :: crossProduct, sign
    real(defReal), dimension(:, :), allocatable  :: vertexCoords, projVertexCoords
    real(defReal), dimension(2)                  :: projIntersectionCoords, diffEdgeCoords, diffIntersectionCoords

    ! Compute dimension to discard and project vertices coordinates.
    vertexIdxs = self % getVertexIdxs()
    nVertices = size(vertexIdxs)
    allocate(vertexCoords(nVertices, 3))
    
    discardDimension = maxloc(abs(self % getNormal()), 1)
    dimensions = pack((/(i, i = 1, 3)/), (/(i, i = 1, 3)/) /= discardDimension)
    do i = 1, nVertices
      vertexCoords(i, :) = vertices % getVertexCoordinates(vertexIdxs(i))

    end do
    projVertexCoords = vertexCoords(:, dimensions)
    projIntersectionCoords = intersectionCoords(dimensions)

    ! Loop through all edges of the projected polygon and check that the projected intersection coordinates
    ! are on the same side of each edge. Note: this works because the polygon is convex.
    do i = 1, nVertices
      ! Pre-compute the difference in coordinates between the intersection point and the current vertex.
      diffIntersectionCoords = projIntersectionCoords - projVertexCoords(i, :)
      nextIdx = mod(i, nVertices) + 1
      
      ! Check if the intersection point is on the current or next vertex and exit if yes.
      if (i == 1 .and. areEqual(diffIntersectionCoords, ZERO)) vertexIdx = vertexIdxs(1)
      if (i < nVertices .and. areEqual(projIntersectionCoords - projVertexCoords(nextIdx, :), ZERO)) vertexIdx = vertexIdxs(nextIdx)
      if (vertexIdx > 0) exit

      ! Compute cross product.
      diffEdgeCoords = projVertexCoords(nextIdx, :) - projVertexCoords(i, :)
      crossProduct = diffIntersectionCoords(1) * diffEdgeCoords(2) - diffIntersectionCoords(2) * diffEdgeCoords(1)

      ! If crossProduct is ZERO, the point may lie on the edge.
      if (areEqual(crossProduct, ZERO)) then
        ! If point actually lies on the edge find the common edge between the two vertices and exit. Return if not.
        if (any(minval(projVertexCoords([i, nextIdx], :), dim = 1) < projIntersectionCoords .and. &
                projIntersectionCoords < maxval(projVertexCoords([i, nextIdx], :), dim = 1))) then
          edgeIdx = vertices % findCommonEdgeIdx(vertexIdxs(i), vertexIdxs(nextIdx))
          exit

        end if
        return

      end if

      ! Initialise sign.
      if (i == 1) then
        sign = crossProduct
        cycle

      end if

      ! If the cross product changes sign the intersection point is outside the polygon and we can return early.
      if (crossProduct * sign < ZERO) return

    end do

    ! If reached here, the intersection point is inside the polygon. Update d.
    d = norm2(diff)

  end subroutine testForInclusion
  
end module polygon_class