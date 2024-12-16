module face_class
  
  use numPrecision,
  use edgeShelf_class,     only : edgeShelf
  use genericProcedures,   only : append, areEqual, computeTriangleArea, computeTriangleCentre, &
                                  computeTriangleNormal, fatalError, findCommon, numToChar
  use vertex_class,        only : vertex
  use vertexShelf_class,   only : vertexShelf
  use triangleShelf_class, only : triangleShelf
  use universalVariables,  only : INF, HALF, THIRD, SURF_TOL, ZERO
  
  implicit none
  private
  
  !!
  !! Face of an OpenFOAM mesh. Consists of a list of vertices indices making the face up and 
  !! face-to-element connectivity information. Also contains a list of triangles into which the face
  !! is decomposed.
  !!
  !! Private members:
  !!   idx            -> Index of the face.
  !!   vertices       -> Array of vertices indices making the face up.
  !!   faceToElements -> Array listing the owner and neighbour elements for the face.
  !!   triangles      -> Array of triangles indices into which the face is decomposed.
  !!   boundaryFace   -> Is the face a boundary face?
  !!   area           -> Area of the face.
  !!   centroid       -> Vector pointing to the centroid of the face.
  !!   normal         -> Normal vector of the face.
  !!
  type, public                                   :: face
    private
    integer(shortInt)                            :: idx = 0
    integer(shortInt), dimension(:), allocatable :: vertices, faceToElements, triangles, edgeIdxs
    logical(defBool)                             :: boundaryFace = .false.
    real(defReal)                                :: area = ZERO
    real(defReal), dimension(3)                  :: centroid = ZERO, normal = ZERO
  contains
    procedure                                    :: addEdgeIdx
    procedure                                    :: addElementToFace
    procedure                                    :: addTriangle
    procedure                                    :: addVertexToFace
    procedure                                    :: computeAreaAndNormal
    procedure                                    :: computeIntersection
    procedure                                    :: getArea
    procedure                                    :: getCentroid
    procedure                                    :: getEdgeIdxs
    procedure                                    :: getFaceToElements
    procedure                                    :: getIdx
    procedure                                    :: getNormal
    procedure                                    :: getTriangles
    procedure                                    :: getVertices
    procedure                                    :: getIsBoundary
    procedure                                    :: kill
    procedure                                    :: setBoundaryFace
    procedure                                    :: setIdx
    procedure                                    :: split
  end type face

contains
  
  !! Subroutine 'addTriangle'
  !!
  !! Basic description:
  !!   Adds the index of a triangle in the face.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the triangle.
  !!
  elemental subroutine addTriangle(self, idx)
    class(face), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx
    
    call append(self % triangles, idx)

  end subroutine addTriangle

  !! Subroutine 'addEdgeIdx'
  !!
  !! Basic description:
  !!   Adds the index of an edge sharing the face.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the edge sharing the face.
  !!
  elemental subroutine addEdgeIdx(self, idx)
    class(face), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx

    call append(self % edgeIdxs, idx)

  end subroutine addEdgeIdx
  
  !! Subroutine 'addElementToFace'
  !!
  !! Basic description:
  !!   Adds the index of an element containing the face.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element containing the face.
  !!
  elemental subroutine addElementToFace(self, idx)
    class(face), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx
    
    call append(self % faceToElements, idx)

  end subroutine addElementToFace
  
  !! Subroutine 'addVertexToFace'
  !!
  !! Basic description:
  !!   Adds the index of a vertex in the face.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the vertex.
  !!
  elemental subroutine addVertexToFace(self, idx)
    class(face), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx
    
    call append(self % vertices, idx)

  end subroutine addVertexToFace
  
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
  !!   vertices -> An array of vertex structures making the face up.
  !!
  !! Errors:
  !!   fatalError if area is negative.
  !!   fatalError if area is infinite.
  !!
  subroutine computeAreaAndNormal(self, vertices)
    class(face), intent(inout)                                 :: self
    type(vertex), dimension(size(self % vertices)), intent(in) :: vertices
    integer(shortInt)                                          :: i, nVertices
    real(defReal)                                              :: norm, sumAreas
    real(defReal), dimension(3)                                :: C, normal, sumNormals, sumAreasCentroid
    real(defReal), dimension(3, 3)                             :: array
    character(100), parameter                                  :: Here = 'computeAreaAndNormal &
                                                                  &(face_class.f90)'
    
    ! Retrieve the number of vertices in the face and compute its centroid by performing 
    ! an arithmetic average of the vertices' coordinates. If the face is a triangle we
    ! perform a direct computation to minimise round-off errors.
    nVertices = size(vertices)
    if (nVertices == 3) then
      do i = 1, nVertices
        array(i, :) = vertices(i) % getCoordinates()

      end do
      self % centroid = computeTriangleCentre(array)
      normal = computeTriangleNormal(array)
      self % area = computeTriangleArea(normal)
      self % normal = normal / norm2(normal)

    else
      ! Initialise C = ZERO and loop over all vertices.
      C = ZERO
      do i = 1, nVertices
        C = C + vertices(i) % getCoordinates()

      end do
      
      ! Normalise C, initialise sumAreas = ZERO and sumAreasCentroid = ZERO and loop over all vertices.
      array(3, :) = C / nVertices
      sumNormals = ZERO
      sumAreas = ZERO
      sumAreasCentroid = ZERO
      do i = 1, nVertices
        ! Set the vectors pointing to the remaining two vertices in the triangle and compute the triangle's
        ! centre and normal vector.
        array(1, :) = vertices(i) % getCoordinates()
        array(2, :) = vertices(mod(i, nVertices) + 1) % getCoordinates()
        
        ! Retrieve current triangle's normal and update sumNormals.
        normal = computeTriangleNormal(array)
        sumNormals = sumNormals + normal
        
        ! Compute Euclidian norm and update sumAreas and sumAreasCentroid.
        norm = norm2(normal)
        sumAreas = sumAreas + norm
        sumAreasCentroid = sumAreasCentroid + norm * sum(array, 1)

      end do
      
      ! Normalise the centroid, compute the face area and normalise its normal vector.
      self % centroid = THIRD * sumAreasCentroid / sumAreas
      self % area = HALF * sumAreas
      self % normal = sumNormals / norm2(sumNormals)

    end if
    
    ! Check that the area is valid.
    if (self % area <= ZERO) call fatalError(Here, 'Negative area for face with index: '//numToChar(self % idx)//'.')
    if (self % area >= INF) call fatalError(Here, 'Infinite area for face with index: '//numToChar(self % idx)//'.')

  end subroutine computeAreaAndNormal

  !! Subroutine 'computeIntersection'
  !!
  !! Basic description:
  !!   Checks whether a line segment intersects the face and if so, computes the distance from
  !!   the segment's origin to the point of intersection.
  !!
  !! Detailed description:
  !!
  !! Arguments:
  !!   startPos [in]          -> 3-D coordinates of the line segment's origin.
  !!   endPos [in]            -> 3-D coordinates of the line segment's end.
  !!   firstVertexCoords [in] -> 3-D coordinates of the face's first vertex.
  !!   isIntersecting [out]   -> .true. if the line segment intersects the triangle.
  !!   d [out]                -> Distance from the line segment's origin to the point of intersection.
  !!
  pure subroutine computeIntersection(self, startPos, endPos, u, vertices, d, edgeIdx, vertexIdx)
    class(face), intent(in)                                    :: self
    real(defReal), dimension(3), intent(in)                    :: startPos, endPos, u
    type(vertex), dimension(size(self % vertices)), intent(in) :: vertices
    real(defReal), intent(out)                                 :: d
    integer(shortInt), intent(out)                             :: edgeIdx, vertexIdx
    real(defReal), dimension(3)                                :: normal, diff, intersectionCoords
    real(defReal)                                              :: denominator, s, crossProduct, sign
    integer(shortInt)                                          :: discardDimension, i, nVertices, nextIdx
    integer(shortInt), dimension(2)                            :: dimensions
    real(defReal), dimension(size(self % vertices), 3)         :: vertexCoords
    real(defReal), dimension(size(self % vertices), 2)         :: projVertexCoords
    real(defReal), dimension(2)                                :: projIntersectionCoords, diffEdgeCoords, &
                                                                  diffIntersectionCoords
    integer(shortInt), dimension(:), allocatable               :: commonEdgeIdxs

    ! Initialise d = INF, edgeIdx = 0 and vertexIdx = 0.
    d = INF
    edgeIdx = 0
    vertexIdx = 0
    
    ! Retrieve the face's normal vector and pre-compute the difference between the line segment's end
    ! and beginning positions.
    normal = self % normal
    diff = endPos - startPos
    denominator = dot_product(normal, diff)

    ! If the denominator is ZERO, return early since the line segment is parallel to the face's plane.
    ! Else, compute the fraction of the line segment required to intersect the face's plane, s.
    if (areEqual(denominator, ZERO)) return
    s = dot_product(normal, self % centroid - startPos) / denominator
    
    ! If s is ZERO, the line segment's origin is on the face. In this case return early if the segment
    ! points in the same direction as the face's normal.
    if (areEqual(s, ZERO) .and. dot_product(normal, u) >= ZERO) return
    
    ! If s < ZERO or s > ONE, return early since the intersection is outside the line segment.
    if (s < ZERO .or. s > ONE) return

    ! Compute the coordinates of the intersection point.
    diff = s * diff
    intersectionCoords = startPos + diff

    ! Compute dimension to discard and project vertices coordinates.
    discardDimension = maxloc(abs(normal), 1)
    dimensions = pack((/(i, i = 1, 3)/), (/(i, i = 1, 3)/) /= discardDimension)
    nVertices = size(self % vertices)
    do i = 1, nVertices
      vertexCoords(i, :) = vertices(i) % getCoordinates()

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
      if (i == 1 .and. areEqual(diffIntersectionCoords, ZERO)) vertexIdx = vertices(i) % getIdx()
      if (i < nVertices .and. areEqual(projIntersectionCoords - projVertexCoords(nextIdx, :), ZERO)) &
      vertexIdx = vertices(nextIdx) % getIdx()
      if (vertexIdx > 0) exit

      ! Compute cross product.
      diffEdgeCoords = projVertexCoords(nextIdx, :) - projVertexCoords(i, :)
      crossProduct = diffIntersectionCoords(1) * diffEdgeCoords(2) - diffIntersectionCoords(2) * diffEdgeCoords(1)

      ! If crossProduct is ZERO, the point may lie on the edge.
      if (areEqual(crossProduct, ZERO)) then
        ! If point actually lies on the edge find the common edge between the two vertices and exit. Return if not.
        if (any(minval(projVertexCoords([i, nextIdx], :), dim = 1) < projIntersectionCoords .and. &
                projIntersectionCoords < maxval(projVertexCoords([i, nextIdx], :), dim = 1))) then
          commonEdgeIdxs = findCommon(vertices(i) % getEdgeIdxs(), vertices(nextIdx) % getEdgeIdxs())
          edgeIdx = commonEdgeIdxs(1)
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

  end subroutine computeIntersection
  
  !! Function 'getArea'
  !!
  !! Basic description:
  !!   Returns the area of the face.
  !!
  !! Result:
  !!   area -> Area of the face.
  !!
  elemental function getArea(self) result(area)
    class(face), intent(in) :: self
    real(defReal)           :: area
    
    area = self % area

  end function getArea
  
  !! Function 'getCentroid'
  !!
  !! Basic description:
  !!   Returns the centroid vector of the face.
  !!
  !! Result:
  !!   centroid -> Vector pointing to the centroid of the face.
  !!
  pure function getCentroid(self) result(centroid)
    class(face), intent(in)     :: self
    real(defReal), dimension(3) :: centroid
    
    centroid = self % centroid

  end function getCentroid

  !! Function 'getEdgeIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the edges in the face.
  !!
  !! Result:
  !!   edgeIdxs -> Indices of the edges in the face.
  !!
  pure function getEdgeIdxs(self) result(edgeIdxs)
    class(face), intent(in)                             :: self
    integer(shortInt), dimension(size(self % edgeIdxs)) :: edgeIdxs

    edgeIdxs = self % edgeIdxs

  end function getEdgeIdxs
  
  !! Function 'getFaceToElements'
  !!
  !! Basic description:
  !!   Returns the elements containing the face.
  !!
  !! Result:
  !!   faceToElements -> Array listing the elements containing the face.
  !!
  pure function getFaceToElements(self) result(faceToElements)
    class(face), intent(in)                                   :: self
    integer(shortInt), dimension(size(self % faceToElements)) :: faceToElements
    
    faceToElements = self % faceToElements

  end function getFaceToElements
  
  !! Function 'getIdx'
  !!
  !! Basic description:
  !!   Returns the index of the face.
  !!
  !! Result:
  !!   idx -> Index of the face.
  !!
  elemental function getIdx(self) result(idx)
    class(face), intent(in) :: self
    integer(shortInt)       :: idx
    
    idx = self % idx

  end function getIdx

  !! Function 'getIsBoundary'
  !!
  !! Basic description:
  !!   Returns .true. if face is a boundary face.
  !!
  elemental function getIsBoundary(self) result(isIt)
    class(face), intent(in) :: self
    logical(defBool)        :: isIt

    isIt = self % boundaryFace

  end function getIsBoundary 
  
  !! Function 'getNormal'
  !!
  !! Basic description:
  !!   Returns the normal vector of the face.
  !!
  !! Result:
  !!   normal -> Normal vector of the face.
  !!
  pure function getNormal(self, idx) result(normal)
    class(face), intent(in)                 :: self
    real(defReal), dimension(3)             :: normal
    integer(shortInt), intent(in), optional :: idx
    
    normal = self % normal
    
    if (.not. present(idx)) return
    if (idx < 0) normal = -normal

  end function getNormal
  
  !! Function 'getTriangles'
  !!
  !! Basic description:
  !!   Returns the indices of the triangles in the face.
  !!
  !! Result:
  !!   triangleIdxs -> Array of triangle indices in the face.
  !!
  pure function getTriangles(self) result(triangleIdxs)
    class(face), intent(in)                              :: self
    integer(shortInt), dimension(size(self % triangles)) :: triangleIdxs
    
    triangleIdxs = self % triangles

  end function getTriangles
  
  !! Function 'getVertices'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in the face.
  !!
  !! Result:
  !!   vertexIdxs -> Array of vertex indices in the face.
  !!
  pure function getVertices(self) result(vertexIdxs)
    class(face), intent(in)                             :: self
    integer(shortInt), dimension(size(self % vertices)) :: vertexIdxs
    
    vertexIdxs = self % vertices

  end function getVertices
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(face), intent(inout) :: self
    
    self % idx = 0
    self % boundaryFace = .false.
    self % area = ZERO
    self % centroid = ZERO
    self % normal = ZERO
    if (allocated(self % edgeIdxs)) deallocate(self % edgeIdxs)
    if (allocated(self % vertices)) deallocate(self % vertices)
    if (allocated(self % triangles)) deallocate(self % triangles)
    if (allocated(self % faceToElements)) deallocate(self % faceToElements)

  end subroutine kill
  
  !! Subroutine 'setBoundaryFace'
  !!
  !! Basic description:
  !!   Sets the face as a boundary face.
  !!
  elemental subroutine setBoundaryFace(self)
    class(face), intent(inout) :: self
    
    self % boundaryFace = .true.

  end subroutine setBoundaryFace
  
  !! Subroutine 'setIdx'
  !!
  !! Basic description:
  !!   Sets the index of the face.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the face.
  !!
  !! Error:
  !!   fatalError if idx < 1.
  !!
  subroutine setIdx(self, idx)
    class(face), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx
    character(100), parameter     :: Here = 'setIdx (face_class.f90)'
    
    ! Catch invalid index.
    if (idx < 1) call fatalError(Here, 'Face index must be +ve. Is: '//numToChar(idx)//'.')
    self % idx = idx

  end subroutine setIdx
  
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
  !! TODO: Refactor this into helper functions.
  subroutine split(self, edges, triangles, vertices, lastEdgeIdx, freeTriangleIdx)
    class(face), intent(inout)                          :: self
    type(edgeShelf), intent(inout)                      :: edges
    type(triangleShelf), intent(inout)                  :: triangles
    type(vertexShelf), intent(inout)                    :: vertices
    integer(shortInt), intent(inout)                    :: lastEdgeIdx, freeTriangleIdx
    integer(shortInt), dimension(size(self % vertices)) :: vertexIdxs
    integer(shortInt)                                   :: nVertices, i, j, minVertexIdx, &
                                                           secondVertexIdx, thirdVertexIdx, edgeIdx
    integer(shortInt), dimension(3)                     :: triangleVertexIdxs
    real(defReal), dimension(3)                         :: minVertexCoords, secondVertexCoords, &
                                                           thirdVertexCoords
    
    ! Retrieve the number of vertices in the face.
    vertexIdxs = self % vertices                                                       
    nVertices = size(vertexIdxs)
    
    ! Pre-compute the location and coordinates of the vertex of smallest index in the face
    ! then loop through all new triangles.
    minVertexIdx = minloc(vertexIdxs, 1)
    minVertexCoords = vertices % shelf(vertexIdxs(minVertexIdx)) % getCoordinates()
    do i = 1, nVertices - 2
      ! Update freeTriangleIdx and add freeTriangleIdx to the indices of triangles in the face.
      freeTriangleIdx = freeTriangleIdx + 1
      call self % addTriangle(freeTriangleIdx)
      
      ! Create a new triangle. Set its index and parent face index.
      call triangles % shelf(freeTriangleIdx) % setIdx(freeTriangleIdx)
      call triangles % shelf(freeTriangleIdx) % setFace(self % idx)
      
      ! If the current face is a boundary face, set the current triangle as boundary triangle.
      if (self % boundaryFace) call triangles % shelf(freeTriangleIdx) % setIsBoundary()
      
      ! If the face is already a triangle there is no need to compute its vertices, centre, area and 
      ! normal vector again: simply copy them from the face and return early.
      if (nVertices == 3) then
        call triangles % shelf(freeTriangleIdx) % setVertices(vertexIdxs)
        call triangles % shelf(freeTriangleIdx) % setAB(vertices % shelf(vertexIdxs(2)) % getCoordinates() - &
                                                        vertices % shelf(vertexIdxs(1)) % getCoordinates())
        call triangles % shelf(freeTriangleIdx) % setAC(vertices % shelf(vertexIdxs(3)) % getCoordinates() - &
                                                        vertices % shelf(vertexIdxs(1)) % getCoordinates())
        call triangles % shelf(freeTriangleIdx) % setCentre(self % centroid)
        call triangles % shelf(freeTriangleIdx) % setArea(self % area)
        call triangles % shelf(freeTriangleIdx) % setNormal(self % normal)
        do j = 1, 3
          edgeIdx = self % edgeIdxs(j)
          call edges % shelf(edgeIdx) % addTriangleIdx(freeTriangleIdx)
          call triangles % shelf(freeTriangleIdx) % addEdgeIdx(edgeIdx)
          call vertices % shelf(vertexIdxs(j)) % addTriangleIdx(freeTriangleIdx)

        end do
        return

      end if

      ! Compute the locations of the second and third vertices in the current triangle and
      ! retrieve their coordinates.
      secondVertexIdx = mod(minVertexIdx + i - 1, nVertices) + 1
      thirdVertexIdx = mod(minVertexIdx + i, nVertices) + 1
      secondVertexCoords = vertices % shelf(vertexIdxs(secondVertexIdx)) % getCoordinates()
      thirdVertexCoords = vertices % shelf(vertexIdxs(thirdVertexIdx)) % getCoordinates()

      ! Set the vertices of the current triangle and compute its centre and normal vector.
      ! Note that because of the vertices' ordering we do not need to check that the normal 
      ! vector points in the correct direction.
      triangleVertexIdxs = vertexIdxs([minVertexIdx, secondVertexIdx, thirdVertexIdx])
      call triangles % shelf(freeTriangleIdx) % setVertices(triangleVertexIdxs)
      call triangles % shelf(freeTriangleIdx) % setAB(secondVertexCoords - minVertexCoords)
      call triangles % shelf(freeTriangleIdx) % setAC(thirdVertexCoords - minVertexCoords)
      call triangles % shelf(freeTriangleIdx) % computeCentre(minVertexCoords, secondVertexCoords, thirdVertexCoords)
      call triangles % shelf(freeTriangleIdx) % computeNormal(minVertexCoords, secondVertexCoords, thirdVertexCoords)
      
      ! Compute the triangle's area, normalise its normal vector and update mesh connectivity information.
      call triangles % shelf(freeTriangleIdx) % computeArea()
      call triangles % shelf(freeTriangleIdx) % normaliseNormal()

      if (i < nVertices - 2) then
        ! Create a new edge.
        lastEdgeIdx = lastEdgeIdx + 1
        call edges % shelf(lastEdgeIdx) % setIdx(lastEdgeIdx)
        call edges % shelf(lastEdgeIdx) % setVertexIdxs(triangleVertexIdxs([1, 3]))
        call vertices % shelf(triangleVertexIdxs(1)) % addEdgeIdx(lastEdgeIdx)
        call vertices % shelf(triangleVertexIdxs(3)) % addEdgeIdx(lastEdgeIdx)

      end if

      ! Update mesh connectivity information.
      do j = 1, 3
        call vertices % shelf(triangleVertexIdxs(j)) % addTriangleIdx(freeTriangleIdx)
        edgeIdx = vertices % findCommonEdgeIdx(triangleVertexIdxs(j), triangleVertexIdxs(mod(j, 3) + 1))
        call triangles % shelf(freeTriangleIdx) % addEdgeIdx(edgeIdx)
        call edges % shelf(edgeIdx) % addTriangleIdx(freeTriangleIdx)

      end do

    end do

  end subroutine split
  
end module face_class