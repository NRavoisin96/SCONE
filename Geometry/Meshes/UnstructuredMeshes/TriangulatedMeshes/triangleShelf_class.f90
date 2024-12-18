module triangleShelf_class
  
  use genericProcedures, only : computeTriangleArea, computeTriangleCentre, computeTriangleNormal, &
                                findCommon, swap
  use numPrecision
  use triangle_class,    only : triangle
  use vertexShelf_class, only : vertexShelf
  
  implicit none
  private
  
  !!
  !! Storage space for triangles in a given OpenFOAM mesh.
  !!
  !! Private members:
  !!   shelf -> An array of triangles.
  !!
  type, public                                :: triangleShelf
    private
    type(triangle), dimension(:), allocatable :: shelf
  contains
    procedure                                 :: addEdgeIdxToTriangle
    procedure                                 :: addTetrahedronIdxToTriangle
    procedure                                 :: allocateShelf
    procedure                                 :: buildTriangle
    procedure                                 :: computeTriangleIntersection
    procedure                                 :: findCommonEdgeIdx
    procedure                                 :: findCommonVertexIdx
    procedure                                 :: getTriangleCentroid
    procedure                                 :: getTriangleEdgeIdxs
    procedure                                 :: getTriangleHasTetrahedra
    procedure                                 :: getTriangleIsBoundary
    procedure                                 :: getTriangleNormal
    procedure                                 :: getTriangleTetrahedronIdxs
    procedure                                 :: getTriangleVertexIdxs
    procedure                                 :: initTriangle
    procedure                                 :: kill
  end type triangleShelf

contains

  !! Subroutine 'addEdgeIdxToTriangle'
  !!
  !! Basic description:
  !!   Adds the index of an edge to a triangle in the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the triangle in the shelf.
  !!   edgeIdx [in] -> Index of the edge in the triangle.
  !!
  elemental subroutine addEdgeIdxToTriangle(self, idx, edgeIdx)
    class(triangleShelf), intent(inout) :: self
    integer(shortInt), intent(in)       :: idx, edgeIdx

    call self % shelf(idx) % addEdgeIdx(edgeIdx)

  end subroutine addEdgeIdxToTriangle

  !! Subroutine 'addTetrahedronIdxToTriangle'
  !!
  !! Basic description:
  !!   Adds the index of a tetrahedron to a triangle in the shelf.
  !!
  !! Arguments:
  !!   idx [in]            -> Index of the triangle in the shelf.
  !!   tetrahedronIdx [in] -> Index of the tetrahedron containing the triangle.
  !!
  elemental subroutine addTetrahedronIdxToTriangle(self, idx, tetrahedronIdx)
    class(triangleShelf), intent(inout) :: self
    integer(shortInt), intent(in)       :: idx, tetrahedronIdx

    call self % shelf(idx) % addTetrahedronIdx(tetrahedronIdx)

  end subroutine addTetrahedronIdxToTriangle

  !! Subroutine 'allocateShelf'
  !!
  !! Basic description:
  !!   Allocates memory in the shelf.
  !!
  !! Arguments:
  !!   nTriangles [in] -> Number of triangles to be included in the shelf.
  !!
  elemental subroutine allocateShelf(self, nTriangles)
    class(triangleShelf), intent(inout) :: self
    integer(shortInt), intent(in)       :: nTriangles

    allocate(self % shelf(nTriangles))

  end subroutine allocateShelf

  !! Subroutine 'buildTriangle'
  !!
  !! Basic description:
  !!   Computes the centre, normal vector and area of a triangle then initialises this triangle in the shelf.
  !!
  !! Arguments:
  !!   idx [in]                     -> Index of the triangle in the shelf.
  !!   vertexIdxs [in]              -> Indices of the vertices of the triangle.
  !!   faceIdx [in]                 -> Index of the face from which the triangle originates (if the face is 
  !!                                   split).
  !!   isBoundary [in]              -> .true. if the triangle is a boundary triangle.
  !!   coordsArray [in]             -> A 3x3 defReal array containing the 3-D coordinates of the vertices of 
  !!                                   the triangle.
  !!   testCentroid [in] [optional] -> 3-D coordinates of a centroid used to test the orientation of the new
  !!                                   triangle's normal vector.
  !!
  pure subroutine buildTriangle(self, idx, vertexIdxs, faceIdx, isBoundary, coordsArray, testCentroid)
    class(triangleShelf), intent(inout)               :: self
    integer(shortInt), intent(in)                     :: idx, faceIdx
    integer(shortInt), dimension(3), intent(inout)    :: vertexIdxs
    logical(defBool), intent(in)                      :: isBoundary
    real(defReal), dimension(3, 3), intent(inout)     :: coordsArray
    real(defReal), dimension(3), intent(in), optional :: testCentroid
    real(defReal), dimension(3)                       :: centroid, normal, tempCoords
    real(defReal)                                     :: area

    ! Compute centroid, normal vector and area of the triangle then normalise its normal vector.
    centroid = computeTriangleCentre(coordsArray)
    normal = computeTriangleNormal(coordsArray)
    area = computeTriangleArea(normal)
    normal = normal / norm2(normal)

    if (present(testCentroid)) then
      ! Check that the triangle's normal vector points in the correct direction (outward). If not,
      ! swap two of the new triangle's vertices and negate its normal vector.
      if (dot_product(centroid - testCentroid, normal) <= ZERO) then
        call swap(vertexIdxs, 1, 2)
        tempCoords = coordsArray(1, :)
        coordsArray(1, :) = coordsArray(2, :)
        coordsArray(2, :) = tempCoords
        normal = -normal

      end if

    end if

    ! Set everything.
    call self % initTriangle(idx, vertexIdxs, faceIdx, isBoundary, coordsArray(2, :) - coordsArray(1, :), &
                             coordsArray(3, :) - coordsArray(1, :), centroid, normal, area)

  end subroutine buildTriangle

  !! Subroutine 'computeTriangleIntersection'
  !!
  !! Basic description:
  !!   Computes the intersection of a line segment of origin r, end rEnd and direction u with a
  !!   face in the shelf.
  !!
  !! Arguments:
  !!   idx [in]               -> Index of the face in the shelf.
  !!   r [in]                 -> Line segment's origin coordinates.
  !!   rEnd [in]              -> Line segment's end coordinates.
  !!   u [in]                 -> Line segment's direction
  !!   firstVertexCoords [in] -> 3-D coordinates of the first triangle vertex.
  !!   vertices [in]          -> A vertexShelf.
  !!   d [out]                -> Distance to intersection with the face.
  !!   edgeIdx [out]          -> Used in case the line segment intersects the face at one of its edges.
  !!   vertexIdx [out]        -> Used in case the line segment intersects the face at one of its vertices.
  !!
  pure subroutine computeTriangleIntersection(self, idx, r, rEnd, u, firstVertexCoords, vertices, d, edgeIdx, vertexIdx)
    class(triangleShelf), intent(in)        :: self
    integer(shortInt), intent(in)           :: idx
    real(defReal), dimension(3), intent(in) :: r, rEnd, u, firstVertexCoords
    type(vertexShelf), intent(in)           :: vertices
    real(defReal), intent(out)              :: d
    integer(shortInt), intent(out)          :: edgeIdx, vertexIdx

    call self % shelf(idx) % computeIntersection(vertices, r, rEnd, u, firstVertexCoords, d, edgeIdx, vertexIdx)

  end subroutine computeTriangleIntersection

  !! Function 'findCommonEdgeIdx'
  !!
  !! Basic description:
  !!   Returns the index of the common edge between two triangles.
  !!
  !! Arguments:
  !!   firstTriangleIdx [in]  -> Index of the first triangle.
  !!   secondTriangleIdx [in] -> Index of the second triangle.
  !!
  !! Result:
  !!   edgeIdx                -> Index of the common edge between the two triangles.
  !!
  elemental function findCommonEdgeIdx(self, firstTriangleIdx, secondTriangleIdx) result(edgeIdx)
    class(triangleShelf), intent(in) :: self
    integer(shortInt), intent(in)    :: firstTriangleIdx, secondTriangleIdx
    integer(shortInt)                :: edgeIdx
    integer(shortInt), dimension(:), allocatable :: commonIdxs

    ! Initialise edgeIdx = 0 then find common edge indices between the two triangles.
    edgeIdx = 0
    commonIdxs = findCommon(self % shelf(firstTriangleIdx) % getEdgeIdxs(), self % shelf(secondTriangleIdx) % getEdgeIdxs())
    
    ! If a common edge has been found update edgeIdx.
    if (size(commonIdxs) > 0) edgeIdx = commonIdxs(1)

  end function findCommonEdgeIdx

  !! Function 'findCommonVertexIdx'
  !!
  !! Basic description:
  !!   Returns the index of the common vertex between a set of triangles in the shelf.
  !!
  !! Arguments:
  !!   idxs [in] -> Indices of the triangles in the shelf.
  !!
  !! Result:
  !!   vertexIdx -> Index of the common vertex between the triangles.
  !!
  pure function findCommonVertexIdx(self, idxs) result(vertexIdx)
    class(triangleShelf), intent(in)             :: self
    integer(shortInt), dimension(:), intent(in)  :: idxs
    integer(shortInt)                            :: vertexIdx, i
    integer(shortInt), dimension(:), allocatable :: commonIdxs

    ! Initialise vertexIdx = 0 then find common vertex indices between the triangles.
    vertexIdx = 0
    commonIdxs = self % shelf(idxs(1)) % getVertices()

    do i = 2, size(idxs)
      commonIdxs = findCommon(commonIdxs, self % shelf(idxs(i)) % getVertices())

    end do

    ! If a common vertex has been found update vertexIdx.
    if (size(commonIdxs) > 0) vertexIdx = commonIdxs(1)

  end function findCommonVertexIdx

  !! Function 'getTriangleCentroid'
  !!
  !! Basic description:
  !!   Returns the centre of a triangle in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the triangle in the shelf.
  !!
  !! Result:
  !!   centroid -> 3-D coordinates of the triangle's centre.
  !!
  pure function getTriangleCentroid(self, idx) result(centroid)
    class(triangleShelf), intent(in) :: self
    integer(shortInt), intent(in)    :: idx
    real(defReal), dimension(3)      :: centroid

    centroid = self % shelf(idx) % getCentre()

  end function getTriangleCentroid

  !! Function 'getTriangleEdgeIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the edges of a triangle in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the triangle in the shelf.
  !!
  !! Result:
  !!   edgeIdxs -> Indices of the edges of the triangle.
  !!
  pure function getTriangleEdgeIdxs(self, idx) result(edgeIdxs)
    class(triangleShelf), intent(in) :: self
    integer(shortInt), intent(in)    :: idx
    integer(shortInt), dimension(3)  :: edgeIdxs

    edgeIdxs = self % shelf(idx) % getEdgeIdxs()

  end function getTriangleEdgeIdxs

  !! Function 'getTriangleHasTetrahedra'
  !!
  !! Basic description:
  !!   Returns .true. if a triangle in the shelf is associated with tetrahedra.
  !!
  !! Arguments:
  !!   idx [in]      -> Index of the triangle in the shelf.
  !!
  !! Result:
  !!   hasTetrahedra -> .true. if the triangle is associated with tetrahedra.
  !!
  elemental function getTriangleHasTetrahedra(self, idx) result(hasTetrahedra)
    class(triangleShelf), intent(in) :: self
    integer(shortInt), intent(in)    :: idx
    logical(defBool)                 :: hasTetrahedra

    hasTetrahedra = self % shelf(idx) % hasTetrahedra()

  end function getTriangleHasTetrahedra

  !! Function 'getTriangleIsBoundary'
  !!
  !! Basic description:
  !!   Returns .true. if a triangle in the shelf is a boundary triangle.
  !!
  !! Arguments:
  !!   idx [in]   -> Index of the triangle in the shelf.
  !!
  !! Result:
  !!   isBoundary -> .true. if the triangle is a boundary triangle.
  !!
  elemental function getTriangleIsBoundary(self, idx) result(isBoundary)
    class(triangleShelf), intent(in) :: self
    integer(shortInt), intent(in)    :: idx
    logical(defBool)                 :: isBoundary

    isBoundary = self % shelf(idx) % getIsBoundary()

  end function getTriangleIsBoundary

  !! Function 'getTriangleNormal'
  !!
  !! Basic description:
  !!   Returns the signed normal of a triangle in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the triangle in the shelf.
  !!
  !! Result:
  !!   normal   -> 3-D coordinates of the triangle's normal vector.
  !!
  pure function getTriangleNormal(self, idx) result(normal)
    class(triangleShelf), intent(in) :: self
    integer(shortInt), intent(in)    :: idx
    real(defReal), dimension(3)      :: normal

    normal = self % shelf(abs(idx)) % getNormal(idx)

  end function getTriangleNormal

  !! Function 'getTriangleTetrahedronIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the tetrahedra containing a triangle in the shelf.
  !!
  !! Arguments:
  !!   idx [in]        -> Index of the triangle in the shelf.
  !!
  !! Result:
  !!   tetrahedronIdxs -> Indices of the tetrahedra containing the triangle.
  !!
  pure function getTriangleTetrahedronIdxs(self, idx) result(tetrahedronIdxs)
    class(triangleShelf), intent(in)             :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: tetrahedronIdxs

    tetrahedronIdxs = self % shelf(idx) % getTetrahedra()

  end function getTriangleTetrahedronIdxs

  !! Function 'getTriangleVertexIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices of a triangle in the shelf.
  !!
  !! Arguments:
  !!   idx [in]   -> Index of the triangle in the shelf.
  !!
  !! Result:
  !!   vertexIdxs -> Indices of the vertices of the triangle.
  !!
  pure function getTriangleVertexIdxs(self, idx) result(vertexIdxs)
    class(triangleShelf), intent(in) :: self
    integer(shortInt), intent(in)    :: idx
    integer(shortInt), dimension(3)  :: vertexIdxs

    vertexIdxs = self % shelf(idx) % getVertices()

  end function getTriangleVertexIdxs

  !! Subroutine 'initTriangle'
  !!
  !! Basic description:
  !!   Initialises a triangle in the shelf.
  !!
  !! Arguments:
  !!   idx [in]        -> Index of the triangle in the shelf.
  !!   vertexIdxs [in] -> Indices of the vertices of the triangle.
  !!   faceIdx [in]    -> Index of the face from which the triangle originates (if the face is split).
  !!   isBoundary [in] -> .true. if the triangle is a boundary triangle.
  !!   AB [in]         -> 3-D coordinates of the first edge vector of the triangle.
  !!   AC [in]         -> 3-D coordinates of the second edge vector of the triangle.
  !!   centroid [in]   -> 3-D coordinates of the centre of the triangle.
  !!   normal [in]     -> 3-D coordinates of the normal vector of the triangle.
  !!   area [in]       -> Area of the triangle.
  !!
  pure subroutine initTriangle(self, idx, vertexIdxs, faceIdx, isBoundary, AB, AC, centroid, normal, area)
    class(triangleShelf), intent(inout)         :: self
    integer(shortInt), intent(in)               :: idx, faceIdx
    integer(shortInt), dimension(3), intent(in) :: vertexIdxs
    logical(defBool), intent(in)                :: isBoundary
    real(defReal), dimension(3), intent(in)     :: AB, AC, centroid, normal
    real(defReal), intent(in)                   :: area

    call self % shelf(idx) % setIdx(idx)
    call self % shelf(idx) % setVertices(vertexIdxs)
    call self % shelf(idx) % setFace(faceIdx)
    if (isBoundary) call self % shelf(idx) % setIsBoundary()
    call self % shelf(idx) % setAB(AB)
    call self % shelf(idx) % setAC(AC)
    call self % shelf(idx) % setCentre(centroid)
    call self % shelf(idx) % setNormal(normal)
    call self % shelf(idx) % setArea(area)

  end subroutine initTriangle
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(triangleShelf), intent(inout) :: self
    
    if (allocated(self % shelf)) deallocate(self % shelf)

  end subroutine kill
  
end module triangleShelf_class