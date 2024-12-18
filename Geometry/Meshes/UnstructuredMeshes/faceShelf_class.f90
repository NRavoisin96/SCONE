module faceShelf_class
  
  use edgeShelf_class,     only : edgeShelf
  use numPrecision
  use genericProcedures,   only : findCommon
  use face_class,          only : face
  use triangleShelf_class, only : triangleShelf
  use vertexShelf_class,   only : vertexShelf
  
  implicit none
  private
  
  !!
  !! Storage space for faces of an OpenFOAM mesh.
  !!
  !! Private members:
  !!   shelf -> Array to store faces.
  !!
  type, public                            :: faceShelf
    private
    type(face), dimension(:), allocatable :: shelf
  contains
    procedure                             :: addEdgeIdxToFace
    procedure                             :: addElementIdxToFace
    procedure                             :: addVertexIdxToFace
    procedure                             :: allocateShelf
    procedure                             :: computeFaceIntersection
    procedure                             :: findCommonEdgeIdx
    procedure                             :: findCommonVertexIdx
    procedure                             :: getFaceArea
    procedure                             :: getFaceCentroid
    procedure                             :: getFaceElementIdxs
    procedure                             :: getFaceIsBoundary
    procedure                             :: getFaceNormal
    procedure                             :: getFaceTriangleIdxs
    procedure                             :: getFaceVertexIdxs
    procedure                             :: getSize
    procedure                             :: initFace
    procedure                             :: kill
    procedure                             :: splitFace
  end type

contains

  !! Subroutine 'addEdgeIdxToFace'
  !!
  !! Basic description:
  !!   Adds the index of an edge to a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the face in the shelf.
  !!   edgeIdx [in] -> Index of the edge in the face.
  !!
  elemental subroutine addEdgeIdxToFace(self, idx, edgeIdx)
    class(faceShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx, edgeIdx

    call self % shelf(idx) % addEdgeIdx(edgeIdx)

  end subroutine addEdgeIdxToFace

  !! Subroutine 'addElementIdxToFace'
  !!
  !! Basic description:
  !!   Adds the index of an element to a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in]        -> Index of the face in the shelf.
  !!   elementIdx [in] -> Index of the element containing the face.
  !!
  elemental subroutine addElementIdxToFace(self, idx, elementIdx)
    class(faceShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx, elementIdx

    call self % shelf(idx) % addElementToFace(elementIdx)

  end subroutine addElementIdxToFace

  !! Subroutine 'addVertexIdxToFace'
  !!
  !! Basic description:
  !!   Adds the index of a vertex to a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in]       -> Index of the face in the shelf.
  !!   vertexIdx [in] -> Index of the vertex in the face.
  !!
  elemental subroutine addVertexIdxToFace(self, idx, vertexIdx)
    class(faceShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx, vertexIdx

    call self % shelf(idx) % addVertexToFace(vertexIdx)

  end subroutine addVertexIdxToFace

  !! Subroutine 'allocateShelf'
  !!
  !! Basic description:
  !!   Allocates memory in the shelf.
  !!
  !! Arguments:
  !!   nFaces [in] -> Number of faces to be included in the shelf.
  !!
  elemental subroutine allocateShelf(self, nFaces)
    class(faceShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: nFaces

    allocate(self % shelf(nFaces))

  end subroutine allocateShelf

  !! Subroutine 'computeFaceIntersection'
  !!
  !! Basic description:
  !!   Computes the intersection of a line segment of origin r, end rEnd and direction u with a
  !!   face in the shelf.
  !!
  !! Arguments:
  !!   idx [in]        -> Index of the face in the shelf.
  !!   r [in]          -> Line segment's origin coordinates.
  !!   rEnd [in]       -> Line segment's end coordinates.
  !!   u [in]          -> Line segment's direction
  !!   vertices [in]   -> A vertexShelf.
  !!   d [out]         -> Distance to intersection with the face.
  !!   edgeIdx [out]   -> Used in case the line segment intersects the face at one of its edges.
  !!   vertexIdx [out] -> Used in case the line segment intersects the face at one of its vertices.
  !!
  pure subroutine computeFaceIntersection(self, idx, r, rEnd, u, vertices, d, edgeIdx, vertexIdx)
    class(faceShelf), intent(in)            :: self
    integer(shortInt), intent(in)           :: idx
    real(defReal), dimension(3), intent(in) :: r, rEnd, u
    type(vertexShelf), intent(in)           :: vertices
    real(defReal), intent(out)              :: d
    integer(shortInt), intent(out)          :: edgeIdx, vertexIdx

    call self % shelf(idx) % computeIntersection(r, rEnd, u, vertices, d, edgeIdx, vertexIdx)

  end subroutine computeFaceIntersection

  !! Function 'findCommonEdgeIdx'
  !!
  !! Basic description:
  !!   Returns the index of the common edge between two faces.
  !!
  !! Arguments:
  !!   firstFaceIdx [in]  -> Index of the first face.
  !!   secondFaceIdx [in] -> Index of the second face.
  !!
  !! Result:
  !!   edgeIdx            -> Index of the common edge between the two faces.
  !!
  elemental function findCommonEdgeIdx(self, firstFaceIdx, secondFaceIdx) result(edgeIdx)
    class(faceShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: firstFaceIdx, secondFaceIdx
    integer(shortInt)                            :: edgeIdx
    integer(shortInt), dimension(:), allocatable :: commonIdxs

    ! Initialise edgeIdx = 0 then find common edge indices between the two faces.
    edgeIdx = 0
    commonIdxs = findCommon(self % shelf(firstFaceIdx) % getEdgeIdxs(), self % shelf(secondFaceIdx) % getEdgeIdxs())

    ! If a common edge has been found update edgeIdx.
    if (size(commonIdxs) > 0) edgeIdx = commonIdxs(1)

  end function findCommonEdgeIdx

  !! Function 'findCommonVertexIdx'
  !!
  !! Basic description:
  !!   Returns the index of the common vertex between a set of faces in the shelf.
  !!
  !! Arguments:
  !!   idxs [in] -> Indices of the faces in the shelf.
  !!
  !! Result:
  !!   vertexIdx -> Index of the common vertex between the faces.
  !!
  pure function findCommonVertexIdx(self, idxs) result(vertexIdx)
    class(faceShelf), intent(in)                 :: self
    integer(shortInt), dimension(:), intent(in)  :: idxs
    integer(shortInt)                            :: vertexIdx, i
    integer(shortInt), dimension(:), allocatable :: commonIdxs

    ! Initialise vertexIdx = 0 then find common vertex indices between the faces.
    vertexIdx = 0
    commonIdxs = self % shelf(idxs(1)) % getVertices()

    do i = 2, size(idxs)
      commonIdxs = findCommon(commonIdxs, self % shelf(idxs(i)) % getVertices())

    end do

    ! If a common vertex has been found update vertexIdx.
    if (size(commonIdxs) > 0) vertexIdx = commonIdxs(1)

  end function findCommonVertexIdx

  !! Function 'getFaceArea'
  !!
  !! Basic description:
  !!   Returns the area of a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the face in the shelf.
  !!
  !! Result:
  !!   area     -> Area of the face's centroid.
  !!
  elemental function getFaceArea(self, idx) result(area)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    real(defReal)                 :: area

    area = self % shelf(idx) % getArea()

  end function getFaceArea

  !! Function 'getFaceCentroid'
  !!
  !! Basic description:
  !!   Returns the 3-D coordinates of the centroid of a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the face in the shelf.
  !!
  !! Result:
  !!   centroid -> 3-D coordinates of the face's centroid.
  !!
  pure function getFaceCentroid(self, idx) result(centroid)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    real(defReal), dimension(3)   :: centroid

    centroid = self % shelf(idx) % getCentroid()

  end function getFaceCentroid

  !! Function 'getFaceElementIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the elements sharing a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in]    -> Index of the face in the shelf.
  !!
  !! Result:
  !!   elementIdxs -> Indices of the elements sharing the face.
  !!
  pure function getFaceElementIdxs(self, idx) result(elementIdxs)
    class(faceShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: elementIdxs

    elementIdxs = self % shelf(idx) % getFaceToElements()

  end function getFaceElementIdxs

  !! Function 'getFaceIsBoundary'
  !!
  !! Basic description:
  !!   Returns .true. if a face in the shelf is a boundary face.
  !!
  !! Arguments:
  !!   idx [in]   -> Index of the face in the shelf.
  !!
  !! Result:
  !!   isBoundary -> .true. if the face in the shelf is a boundary face.
  !!
  elemental function getFaceIsBoundary(self, idx) result(isBoundary)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    logical(defBool)              :: isBoundary

    isBoundary = self % shelf(idx) % getIsBoundary()

  end function getFaceIsBoundary

  !! Function 'getFaceNormal'
  !!
  !! Basic description:
  !!   Returns the signed normal vector of a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the face in the shelf. Can be negative if the normal vector needs to be flipped.
  !!
  !! Result:
  !!   normal   -> 3-D coordinates of the face's signed normal vector.
  !!
  pure function getFaceNormal(self, idx) result(centroid)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    real(defReal), dimension(3)   :: centroid

    centroid = self % shelf(abs(idx)) % getNormal(idx)

  end function getFaceNormal

  !! Function 'getFaceTriangleIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the triangles in a face of the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the face in the shelf.
  !!
  !! Result:
  !!   triangleIdxs -> Indices of the triangles in the face.
  !!
  pure function getFaceTriangleIdxs(self, idx) result(triangleIdxs)
    class(faceShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: triangleIdxs

    triangleIdxs = self % shelf(idx) % getTriangles()

  end function getFaceTriangleIdxs

  !! Function 'getFaceVertexIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in a face of the shelf.
  !!
  !! Arguments:
  !!   idx [in]   -> Index of the face in the shelf.
  !!
  !! Result:
  !!   vertexIdxs -> Indices of the vertices in the face.
  !!
  pure function getFaceVertexIdxs(self, idx) result(vertexIdxs)
    class(faceShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: vertexIdxs

    vertexIdxs = self % shelf(idx) % getVertices()

  end function getFaceVertexIdxs
  
  !! Function 'getSize'
  !!
  !! Basic description:
  !!   Returns the size of the faceShelf.
  !!
  !! Result:
  !!   size -> Size of the faceShelf.
  !!
  elemental function getSize(self) result(nFaces)
    class(faceShelf), intent(in) :: self
    integer(shortInt)            :: nFaces
    
    nFaces = size(self % shelf)
    
  end function getSize

  !! Subroutine 'initFace'
  !!
  !! Basic description:
  !!   Initialises a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in]            -> Index of the face.
  !!   nInternalFaces [in] -> Number of internal faces in the shelf.
  !!
  subroutine initFace(self, idx, nInternalFaces, vertexIdxs, vertices)
    class(faceShelf), intent(inout)             :: self
    integer(shortInt), intent(in)               :: idx, nInternalFaces
    integer(shortInt), dimension(:), intent(in) :: vertexIdxs
    type(vertexShelf), intent(in)               :: vertices
    integer(shortInt)                           :: i
    
    ! Set the face index. If idx > nInternalFaces set this face as a boundary face.
    call self % shelf(idx) % setIdx(idx)
    if (idx > nInternalFaces) call self % shelf(idx) % setBoundaryFace()

    ! Set vertex connectivity information then compute face area and normal vector.
    do i = 1, size(vertexIdxs)
      call self % shelf(idx) % addVertexToFace(vertexIdxs(i))

    end do
    call self % shelf(idx) % computeAreaAndNormal(vertices)

  end subroutine initFace
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(faceShelf), intent(inout) :: self
    
    if (allocated(self % shelf)) deallocate(self % shelf)

  end subroutine kill

  !! Subroutine 'splitFace'
  !!
  !! Basic description:
  !!   Splits a face in the shelf into triangles.
  !!
  !! Arguments:
  !!   idx [in]                -> Index of the face in the shelf.
  !!   edges [inout]           -> An edgeShelf.
  !!   triangles [inout]       -> A triangleShelf.
  !!   vertices [inout]        -> A vertexShelf.
  !!   lastEdgeIdx [inout]     -> Index of the last edge in the edgeShelf.
  !!   lastTriangleIdx [inout] -> Index of the last triangle in the triangleShelf.
  !!
  elemental subroutine splitFace(self, idx, edges, triangles, vertices, lastEdgeIdx, lastTriangleIdx)
    class(faceShelf), intent(inout)    :: self
    integer(shortInt), intent(in)      :: idx
    type(edgeShelf), intent(inout)     :: edges
    type(triangleShelf), intent(inout) :: triangles
    type(vertexShelf), intent(inout)   :: vertices
    integer(shortInt), intent(inout)   :: lastEdgeIdx, lastTriangleIdx

    call self % shelf(idx) % split(edges, triangles, vertices, lastEdgeIdx, lastTriangleIdx)

  end subroutine splitFace
  
end module faceShelf_class