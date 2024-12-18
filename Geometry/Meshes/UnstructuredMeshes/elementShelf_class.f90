module elementShelf_class
  
  use edgeShelf_class,     only : edgeShelf
  use element_class,       only : element
  use faceShelf_class,     only : faceShelf
  use numPrecision
  use pyramidShelf_class,  only : pyramidShelf
  use triangleShelf_class, only : triangleShelf
  use vertexShelf_class,   only : vertexShelf
  
  implicit none
  private
  
  !!
  !! Storage space for elements in a given OpenFOAM mesh.
  !!
  !! Private members:
  !!   shelf -> Array to store elements.
  !!
  type, public                               :: elementShelf
    private
    type(element), dimension(:), allocatable :: shelf
  contains
    procedure                                :: addEdgeIdxToElement
    procedure                                :: addFaceIdxToElement
    procedure                                :: addVertexIdxToElement
    procedure                                :: allocateShelf
    procedure                                :: computeFaceIntersection
    procedure                                :: computePotentialFaceIdxs
    procedure                                :: getElementCentroid
    procedure                                :: getElementFaceIdxs
    procedure                                :: getElementVertexIdxs
    procedure                                :: getElementVolume
    procedure                                :: initElement
    procedure                                :: isConvex
    procedure                                :: kill
    procedure                                :: splitElement
    procedure                                :: testForInclusion
  end type elementShelf

contains

  !! Subroutine 'addEdgeIdxToElement'
  !!
  !! Basic description:
  !!   Adds the index of an edge to an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the element in the shelf.
  !!   edgeIdx [in] -> Index of the edge in the element.
  !!
  elemental subroutine addEdgeIdxToElement(self, idx, edgeIdx)
    class(elementShelf), intent(inout) :: self
    integer(shortInt), intent(in)      :: idx, edgeIdx

    call self % shelf(idx) % addEdgeIdx(edgeIdx)

  end subroutine addEdgeIdxToElement

  !! Subroutine 'addFaceIdxToElement'
  !!
  !! Basic description:
  !!   Adds the index of a face to an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the element in the shelf.
  !!   faceIdx [in] -> Index of the face in the element.
  !!
  elemental subroutine addFaceIdxToElement(self, idx, faceIdx)
    class(elementShelf), intent(inout) :: self
    integer(shortInt), intent(in)      :: idx, faceIdx

    call self % shelf(idx) % addFaceToElement(faceIdx)

  end subroutine addFaceIdxToElement

  !! Subroutine 'addVertexIdxToElement'
  !!
  !! Basic description:
  !!   Adds the index of a vertex to an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in]       -> Index of the element in the shelf.
  !!   vertexIdx [in] -> Index of the vertex in the element.
  !!
  elemental subroutine addVertexIdxToElement(self, idx, vertexIdx)
    class(elementShelf), intent(inout) :: self
    integer(shortInt), intent(in)      :: idx, vertexIdx

    call self % shelf(idx) % addVertexToElement(vertexIdx)

  end subroutine addVertexIdxToElement

  !! Subroutine 'allocateShelf'
  !!
  !! Basic description:
  !!   Allocates memory in the shelf.
  !!
  !! Arguments:
  !!   nElements [in] -> Number of elements to be included in the shelf.
  !!
  elemental subroutine allocateShelf(self, nElements)
    class(elementShelf), intent(inout) :: self
    integer(shortInt), intent(in)      :: nElements

    allocate(self % shelf(nElements))

  end subroutine allocateShelf

  !! Subroutine 'computeFaceIntersection'
  !!
  !! Basic description:
  !!   Computes the intersection of a line segment with the faces of an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in]                 -> Index of the element in the shelf.
  !!   r [in]                   -> Line segment's origin coordinates.
  !!   rEnd [in]                -> Line segment's end coordinates.
  !!   potentialFaceIdxs [in]   -> Indices of the element faces potentially intersected by the line segment.
  !!   faces [in]               -> A faceShelf.
  !!   intersectedFaceIdx [out] -> Index of the intersected element face.
  !!   lambda [out]             -> Fraction of the line segment to intersection.
  !!
  pure subroutine computeFaceIntersection(self, idx, r, rEnd, potentialFaceIdxs, faces, intersectedFaceIdx, lambda)
    class(elementShelf), intent(in)             :: self
    integer(shortInt), intent(in)               :: idx
    real(defReal), dimension(3), intent(in)     :: r, rEnd
    integer(shortInt), dimension(:), intent(in) :: potentialFaceIdxs
    type(faceShelf), intent(in)                 :: faces
    integer(shortInt), intent(out)              :: intersectedFaceIdx
    real(defReal), intent(out)                  :: lambda

    call self % shelf(idx) % computeIntersectedFace(r, rEnd, potentialFaceIdxs, intersectedFaceIdx, lambda, faces)

  end subroutine computeFaceIntersection

  !! Function 'computePotentialFaceIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the potentially intersected faces of an element in the shelf by a line segment.
  !!
  !! Argument:
  !!   idx [in]          -> Index of the element in the shelf.
  !!   rEnd [in]         -> Line segment's end coordinates.
  !!   faces [in]        -> A faceShelf.
  !!
  !! Result:
  !!   potentialFaceIdxs -> Indices of the element's faces potentially intersected by the line segment.
  !!
  pure function computePotentialFaceIdxs(self, idx, rEnd, faces) result(potentialFaceIdxs)
    class(elementShelf), intent(in)              :: self
    integer(shortInt), intent(in)                :: idx
    real(defReal), dimension(3), intent(in)      :: rEnd
    type(faceShelf), intent(in)                  :: faces
    integer(shortInt), dimension(:), allocatable :: potentialFaceIdxs

    potentialFaceIdxs = self % shelf(idx) % computePotentialFaces(rEnd, faces)

  end function computePotentialFaceIdxs

  !! Function 'getElementCentroid'
  !!
  !! Basic description:
  !!   Returns the centroid of an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element in the shelf.
  !!
  !! Result:
  !!   centroid -> 3-D coordinates of the element's centroid.
  !!
  pure function getElementCentroid(self, idx) result(centroid)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    real(defReal), dimension(3)     :: centroid

    centroid = self % shelf(idx) % getCentroid()

  end function getElementCentroid

  !! Function 'getElementFaceIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the faces in an element of the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element in the shelf.
  !!
  !! Result:
  !!   faceIdxs -> Indices of the faces in the element.
  !!
  pure function getElementFaceIdxs(self, idx) result(faceIdxs)
    class(elementShelf), intent(in)              :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: faceIdxs

    faceIdxs = self % shelf(idx) % getFaces()

  end function getElementFaceIdxs

  !! Function 'getElementVertexIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in an element of the shelf.
  !!
  !! Arguments:
  !!   idx [in]   -> Index of the element in the shelf.
  !!
  !! Result:
  !!   vertexIdxs -> Indices of the vertices in the element.
  !!
  pure function getElementVertexIdxs(self, idx) result(vertexIdxs)
    class(elementShelf), intent(in)              :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: vertexIdxs

    vertexIdxs = self % shelf(idx) % getVertices()

  end function getElementVertexIdxs

  !! Function 'getElementVolume'
  !!
  !! Basic description:
  !!   Returns the volume of an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element in the shelf.
  !!
  !! Result:
  !!   volume   -> Volume of the element.
  !!
  elemental function getElementVolume(self, idx) result(volume)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    real(defReal)                   :: volume

    volume = self % shelf(idx) % getVolume()

  end function getElementVolume

  !! Subroutine 'initElement'
  !!
  !! Basic description:
  !!   Initialises an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in]      -> Index of the element in the shelf.
  !!   faces [in]    -> A faceShelf.
  !!   vertices [in] -> A vertexShelf.
  !!
  subroutine initElement(self, idx, faces, vertices)
    class(elementShelf), intent(inout) :: self
    integer(shortInt), intent(in)      :: idx
    type(faceShelf), intent(in)        :: faces
    type(vertexShelf), intent(in)      :: vertices

    call self % shelf(idx) % setIdx(idx)
    call self % shelf(idx) % computeVolumeAndCentroid(vertices, faces)

  end subroutine initElement

  !! Function 'isConvex'
  !!
  !! Basic description:
  !!   Returns .true. if an element in the shelf is convex.
  !!
  !! Arguments:
  !!   idx [in]      -> Index of the element in the shelf.
  !!   faces [in]    -> A faceShelf.
  !!   vertices [in] -> A vertexShelf.
  !!
  !! Result:
  !!   isIt          -> .true. if the element is convex.
  !!
  elemental function isConvex(self, idx, faces, vertices) result(isIt)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    type(faceShelf), intent(in)     :: faces
    type(vertexShelf), intent(in)   :: vertices
    logical(defBool)                :: isIt

    isIt = self % shelf(idx) % isConvex(vertices, faces)

  end function isConvex
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(elementShelf), intent(inout) :: self
    
    if (allocated(self % shelf)) deallocate(self % shelf)

  end subroutine kill

  !! Subroutine 'splitElement'
  !!
  !! Basic description:
  !!   Splits an element in the shelf into pyramids.
  !!
  !! Arguments:
  !!   idx [in]                -> Index of the element in the shelf.
  !!   faces [in]              -> A faceShelf.
  !!   edges [inout]           -> An edgeShelf.
  !!   vertices [inout]        -> A vertexShelf.
  !!   triangles [inout]       -> A triangleShelf.
  !!   pyramids [inout]        -> A pyramidShelf.
  !!   lastEdgeIdx [inout]     -> Index of the last edge in the edgeShelf.
  !!   lastTriangleIdx [inout] -> Index of the last triangle in the triangleShelf.
  !!   lastPyramidIdx [inout]  -> Index of the last pyramid in the pyramidShelf.
  !!   lastVertexIdx [in]      -> Index of the last vertex in the vertexIdx.
  !!   
  elemental subroutine splitElement(self, idx, faces, edges, vertices, triangles, pyramids, &
                                    lastEdgeIdx, lastTriangleIdx, lastPyramidIdx, lastVertexIdx)
    class(elementShelf), intent(inout) :: self
    integer(shortInt), intent(in)      :: idx
    type(faceShelf), intent(in)        :: faces
    type(edgeShelf), intent(inout)     :: edges
    type(vertexShelf), intent(inout)   :: vertices
    type(triangleShelf), intent(inout) :: triangles
    type(pyramidShelf), intent(inout)  :: pyramids
    integer(shortInt), intent(inout)   :: lastEdgeIdx, lastTriangleIdx, lastPyramidIdx
    integer(shortInt), intent(in)      :: lastVertexIdx

    call self % shelf(idx) % split(faces, edges, vertices, triangles, pyramids, lastEdgeIdx, &
                                   lastTriangleIdx, lastPyramidIdx, lastVertexIdx)

  end subroutine splitElement

  !! Subroutine 'testForInclusion'
  !!
  !! Basic description:
  !!   Tests whether a given element in the shelf contains a point.
  !!
  !! Arguments:
  !!   idx [in]              -> Index of the element in the shelf.
  !!   r [in]                -> 3-D coordinates of the point.
  !!   faces [in]            -> A faceShelf.
  !!   failedFaceIdx [out]   -> Index of the first element's face for which the test fails.
  !!   surfTolFaceIdxs [out] -> Indices of the element's faces on which the point lies.
  !!
  pure subroutine testForInclusion(self, idx, r, faces, failedFaceIdx, surfTolFaceIdxs)
    class(elementShelf), intent(in)                           :: self
    integer(shortInt), intent(in)                             :: idx
    real(defReal), dimension(3), intent(in)                   :: r
    type(faceShelf), intent(in)                               :: faces
    integer(shortInt), intent(out)                            :: failedFaceIdx
    integer(shortInt), dimension(:), allocatable, intent(out) :: surfTolFaceIdxs

    call self % shelf(idx) % testForInclusion(faces, r, failedFaceIdx, surfTolFaceIdxs)

  end subroutine testForInclusion
  
end module elementShelf_class