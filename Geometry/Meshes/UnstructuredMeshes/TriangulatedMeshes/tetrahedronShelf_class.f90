module tetrahedronShelf_class
  
  use genericProcedures,   only : computeTetrahedronCentre
  use numPrecision
  use tetrahedron_class,   only : tetrahedron
  use triangleShelf_class, only : triangleShelf
  
  implicit none
  private
  
  !!
  !! Storage space for tetrahedra in a given OpenFOAM mesh.
  !!
  !! Private members:
  !!   shelf -> An array of tetrahedra.
  !!
  type, public                                   :: tetrahedronShelf
    private
    type(tetrahedron), dimension(:), allocatable :: shelf
  contains
    procedure                                    :: addEdgeIdxToTetrahedron
    procedure                                    :: addTriangleIdxToTetrahedron
    procedure                                    :: allocateShelf
    procedure                                    :: buildTetrahedron
    procedure                                    :: computePotentialTriangleIdxs
    procedure                                    :: computeTriangleIntersection
    procedure                                    :: getTetrahedronCentroid
    procedure                                    :: getTetrahedronEdgeIdxs
    procedure                                    :: getTetrahedronElementIdx
    procedure                                    :: getTetrahedronTriangleIdxs
    procedure                                    :: initTetrahedron
    procedure                                    :: kill
    procedure                                    :: testForInclusion
  end type tetrahedronShelf

contains

  !! Subroutine 'addEdgeIdxToTetrahedron'
  !!
  !! Basic description:
  !!   Adds the index of an edge to a tetrahedron in the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the tetrahedron in the shelf.
  !!   edgeIdx [in] -> Index of the edge in the tetrahedron.
  !!
  elemental subroutine addEdgeIdxToTetrahedron(self, idx, edgeIdx)
    class(tetrahedronShelf), intent(inout) :: self
    integer(shortInt), intent(in)          :: idx, edgeIdx

    call self % shelf(idx) % addEdgeIdx(edgeIdx)

  end subroutine addEdgeIdxToTetrahedron

  !! Subroutine 'addTriangleIdxToTetrahedron'
  !!
  !! Basic description:
  !!   Adds the index of a triangle to a tetrahedron in the shelf.
  !!
  !! Arguments:
  !!   idx [in]         -> Index of the tetrahedron in the shelf.
  !!   triangleIdx [in] -> Index of the triangle in the tetrahedron.
  !!
  elemental subroutine addTriangleIdxToTetrahedron(self, idx, triangleIdx)
    class(tetrahedronShelf), intent(inout) :: self
    integer(shortInt), intent(in)          :: idx, triangleIdx

    call self % shelf(idx) % addTriangle(triangleIdx)

  end subroutine addTriangleIdxToTetrahedron

  !! Subroutine 'allocateShelf'
  !!
  !! Basic description:
  !!   Allocates memory in the shelf.
  !!
  !! Arguments:
  !!   nTetrahedra [in] -> Number of tetrahedra to be included in the shelf.
  !!
  elemental subroutine allocateShelf(self, nTetrahedra)
    class(tetrahedronShelf), intent(inout) :: self
    integer(shortInt), intent(in)          :: nTetrahedra

    allocate(self % shelf(nTetrahedra))

  end subroutine allocateShelf

  !! Subroutine 'buildTetrahedron'
  !!
  !! Basic description:
  !!   Computes the centroid of a tetrahedron then initialises this tetrahedron in the shelf.
  !!
  !! Arguments:
  !!   idx [in]         -> Index of the tetrahedron in the shelf.
  !!   elementIdx [in]  -> Index of the element from which the tetrahedron originates (if the face is split).
  !!   vertexIdxs [in]  -> Indices of the vertices of the tetrahedron.
  !!   coordsArray [in] -> A 4x3 defReal array containing the 3-D coordinates of the vertices of the tetrahedron.
  !!
  pure subroutine buildTetrahedron(self, idx, elementIdx, vertexIdxs, coordsArray)
    class(tetrahedronShelf), intent(inout)      :: self
    integer(shortInt), intent(in)               :: idx, elementIdx
    integer(shortInt), dimension(4), intent(in) :: vertexIdxs
    real(defReal), dimension(4, 3), intent(in)  :: coordsArray
    real(defReal), dimension(3)                 :: centroid

    ! Compute tetrahedron centre.
    centroid = computeTetrahedronCentre(coordsArray)

    ! Set everything.
    call self % initTetrahedron(idx, elementIdx, vertexIdxs, centroid)

  end subroutine buildTetrahedron

  !! Function 'computePotentialTriangleIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the potentially intersected triangles of a tetrahedron in the shelf by a line segment.
  !!
  !! Argument:
  !!   idx [in]              -> Index of the element in the shelf.
  !!   rEnd [in]             -> Line segment's end coordinates.
  !!   triangles [in]        -> A triangleShelf.
  !!
  !! Result:
  !!   potentialTriangleIdxs -> Indices of the tetrahedron's triangles potentially intersected by the line segment.
  !!
  pure function computePotentialTriangleIdxs(self, idx, rEnd, triangles) result(potentialTriangleIdxs)
    class(tetrahedronShelf), intent(in)          :: self
    integer(shortInt), intent(in)                :: idx
    real(defReal), dimension(3), intent(in)      :: rEnd
    type(triangleShelf), intent(in)              :: triangles
    integer(shortInt), dimension(:), allocatable :: potentialTriangleIdxs

    potentialTriangleIdxs = self % shelf(idx) % computePotentialTriangles(rEnd, triangles)

  end function computePotentialTriangleIdxs

  !! Subroutine 'computeTriangleIntersection'
  !!
  !! Basic description:
  !!   Computes the intersection of a line segment with the triangles of a tetrahedron in the shelf.
  !!
  !! Arguments:
  !!   idx [in]                     -> Index of the tetrahedron in the shelf.
  !!   r [in]                       -> Line segment's origin coordinates.
  !!   rEnd [in]                    -> Line segment's end coordinates.
  !!   potentialTriangleIdxs [in]   -> Indices of the tetrahedron triangles potentially intersected by the line segment.
  !!   triangles [in]               -> A triangleShelf.
  !!   intersectedTriangleIdx [out] -> Index of the intersected tetrahedron triangle.
  !!   lambda [out]                 -> Fraction of the line segment to intersection.
  !!
  pure subroutine computeTriangleIntersection(self, idx, r, rEnd, potentialTriangleIdxs, triangles, intersectedTriangleIdx, lambda)
    class(tetrahedronShelf), intent(in)         :: self
    integer(shortInt), intent(in)               :: idx
    real(defReal), dimension(3), intent(in)     :: r, rEnd
    integer(shortInt), dimension(:), intent(in) :: potentialTriangleIdxs
    type(triangleShelf), intent(in)             :: triangles
    integer(shortInt), intent(out)              :: intersectedTriangleIdx
    real(defReal), intent(out)                  :: lambda

    call self % shelf(idx) % computeIntersectedTriangle(r, rEnd, potentialTriangleIdxs, intersectedTriangleIdx, lambda, triangles)

  end subroutine computeTriangleIntersection

  !! Function 'getTetrahedronCentroid'
  !!
  !! Basic description:
  !!   Returns the centroid of a tetrahedron in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the tetrahedron in the shelf.
  !!
  !! Result:
  !!   centroid -> 3-D coordinates of the tetrahedron's centroid.
  !!
  pure function getTetrahedronCentroid(self, idx) result(centroid)
    class(tetrahedronShelf), intent(in) :: self
    integer(shortInt), intent(in)       :: idx
    real(defReal), dimension(3)         :: centroid

    centroid = self % shelf(idx) % getCentroid()

  end function getTetrahedronCentroid

  !! Function 'getTetrahedronEdgeIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the edges in a tetrahedron of the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the tetrahedron in the shelf.
  !!
  !! Result:
  !!   edgeIdxs -> Indices of the edges in the tetrahedron.
  !!
  pure function getTetrahedronEdgeIdxs(self, idx) result(edgeIdxs)
    class(tetrahedronShelf), intent(in) :: self
    integer(shortInt), intent(in)       :: idx
    integer(shortInt), dimension(6)     :: edgeIdxs

    edgeIdxs = self % shelf(idx) % getEdgeIdxs()

  end function getTetrahedronEdgeIdxs

  !! Function 'getTetrahedronElementIdx'
  !!
  !! Basic description:
  !!   Returns the index of the parent element of a tetrahedron of the shelf.
  !!
  !! Arguments:
  !!   idx [in]   -> Index of the tetrahedron in the shelf.
  !!
  !! Result:
  !!   elementIdx -> Index of the parent element of the tetrahedron.
  !!
  elemental function getTetrahedronElementIdx(self, idx) result(elementIdx)
    class(tetrahedronShelf), intent(in) :: self
    integer(shortInt), intent(in)       :: idx
    integer(shortInt)                   :: elementIdx

    elementIdx = self % shelf(idx) % getElement()

  end function getTetrahedronElementIdx

  !! Function 'getTetrahedronTriangleIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the triangles in a tetrahedron of the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the tetrahedron in the shelf.
  !!
  !! Result:
  !!   triangleIdxs -> Indices of the triangles in the tetrahedron.
  !!
  pure function getTetrahedronTriangleIdxs(self, idx) result(triangleIdxs)
    class(tetrahedronShelf), intent(in) :: self
    integer(shortInt), intent(in)       :: idx
    integer(shortInt), dimension(4)     :: triangleIdxs

    triangleIdxs = self % shelf(idx) % getTriangles()

  end function getTetrahedronTriangleIdxs

  !! Subroutine 'initTetrahedron'
  !!
  !! Basic description:
  !!   Initialises a tetrahedron in the shelf.
  !!
  !! Arguments:
  !!   idx [in]        -> Index of the tetrahedron in the shelf.
  !!   elementIdx [in] -> Index of the element from which the tetrahedron originates (if the element is split).
  !!   vertexIdxs [in] -> Indices of the vertices of the tetrahedron.
  !!   centroid [in]   -> 3-D coordinates of the centroid of the tetrahedron.
  !!
  pure subroutine initTetrahedron(self, idx, elementIdx, vertexIdxs, centroid)
    class(tetrahedronShelf), intent(inout)      :: self
    integer(shortInt), intent(in)               :: idx, elementIdx
    integer(shortInt), dimension(4), intent(in) :: vertexIdxs
    real(defReal), dimension(3), intent(in)     :: centroid

    call self % shelf(idx) % setIdx(idx)
    call self % shelf(idx) % setElement(elementIdx)
    call self % shelf(idx) % setVertices(vertexIdxs)
    call self % shelf(idx) % setCentroid(centroid)

  end subroutine initTetrahedron
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(tetrahedronShelf), intent(inout) :: self
    
    if (allocated(self % shelf)) deallocate(self % shelf)

  end subroutine kill

  !! Subroutine 'testForInclusion'
  !!
  !! Basic description:
  !!   Tests whether a given tetrahedron in the shelf contains a point.
  !!
  !! Arguments:
  !!   idx [in]                  -> Index of the tetrahedron in the shelf.
  !!   r [in]                    -> 3-D coordinates of the point.
  !!   triangles [in]            -> A triangleShelf.
  !!   failedTriangleIdx [out]   -> Index of the first tetrahedron's triangle for which the test fails.
  !!   surfTolTriangleIdxs [out] -> Indices of the tetrahedron's triangles on which the point lies.
  !!
  pure subroutine testForInclusion(self, idx, r, triangles, failedTriangleIdx, surfTolTriangleIdxs)
    class(tetrahedronShelf), intent(in)                       :: self
    integer(shortInt), intent(in)                             :: idx
    real(defReal), dimension(3), intent(in)                   :: r
    type(triangleShelf), intent(in)                           :: triangles
    integer(shortInt), intent(out)                            :: failedTriangleIdx
    integer(shortInt), dimension(:), allocatable, intent(out) :: surfTolTriangleIdxs

    call self % shelf(idx) % testForInclusion(triangles, r, failedTriangleIdx, surfTolTriangleIdxs)

  end subroutine testForInclusion
  
end module tetrahedronShelf_class