module vertexShelf_class
  
  use numPrecision
  use genericProcedures, only : append, findCommon
  use vertex_class,      only : vertex
  
  implicit none
  private
  
  !!
  !! Storage space for vertices of a given OpenFOAM mesh.
  !!
  !! Public members:
  !!   shelf               -> Array to store vertices.
  !!   offset              -> User-supplied 3-D offset applied to the coordinates 
  !!                          of all the vertices in the shelf.
  !!   extremalCoordinates -> Array of minimum and maximum x-, y- and z-
  !!                          coordinates in the shelf.
  !!
  type, public :: vertexShelf
    private
    type(vertex), dimension(:), allocatable :: shelf
    real(defReal), dimension(3)             :: offset = ZERO
    real(defReal), dimension(6)             :: extremalCoordinates = ZERO
  contains
    procedure                               :: addEdgeIdxToVertex
    procedure                               :: addElementIdxToVertex
    procedure                               :: addFaceIdxToVertex
    procedure                               :: addTetrahedronIdxToVertex
    procedure                               :: addTriangleIdxToVertex
    procedure                               :: allocateShelf
    procedure                               :: expandShelf
    procedure                               :: findCommonEdgeIdx
    procedure                               :: findCommonTriangleIdx
    procedure                               :: getAllCoordinates
    procedure                               :: getExtremalCoordinates
    procedure                               :: getOffset
    procedure                               :: getSize
    procedure                               :: getVertexCoordinates
    procedure                               :: getVertexElementIdxs
    generic                                 :: getVertexFaceIdxs => getVertexFaceIdxs_shortInt, &
                                                                    getVertexFaceIdxs_shortIntArray
    procedure, private                      :: getVertexFaceIdxs_shortInt
    procedure, private                      :: getVertexFaceIdxs_shortIntArray
    procedure                               :: getVertexTetrahedronIdxs
    generic                                 :: getVertexTriangleIdxs => getVertexTriangleIdxs_shortInt, &
                                                                        getVertexTriangleIdxs_shortIntArray
    procedure, private                      :: getVertexTriangleIdxs_shortInt
    procedure, private                      :: getVertexTriangleIdxs_shortIntArray
    procedure                               :: initVertex
    procedure                               :: kill
    procedure                               :: setExtremalCoordinates
    procedure                               :: setOffset
  end type 

contains

  !! Subroutine 'addEdgeIdxToVertex'
  !!
  !! Basic description:
  !!   Adds the index of an edge to a vertex in the shelf.
  !!
  !! Arguments:
  !!   vertexIdx [in] -> Index of the vertex in the shelf.
  !!   edgeIdx [in]   -> Index of the edge containing the vertex.
  !!
  elemental subroutine addEdgeIdxToVertex(self, vertexIdx, edgeIdx)
    class(vertexShelf), intent(inout) :: self
    integer(shortInt), intent(in)     :: vertexIdx, edgeIdx

    call self % shelf(vertexIdx) % addEdgeIdx(edgeIdx)

  end subroutine addEdgeIdxToVertex

  !! Subroutine 'addElementIdxToVertex'
  !!
  !! Basic description:
  !!   Adds the index of an element to a vertex in the shelf.
  !!
  !! Arguments:
  !!   vertexIdx [in]  -> Index of the vertex in the shelf.
  !!   elementIdx [in] -> Index of the element containing the vertex.
  !!
  elemental subroutine addElementIdxToVertex(self, vertexIdx, elementIdx)
    class(vertexShelf), intent(inout) :: self
    integer(shortInt), intent(in)     :: vertexIdx, elementIdx

    call self % shelf(vertexIdx) % addElementIdx(elementIdx)

  end subroutine addElementIdxToVertex

  !! Subroutine 'addFaceIdxToVertex'
  !!
  !! Basic description:
  !!   Adds the index of a face to a vertex in the shelf.
  !!
  !! Arguments:
  !!   vertexIdx [in] -> Index of the vertex in the shelf.
  !!   faceIdx [in]   -> Index of the face containing the vertex.
  !!
  elemental subroutine addFaceIdxToVertex(self, vertexIdx, faceIdx)
    class(vertexShelf), intent(inout) :: self
    integer(shortInt), intent(in)     :: vertexIdx, faceIdx

    call self % shelf(vertexIdx) % addFaceIdx(faceIdx)

  end subroutine addFaceIdxToVertex

  !! Subroutine 'addTetrahedronIdxToVertex'
  !!
  !! Basic description:
  !!   Adds the index of a tetrahedron to a vertex in the shelf.
  !!
  !! Arguments:
  !!   vertexIdx [in]      -> Index of the vertex in the shelf.
  !!   tetrahedronIdx [in] -> Index of the tetrahedron containing the vertex.
  !!
  elemental subroutine addTetrahedronIdxToVertex(self, vertexIdx, tetrahedronIdx)
    class(vertexShelf), intent(inout) :: self
    integer(shortInt), intent(in)     :: vertexIdx, tetrahedronIdx

    call self % shelf(vertexIdx) % addTetrahedronIdx(tetrahedronIdx)

  end subroutine addTetrahedronIdxToVertex

  !! Subroutine 'addTriangleIdxToVertex'
  !!
  !! Basic description:
  !!   Adds the index of a triangle to a vertex in the shelf.
  !!
  !! Arguments:
  !!   vertexIdx [in]   -> Index of the vertex in the shelf.
  !!   triangleIdx [in] -> Index of the triangle containing the vertex.
  !!
  elemental subroutine addTriangleIdxToVertex(self, vertexIdx, triangleIdx)
    class(vertexShelf), intent(inout) :: self
    integer(shortInt), intent(in)     :: vertexIdx, triangleIdx

    call self % shelf(vertexIdx) % addTriangleIdx(triangleIdx)

  end subroutine addTriangleIdxToVertex

  !! Subroutine 'allocateShelf'
  !!
  !! Basic description:
  !!   Allocates memory in the shelf.
  !!
  !! Arguments:
  !!   nVertices [in] -> Number of vertices to be included in the shelf.
  !!
  elemental subroutine allocateShelf(self, nVertices)
    class(vertexShelf), intent(inout) :: self
    integer(shortInt), intent(in)     :: nVertices

    allocate(self % shelf(nVertices))

  end subroutine allocateShelf

  !! Subroutine 'expandShelf'
  !!
  !! Basic description:
  !!   Expands the shelf by a specified number of additional vertices. Copies elements
  !!   already present. Allocates the shelf if it is not allocated yet.
  !!
  !! Arguments:
  !!   nAdditionalVertices [in] -> Number of additional vertices to be included in the shelf.
  !!
  elemental subroutine expandShelf(self, nAdditionalVertices)
    class(vertexShelf), intent(inout)       :: self
    integer(shortInt), intent(in)           :: nAdditionalVertices
    integer(shortInt)                       :: nVertices
    type(vertex), dimension(:), allocatable :: shelf

    if (allocated(self % shelf)) then
      ! If shelf is already allocated, compute the number of vertices in the shelf to be expanded
      ! and copy elements already present.
      nVertices = size(self % shelf)
      shelf = self % shelf
      
      ! Deallocate shelf and reallocate to new size then copy original elements.
      deallocate(self % shelf)
      allocate(self % shelf(nVertices + nAdditionalVertices))
      self % shelf(1:nVertices) = shelf

    else
      allocate(self % shelf(nAdditionalVertices))

    end if

  end subroutine expandShelf

  !! Function 'findCommonEdgeIdx'
  !!
  !! Basic description:
  !!   Finds the index of the edge containing two vertices.
  !!
  !! Arguments:
  !!   firstVertexIdx [in]  -> Index of the first vertex in the edge.
  !!   secondVertexIdx [in] -> Index of the second vertex in the edge.
  !!
  !! Result:
  !!   edgeIdx              -> Index of the edge containing the two vertices.
  !!
  elemental function findCommonEdgeIdx(self, firstVertexIdx, secondVertexIdx) result(edgeIdx)
    class(vertexShelf), intent(in)               :: self
    integer(shortInt), intent(in)                :: firstVertexIdx, secondVertexIdx
    integer(shortInt)                            :: edgeIdx
    integer(shortInt), dimension(:), allocatable :: commonIdxs

    ! Initialise edgeIdx = 0 and return immediately if any vertices are not associated with edges.
    edgeIdx = 0
    if (.not. self % shelf(firstVertexIdx) % hasEdges() .or. .not. self % shelf(secondVertexIdx) % hasEdges()) return
    
    ! Find common edge indices. Update edgeIdx only if common indices have been found.
    commonIdxs = findCommon(self % shelf(firstVertexIdx) % getEdgeIdxs(), self % shelf(secondVertexIdx) % getEdgeIdxs())
    if (size(commonIdxs) > 0) edgeIdx = commonIdxs(1)

  end function findCommonEdgeIdx

  !! Function 'findCommonTriangleIdx'
  !!
  !! Basic description:
  !!   Finds the index of the triangle containing three vertices.
  !!
  !! Arguments:
  !!   vertexIdxs [in] -> Indices of the vertices.
  !!
  !! Result:
  !!   triangleIdx     -> Index of the triangle containing the two vertices.
  !!
  pure function findCommonTriangleIdx(self, vertexIdxs) result(triangleIdx)
    class(vertexShelf), intent(in)               :: self
    integer(shortInt), dimension(3), intent(in)  :: vertexIdxs
    integer(shortInt)                            :: triangleIdx, i
    integer(shortInt), dimension(:), allocatable :: commonIdxs

    ! Initialise triangleIdx = 0 and return immediately if any vertices are not associated with triangles.
    triangleIdx = 0
    if (any(.not. self % shelf(vertexIdxs) % hasTriangles(), 1)) return
    
    ! Find common edge indices. Update edgeIdx only if common indices have been found.
    commonIdxs = self % shelf(vertexIdxs(1)) % getVertexToTriangles()
    do i = 2, 3
      commonIdxs = findCommon(commonIdxs, self % shelf(vertexIdxs(i)) % getVertexToTriangles())

    end do
    if (size(commonIdxs) > 0) triangleIdx = commonIdxs(1)

  end function findCommonTriangleIdx

  !! Function 'getAllCoordinates'
  !!
  !! Basic description:
  !!   Returns the 3-D coordinates of all the vertices in the shelf.
  !!
  !! Result:
  !!   allCoordinates -> Array listing the 3-D coordinates of all the vertices.
  !!
  pure function getAllCoordinates(self) result(allCoordinates)
    class(vertexShelf), intent(in)                :: self
    real(defReal), dimension(self % getSize(), 3) :: allCoordinates
    integer(shortInt)                             :: i

    do i = 1, self % getSize()
      allCoordinates(i, :) = self % shelf(i) % getCoordinates()

    end do

  end function getAllCoordinates
  
  !! Function 'getExtremalCoordinates'
  !!
  !! Basic description:
  !!   Returns the minimum and maximum x-, y- and z- coordinates of the vertices in the shelf.
  !!
  !! Result:
  !!   extremalCoordinates -> Array containing six entries: the first three list the minimum x-, y-
  !!   and z-coordinates, while the last three list the maximum x-, y- and z-coordinates.
  !!
  pure function getExtremalCoordinates(self) result(extremalCoordinates)
    class(vertexShelf), intent(in)                :: self
    real(defReal), dimension(6)                   :: extremalCoordinates
    
    extremalCoordinates = self % extremalCoordinates

  end function getExtremalCoordinates

  !! Function 'getOffset'
  !!
  !! Basic description:
  !!   Returns the offset of the vertices in the shelf.
  !!
  !! Result:
  !!   offset -> Offset of all the vertices in the shelf.
  !!
  pure function getOffset(self) result(offset)
    class(vertexShelf), intent(in) :: self
    real(defReal), dimension(3)    :: offset

    offset = self % offset

  end function getOffset
  
  !! Function 'getSize'
  !!
  !! Basic description:
  !!   Returns the size of the shelf.
  !!
  !! Result:
  !!   nVertices -> Size of the shelf.
  !!
  elemental function getSize(self) result(nVertices)
    class(vertexShelf), intent(in) :: self
    integer(shortInt)              :: nVertices
    
    nVertices = size(self % shelf)

  end function getSize

  !! Function 'getVertexCoordinates'
  !!
  !! Basic description:
  !!   Returns the 3-D coordinates of a vertex in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the vertex in the shelf.
  !!
  !! Result:
  !!   coords   -> 3-D coordinates of the vertex.
  !!
  pure function getVertexCoordinates(self, idx) result(coords)
    class(vertexShelf), intent(in) :: self
    integer(shortInt), intent(in)  :: idx
    real(defReal), dimension(3)    :: coords

    coords = self % shelf(idx) % getCoordinates()

  end function getVertexCoordinates

  !! Function 'getVertexElementIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of all the elements containing a vertex in the shelf.
  !!
  !! Arguments:
  !!   idx [in]    -> Index of the vertex in the shelf.
  !!
  !! Result:
  !!   elementIdxs -> Indices of all the elements containing the vertex.
  !!
  pure function getVertexElementIdxs(self, idx) result(elementIdxs)
    class(vertexShelf), intent(in)               :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: elementIdxs

    elementIdxs = self % shelf(idx) % getVertexToElements()

  end function getVertexElementIdxs

  !! Function 'getVertexFaceIdxs_shortInt'
  !!
  !! Basic description:
  !!   Returns the indices of all the faces containing a vertex in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the vertex in the shelf.
  !!
  !! Result:
  !!   faceIdxs -> Indices of all the faces containing the vertex.
  !!
  pure function getVertexFaceIdxs_shortInt(self, idx) result(faceIdxs)
    class(vertexShelf), intent(in)               :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: faceIdxs

    faceIdxs = self % shelf(idx) % getVertexToFaces()

  end function getVertexFaceIdxs_shortInt

  !! Function 'getVertexFaceIdxs_shortIntArray'
  !!
  !! Basic description:
  !!   Returns the unique indices of all the faces containing a set of vertices
  !!   in the shelf.
  !!
  !! Arguments:
  !!   idxs [in] -> Indices of the vertices in the shelf.
  !!
  !! Result:
  !!   faceIdxs  -> Unique indices of all the faces containing the vertices.
  !!
  pure function getVertexFaceIdxs_shortIntArray(self, idxs) result(faceIdxs)
    class(vertexShelf), intent(in)               :: self
    integer(shortInt), dimension(:), intent(in)  :: idxs
    integer(shortInt), dimension(:), allocatable :: faceIdxs
    integer(shortInt)                            :: i

    do i = 1, size(idxs)
      call append(faceIdxs, self % shelf(idxs(i)) % getVertexToFaces(), .true.)

    end do

  end function getVertexFaceIdxs_shortIntArray

  !! Function 'getVertexTetrahedronIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of all the tetrahedra containing a vertex in the shelf.
  !!
  !! Arguments:
  !!   idx [in]        -> Index of the vertex in the shelf.
  !!
  !! Result:
  !!   tetrahedronIdxs -> Indices of all the tetrahedra containing the vertex.
  !!
  pure function getVertexTetrahedronIdxs(self, idx) result(tetrahedronIdxs)
    class(vertexShelf), intent(in)               :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: tetrahedronIdxs

    tetrahedronIdxs = self % shelf(idx) % getVertexToTetrahedra()

  end function getVertexTetrahedronIdxs

  !! Function 'getVertexTriangleIdxs_shortInt'
  !!
  !! Basic description:
  !!   Returns the indices of all the triangles containing a vertex in the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the vertex in the shelf.
  !!
  !! Result:
  !!   triangleIdxs -> Indices of all the triangles containing the vertex.
  !!
  pure function getVertexTriangleIdxs_shortInt(self, idx) result(triangleIdxs)
    class(vertexShelf), intent(in)               :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: triangleIdxs

    triangleIdxs = self % shelf(idx) % getVertexToTriangles()

  end function getVertexTriangleIdxs_shortInt

  !! Function 'getVertexTriangleIdxs_shortIntArray'
  !!
  !! Basic description:
  !!   Returns the unique indices of all the triangles containing a set of vertices
  !!   in the shelf.
  !!
  !! Arguments:
  !!   idxs [in]    -> Indices of the vertices in the shelf.
  !!
  !! Result:
  !!   triangleIdxs -> Unique indices of all the triangles containing the vertices.
  !!
  pure function getVertexTriangleIdxs_shortIntArray(self, idxs) result(triangleIdxs)
    class(vertexShelf), intent(in)               :: self
    integer(shortInt), dimension(:), intent(in)  :: idxs
    integer(shortInt), dimension(:), allocatable :: triangleIdxs
    integer(shortInt)                            :: i

    do i = 1, size(idxs)
      call append(triangleIdxs, self % shelf(idxs(i)) % getVertexToTriangles(), .true.)

    end do

  end function getVertexTriangleIdxs_shortIntArray

  !! Subroutine 'initVertex'
  !!
  !! Basic description:
  !!   Initialises a vertex in the shelf.
  !!
  !! Arguments:
  !!   idx [in]    -> Index of the vertex.
  !!   coords [in] -> 3-D coordinates of the vertex.
  !!
  pure subroutine initVertex(self, idx, coords)
    class(vertexShelf), intent(inout)       :: self
    integer(shortInt), intent(in)           :: idx
    real(defReal), dimension(3), intent(in) :: coords

    call self % shelf(idx) % setIdx(idx)
    call self % shelf(idx) % setCoordinates(coords)

  end subroutine initVertex
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(vertexShelf), intent(inout) :: self

    self % offset = ZERO
    self % extremalCoordinates = ZERO
    if (allocated(self % shelf)) deallocate(self % shelf)

  end subroutine kill

  !! Subroutine 'setExtremalCoordinates'
  !!
  !! Basic description:
  !!   Sets the extremal coordinates of the vertices in the mesh.
  !!
  !! Arguments:
  !!   coords [in] -> defReal array of extremal coordinates (x_min, y_min, z_min, x_max, y_max and z_max).
  !!
  pure subroutine setExtremalCoordinates(self, coords)
    class(vertexShelf), intent(inout)       :: self
    real(defReal), dimension(6), intent(in) :: coords

    self % extremalCoordinates = coords

  end subroutine setExtremalCoordinates

  !! Subroutine 'setOffset'
  !!
  !! Basic description:
  !!   Sets the offset of the shelf. This is a 3-D translation vector applied to
  !!   the coordinates of all the vertices in the shelf.
  !!
  !! Arguments:
  !!   offset [in] -> 3-D offset.
  !!
  pure subroutine setOffset(self, offset)
    class(vertexShelf), intent(inout)       :: self
    real(defReal), dimension(3), intent(in) :: offset
    
    self % offset = offset

  end subroutine setOffset

end module vertexShelf_class