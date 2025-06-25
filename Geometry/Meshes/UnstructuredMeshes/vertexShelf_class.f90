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
    procedure                               :: allocateShelf
    procedure                               :: expandShelf
    procedure                               :: findCommonEdgeIdx
    procedure                               :: findCommonFaceIdx
    procedure                               :: getAllCoordinates
    procedure                               :: getExtremalCoordinates
    procedure                               :: getOffset
    procedure                               :: getSize
    generic                                 :: getVertexCoordinates => getVertexCoordinates_shortInt, &
                                                                       getVertexCoordinates_shortIntArray
    procedure, private                      :: getVertexCoordinates_shortInt
    procedure, private                      :: getVertexCoordinates_shortIntArray
    procedure                               :: getVertexEdgeIdxs
    procedure                               :: getVertexElementIdxs
    generic                                 :: getVertexFaceIdxs => getVertexFaceIdxs_shortInt, &
                                                                    getVertexFaceIdxs_shortIntArray
    procedure, private                      :: getVertexFaceIdxs_shortInt
    procedure, private                      :: getVertexFaceIdxs_shortIntArray
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
  pure function findCommonFaceIdx(self, vertexIdxs) result(faceIdx)
    class(vertexShelf), intent(in)               :: self
    integer(shortInt), dimension(3), intent(in)  :: vertexIdxs
    integer(shortInt)                            :: faceIdx, i
    integer(shortInt), dimension(:), allocatable :: commonIdxs

    ! Initialise triangleIdx = 0 and return immediately if any vertices are not associated with triangles.
    faceIdx = 0
    if (any(.not. self % shelf(vertexIdxs) % hasFaces(), 1)) return
    
    ! Find common edge indices. Update edgeIdx only if common indices have been found.
    commonIdxs = self % shelf(vertexIdxs(1)) % getFaceIdxs()
    do i = 2, 3
      commonIdxs = findCommon(commonIdxs, self % shelf(vertexIdxs(i)) % getFaceIdxs())

    end do
    if (size(commonIdxs) > 0) faceIdx = commonIdxs(1)

  end function findCommonFaceIdx

  !! Function 'getAllCoordinates'
  !!
  !! Basic description:
  !!   Returns the 3-D coordinates of all the vertices in the shelf.
  !!
  !! Result:
  !!   allCoordinates -> Array listing the 3-D coordinates of all the vertices.
  !!
  pure function getAllCoordinates(self) result(allCoordinates)
    class(vertexShelf), intent(in)                  :: self
    real(defReal), dimension(3, size(self % shelf)) :: allCoordinates
    integer(shortInt)                               :: i

    do i = 1, size(self % shelf)
      allCoordinates(:, i) = self % shelf(i) % getCoordinates()

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

  !! Function 'getVertexCoordinates_shortInt'
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
  pure function getVertexCoordinates_shortInt(self, idx) result(coords)
    class(vertexShelf), intent(in) :: self
    integer(shortInt), intent(in)  :: idx
    real(defReal), dimension(3)    :: coords

    coords = self % shelf(idx) % getCoordinates()

  end function getVertexCoordinates_shortInt

  !! Function 'getVertexCoordinates_shortInt'
  !!
  !! Basic description:
  !!   Returns the 3-D coordinates of vertices in the shelf.
  !!
  !! Arguments:
  !!   idxs [in] -> Indices of the vertices in the shelf.
  !!
  !! Result:
  !!   coords   -> 3-D coordinates of the vertices.
  !!
  pure function getVertexCoordinates_shortIntArray(self, idxs) result(coords)
    class(vertexShelf), intent(in)              :: self
    integer(shortInt), dimension(:), intent(in) :: idxs
    real(defReal), dimension(3, size(idxs))     :: coords
    integer(shortInt)                           :: i

    do i = 1, size(idxs)
      coords(:, i) = self % shelf(idxs(i)) % getCoordinates()

    end do

  end function getVertexCoordinates_shortIntArray

  !! Function 'getVertexEdgeIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of all the edges containing a vertex in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the vertex in the shelf.
  !!
  !! Result:
  !!   edgeIdxs -> Indices of all the edges containing the vertex.
  !!
  pure function getVertexEdgeIdxs(self, idx) result(edgeIdxs)
    class(vertexShelf), intent(in)               :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: edgeIdxs

    edgeIdxs = self % shelf(idx) % getEdgeIdxs()

  end function getVertexEdgeIdxs

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

    elementIdxs = self % shelf(idx) % getElementIdxs()

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

    faceIdxs = self % shelf(idx) % getFaceIdxs()

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
      call append(faceIdxs, self % shelf(idxs(i)) % getFaceIdxs(), .true.)

    end do

  end function getVertexFaceIdxs_shortIntArray

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