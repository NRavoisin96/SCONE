module edgeShelf_class
  
  use edge_class,        only : edge
  use genericProcedures, only : fatalError, numToChar
  use numPrecision
  
  implicit none
  private
  
  type, public                            :: edgeShelf
    private
    type(edge), dimension(:), allocatable :: shelf
  contains
    procedure                             :: addElementIdxToEdge
    procedure                             :: addFaceIdxToEdge
    procedure                             :: addTetrahedronIdxToEdge
    procedure                             :: addTriangleIdxToEdge
    procedure                             :: allocateShelf
    procedure                             :: collapseShelf
    procedure                             :: expandShelf
    procedure                             :: getEdgeElementIdxs
    procedure                             :: getEdgeFaceIdxs
    procedure                             :: getEdgeTetrahedronIdxs
    procedure                             :: getEdgeTriangleIdxs
    procedure                             :: initEdge
    procedure                             :: kill
  end type edgeShelf

contains

  !! Subroutine 'addElementIdxToEdge'
  !!
  !! Basic description:
  !!   Adds the index of an element to an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]        -> Index of the edge in the shelf.
  !!   elementIdx [in] -> Index of the element containing the edge.
  !!
  elemental subroutine addElementIdxToEdge(self, idx, elementIdx)
    class(edgeShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx, elementIdx

    call self % shelf(idx) % addElementIdx(elementIdx)

  end subroutine addElementIdxToEdge

  !! Subroutine 'addFaceIdxToEdge'
  !!
  !! Basic description:
  !!   Adds the index of a face to an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the edge in the shelf.
  !!   faceIdx [in] -> Index of the face containing the edge.
  !!
  elemental subroutine addFaceIdxToEdge(self, idx, faceIdx)
    class(edgeShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx, faceIdx

    call self % shelf(idx) % addFaceIdx(faceIdx)

  end subroutine addFaceIdxToEdge

  !! Subroutine 'addTetrahedronIdxToEdge'
  !!
  !! Basic description:
  !!   Adds the index of a tetrahedron to an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]            -> Index of the edge in the shelf.
  !!   tetrahedronIdx [in] -> Index of the tetrahedron containing the edge.
  !!
  elemental subroutine addTetrahedronIdxToEdge(self, idx, tetrahedronIdx)
    class(edgeShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx, tetrahedronIdx

    call self % shelf(idx) % addTetrahedronIdx(tetrahedronIdx)

  end subroutine addTetrahedronIdxToEdge

  !! Subroutine 'addTriangleIdxToEdge'
  !!
  !! Basic description:
  !!   Adds the index of a triangle to an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]         -> Index of the edge in the shelf.
  !!   triangleIdx [in] -> Index of the triangle containing the edge.
  !!
  elemental subroutine addTriangleIdxToEdge(self, idx, triangleIdx)
    class(edgeShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx, triangleIdx

    call self % shelf(idx) % addTriangleIdx(triangleIdx)

  end subroutine addTriangleIdxToEdge

  !! Subroutine 'allocateShelf'
  !!
  !! Basic description:
  !!   Allocates memory in the shelf.
  !!
  !! Arguments:
  !!   nEdges [in] -> Number of edges to be included in the shelf.
  !!
  elemental subroutine allocateShelf(self, nEdges)
    class(edgeShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: nEdges

    allocate(self % shelf(nEdges))

  end subroutine allocateShelf

  !! Subroutine 'collapseShelf'
  !!
  !! Basic description:
  !!   Reduces the size of the shelf to lastIdx.
  !!
  !! Arguments:
  !!   lastIdx [in] -> Index of the last edge in the shelf to be collapsed.
  !!
  elemental subroutine collapseShelf(self, lastIdx)
    class(edgeShelf), intent(inout)       :: self
    integer(shortInt), intent(in)         :: lastIdx
    type(edge), dimension(:), allocatable :: shelf
    
    if (allocated(self % shelf)) then
      ! Create a temporary shelf and copy all the elements up to lastIdx from the 
      ! original shelf.
      shelf = self % shelf(1:lastIdx)
      
      ! Deallocate and reallocate shelf then copy elements back.
      deallocate(self % shelf)
      allocate(self % shelf(lastIdx))
      self % shelf = shelf

    else
      allocate(self % shelf(lastIdx))

    end if

  end subroutine collapseShelf

  !! Subroutine 'expandShelf'
  !!
  !! Basic description:
  !!   Expands the shelf by a specified number of additional edges. Copies elements
  !!   already present. Allocates the shelf if it is not allocated yet.
  !!
  !! Arguments:
  !!   nAdditionalEdges [in] -> Number of additional edges to be included in the shelf.
  !!
  elemental subroutine expandShelf(self, nAdditionalEdges)
    class(edgeShelf), intent(inout)       :: self
    integer(shortInt), intent(in)         :: nAdditionalEdges
    integer(shortInt)                     :: nEdges
    type(edge), dimension(:), allocatable :: shelf

    if (allocated(self % shelf)) then
      ! If shelf is already allocated, compute the number of edges in the shelf to be expanded
      ! and copy elements already present.
      nEdges = size(self % shelf)
      shelf = self % shelf
      
      ! Deallocate shelf and reallocate to new size then copy original elements.
      deallocate(self % shelf)
      allocate(self % shelf(nEdges + nAdditionalEdges))
      self % shelf(1:nEdges) = shelf

    else
      allocate(self % shelf(nAdditionalEdges))

    end if

  end subroutine expandShelf

  !! Function 'getEdgeElementIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the elements sharing an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]    -> Index of the edge in the shelf.
  !!
  !! Result:
  !!   elementIdxs -> Indices of the elements sharing the edge.
  !!
  pure function getEdgeElementIdxs(self, idx) result(elementIdxs)
    class(edgeShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: elementIdxs

    elementIdxs = self % shelf(idx) % getElementIdxs()

  end function getEdgeElementIdxs

  !! Function 'getEdgeFaceIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the faces sharing an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the edge in the shelf.
  !!
  !! Result:
  !!   faceIdxs -> Indices of the faces sharing the edge.
  !!
  pure function getEdgeFaceIdxs(self, idx) result(faceIdxs)
    class(edgeShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: faceIdxs

    faceIdxs = self % shelf(idx) % getFaceIdxs()

  end function getEdgeFaceIdxs

  !! Function 'getEdgeTetrahedronIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the tetrahedra sharing an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]        -> Index of the edge in the shelf.
  !!
  !! Result:
  !!   tetrahedronIdxs -> Indices of the tetrahedra sharing the edge.
  !!
  pure function getEdgeTetrahedronIdxs(self, idx) result(tetrahedronIdxs)
    class(edgeShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: tetrahedronIdxs

    tetrahedronIdxs = self % shelf(idx) % getTetrahedronIdxs()

  end function getEdgeTetrahedronIdxs

  !! Function 'getEdgeTriangleIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the triangles sharing an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the edge in the shelf.
  !!
  !! Result:
  !!   triangleIdxs -> Indices of the triangles sharing the edge.
  !!
  pure function getEdgeTriangleIdxs(self, idx) result(triangleIdxs)
    class(edgeShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: triangleIdxs

    triangleIdxs = self % shelf(idx) % getTriangleIdxs()

  end function getEdgeTriangleIdxs

  !! Subroutine 'initEdge'
  !!
  !! Basic description:
  !!   Initialises an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]        -> Index of the edge in the shelf.
  !!   vertexIdxs [in] -> Indices of the vertices in the edge.
  !!
  pure subroutine initEdge(self, idx, vertexIdxs)
    class(edgeShelf), intent(inout)             :: self
    integer(shortInt), intent(in)               :: idx
    integer(shortInt), dimension(2), intent(in) :: vertexIdxs

    call self % shelf(idx) % setIdx(idx)
    call self % shelf(idx) % setVertexIdxs(vertexIdxs)

  end subroutine initEdge

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an unitialised state.
  !!
  elemental subroutine kill(self)
    class(edgeShelf), intent(inout) :: self

    if (allocated(self % shelf)) deallocate(self % shelf)

  end subroutine kill
  
end module edgeShelf_class