module edgeShelf_class
  
  use edge_class,        only : edge
  use genericProcedures, only : numToChar
  use numPrecision
  use vertexShelf_class, only : vertexShelf
  
  implicit none
  private
  
  type, public                            :: edgeShelf
    private
    type(edge), dimension(:), allocatable :: shelf
  contains
    procedure                             :: addElementIdxToEdge
    procedure                             :: addFaceIdxToEdge
    procedure                             :: allocateShelf
    procedure                             :: collapseShelf
    procedure                             :: expandShelf
    procedure                             :: findElementIdxFromEdgeAngularSectorSearch
    procedure                             :: getEdgeElementIdxs
    procedure                             :: getEdgeFaceIdxs
    procedure                             :: getEdgeVertexIdxs
    procedure                             :: getSize
    procedure                             :: getEdgeUnitVector
    procedure                             :: getEdgeVector
    procedure                             :: getEdgeLocalBasis1
    procedure                             :: getEdgeLocalBasis2
    procedure                             :: getEdgeLength
    procedure                             :: getEdgeDotProductOfVector
    procedure                             :: getEdgeAnglesArray
    procedure                             :: getEdgeElementIdxsArray
    procedure                             :: getEdgeIsBoundary
    procedure                             :: setEdgeUnitVector
    procedure                             :: setEdgeVector ! This can be removed and retrieved from UnitVector (memory vs time)
    procedure                             :: setEdgeLocalBasis1
    procedure                             :: setEdgeLocalBasis2
    procedure                             :: setEdgeLength
    procedure                             :: setEdgeDotProductOfVector
    procedure                             :: setEdgeAnglesArray
    procedure                             :: setEdgeElementIdxsArray
    procedure                             :: setEdgeIsboundary
    procedure                             :: isAllocatedEdgeAnglesArray
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

  !!
  !!
  !!
  pure subroutine findElementIdxFromEdgeAngularSectorSearch(self, idx, r, vertices, elementIdx)
    class(edgeShelf), intent(in)            :: self
    integer(shortInt), intent(in)           :: idx
    real(defReal), dimension(3), intent(in) :: r
    type(vertexShelf), intent(in)           :: vertices
    integer(shortInt), intent(inout)        :: elementIdx

    call self % shelf(idx) % findElementIdxFromAngularSectorSearch(r, vertices, elementIdx)

  end subroutine findElementIdxFromEdgeAngularSectorSearch

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

  !! Function 'getEdgeVertexIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in an edge of the shelf.
  !!
  !! Arguments:
  !!   idx [in]   -> Index of the edge in the shelf.
  !!
  !! Result:
  !!   vertexIdxs -> Indices of the vertices in the edge.
  !!
  pure function getEdgeVertexIdxs(self, idx) result(vertexIdxs)
    class(edgeShelf), intent(in)    :: self
    integer(shortInt), intent(in)   :: idx
    integer(shortInt), dimension(2) :: vertexIdxs

    vertexIdxs = self % shelf(idx) % getVertexIdxs()

  end function getEdgeVertexIdxs

  !! Function 'getSize'
  !!
  !! Basic description:
  !!   Returns the number of edges in the shelf.
  !!
  !! Result:
  !!   nEdges -> Number of edges in the shelf.
  !!
  elemental function getSize(self) result(nEdges)
    class(edgeShelf), intent(in) :: self
    integer(shortInt)            :: nEdges

    nEdges = size(self % shelf)

  end function getSize

  !!
  !!
  !!
  pure function getEdgeUnitVector(self, idx) result(unitVector)
    class(edgeShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    real(defReal), dimension(3)   :: unitVector

    unitVector = self % shelf(idx) % getUnitVector()

  end function getEdgeUnitVector

  !!
  !!
  !!
  pure function getEdgeVector(self, idx) result(vector)
    class(edgeShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    real(defReal), dimension(3)   :: vector

    vector = self % shelf(idx) % getVector()

  end function getEdgevector

  !!
  !!
  !!
  pure function getEdgeLocalBasis1(self, idx) result(localBasis1)
    class(edgeShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    real(defReal), dimension(3)   :: localBasis1

    localBasis1 = self % shelf(idx) % getLocalBasis1()

  end function getEdgeLocalBasis1

  !!
  !!
  !!
  pure function getEdgeLocalBasis2(self, idx) result(localBasis2)
    class(edgeShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    real(defReal), dimension(3)   :: localBasis2

    localBasis2 = self % shelf(idx) % getLocalBasis2()

  end function getEdgeLocalBasis2

  !!
  !!
  !!
  elemental function getEdgeLength(self, idx) result(length)
    class(edgeShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    real(defReal)                 :: length

    length = self % shelf(idx) % getLength()

  end function getEdgeLength

  !!
  !!
  !!
  elemental function getEdgeDotProductOfVector(self, idx) result(dotProductOfVector)
    class(edgeShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    real(defReal)                 :: dotProductOfVector

    dotProductOfVector = self % shelf(idx) % getDotProductOfVector()

  end function getEdgeDotProductOfVector

  !!
  !!
  !!
  pure function getEdgeAnglesArray(self, idx) result(anglesArray)
    class(edgeShelf), intent(in)                :: self
    integer(shortInt), intent(in)               :: idx
    real(defReal), dimension(:, :), allocatable :: anglesArray

    anglesArray = self % shelf(idx) % getAnglesArray()

  end function getEdgeAnglesArray

  !!
  !!
  !!
  pure function getEdgeElementIdxsArray(self, idx) result(elementIdxsArray)
    class(edgeShelf), intent(in)                   :: self
    integer(shortInt), intent(in)                  :: idx
    integer(shortInt), dimension(:), allocatable   :: elementIdxsArray

    elementIdxsArray = self % shelf(idx) % getElementIdxsArray()

  end function getEdgeElementIdxsArray

  !!
  !!
  !!
  elemental function getEdgeIsBoundary(self, idx) result(isBoundary)
    class(edgeShelf), intent(in)                   :: self
    integer(shortInt), intent(in)                  :: idx
    logical                                        :: isBoundary

    isBoundary = self % shelf(idx) % getIsBoundary()

  end function getEdgeIsBoundary

  !!
  !!
  !!
  pure subroutine setEdgeUnitVector(self, idx, unitVector)
    class(edgeShelf), intent(inout)             :: self
    integer(shortInt), intent(in)               :: idx
    real(defReal), dimension(3), intent(in)     :: unitVector

    call self % shelf(idx) % setUnitVector(unitVector)

  end subroutine setEdgeUnitVector

  !!
  !!
  !!
  pure subroutine setEdgeVector(self, idx, vector)
    class(edgeShelf), intent(inout)             :: self
    integer(shortInt), intent(in)               :: idx
    real(defReal), dimension(3), intent(in)     :: vector

    call self % shelf(idx) % setVector(vector)

  end subroutine setEdgeVector

  !!
  !!
  !!
  pure subroutine setEdgeLocalBasis1(self, idx, localBasis1)
    class(edgeShelf), intent(inout)             :: self
    integer(shortInt), intent(in)               :: idx
    real(defReal), dimension(3), intent(in)     :: localBasis1

    call self % shelf(idx) % setLocalBasis1(localBasis1)

  end subroutine setEdgeLocalBasis1

  !!
  !!
  !!
  pure subroutine setEdgeLocalBasis2(self, idx, localBasis2)
    class(edgeShelf), intent(inout)             :: self
    integer(shortInt), intent(in)               :: idx
    real(defReal), dimension(3), intent(in)     :: localBasis2

    call self % shelf(idx) % setLocalBasis2(localBasis2)

  end subroutine setEdgeLocalBasis2

  !!
  !!
  !!
  pure subroutine setEdgeLength(self, idx, Length)
    class(edgeShelf), intent(inout)             :: self
    integer(shortInt), intent(in)               :: idx
    real(defReal), intent(in)                   :: Length

    call self % shelf(idx) % setLength(Length)

  end subroutine setEdgeLength

  !!
  !!
  !!
  pure subroutine setEdgeDotProductOfVector(self, idx, dotProductOfVector)
    class(edgeShelf), intent(inout)             :: self
    integer(shortInt), intent(in)               :: idx
    real(defReal), intent(in)                   :: dotProductOfVector

    call self % shelf(idx) % setDotProductOfVector(dotProductOfVector)

  end subroutine setEdgeDotProductOfVector

  !!
  !!
  !!
  pure subroutine setEdgeAnglesArray(self, idx, anglesArray)
    class(edgeShelf), intent(inout)             :: self
    integer(shortInt), intent(in)               :: idx
    real(defReal), dimension(:, :), intent(in)  :: anglesArray

    call self % shelf(idx) % setAnglesArray(anglesArray)

  end subroutine setEdgeAnglesArray

  !!
  !!
  !!
  pure subroutine setEdgeElementIdxsArray(self, idx, elementIdxsArray)
    class(edgeShelf), intent(inout)                 :: self
    integer(shortInt), intent(in)                   :: idx
    integer(shortInt), dimension(:), intent(in)     :: elementIdxsArray

    call self % shelf(idx) % setElementIdxsArray(elementIdxsArray)

  end subroutine setEdgeElementIdxsArray

  !!
  !!
  !!
  elemental subroutine setEdgeIsBoundary(self, idx, isBoundary)
    class(edgeShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx
    logical(defBool), intent(in)    :: isBoundary

    call self % shelf(idx) % setIsBoundary(isBoundary)

   end subroutine setEdgeIsBoundary

  !!
  !!
  !!
  elemental function isAllocatedEdgeAnglesArray(self, idx) result(isAllocated)
    class(edgeShelf), intent(in)                    :: self
    integer(shortInt), intent(in)                   :: idx
    logical                                         :: isAllocated

    isAllocated = self % shelf(idx) % isAllocatedAnglesArray()

  end function isAllocatedEdgeAnglesArray


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