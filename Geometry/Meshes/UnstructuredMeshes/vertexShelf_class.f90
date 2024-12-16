module vertexShelf_class
  
  use numPrecision
  use genericProcedures, only : findCommon
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
    type(vertex), dimension(:), allocatable, public :: shelf
    real(defReal), dimension(3)                     :: offset = ZERO
    real(defReal), dimension(6)                     :: extremalCoordinates = ZERO
  contains
    procedure                                       :: findCommonEdgeIdx
    procedure                                       :: getAllCoordinates
    procedure                                       :: getExtremalCoordinates
    procedure                                       :: getOffset
    procedure                                       :: getSize
    procedure                                       :: kill
    procedure                                       :: setExtremalCoordinates
    procedure                                       :: setOffset
  end type

contains
  
  !! Subroutine 'applyOffset'
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

    ! Initialise edgeIdx = 0 and return immediately if either vertex are not associated with edges.
    edgeIdx = 0
    if (.not. self % shelf(firstVertexIdx) % hasEdges() .or. .not. self % shelf(secondVertexIdx) % hasEdges()) return
    
    ! Find common edge indices. Update edgeIdx only if common indices have been found.
    commonIdxs = findCommon(self % shelf(firstVertexIdx) % getEdgeIdxs(), self % shelf(secondVertexIdx) % getEdgeIdxs())
    if (size(commonIdxs) > 0) edgeIdx = commonIdxs(1)

  end function findCommonEdgeIdx

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

  pure subroutine setExtremalCoordinates(self, coords)
    class(vertexShelf), intent(inout)       :: self
    real(defReal), dimension(6), intent(in) :: coords

    self % extremalCoordinates = coords

  end subroutine setExtremalCoordinates

end module vertexShelf_class