module vertexShelf_class
  
  use numPrecision
  use genericProcedures, only : fatalError, numToChar, openToRead
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
    procedure                                       :: setOffset
    procedure                                       :: collapse
    procedure                                       :: getAllCoordinates
    procedure                                       :: getExtremalCoordinates
    procedure                                       :: getOffset
    procedure                                       :: getSize
    procedure                                       :: kill
    procedure                                       :: setExtremalCoordinates
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
  
  !! Subroutine 'collapse'
  !!
  !! Basic description:
  !!   Reduces the size of the shelf to lastVertexIdx.
  !!
  !! Arguments:
  !!   lastVertexIdx [in] -> Index of the last vertex in the collapsed shelf.
  !!
  !! Errors:
  !!   fatalError if lastVertexIdx < 1.
  !!
  subroutine collapse(self, lastVertexIdx)
    class(vertexShelf), intent(inout) :: self
    integer(shortInt), intent(in)     :: lastVertexIdx
    type(vertexShelf)                 :: tempShelf
    character(100), parameter         :: Here = 'collapse (vertexShelf_class.f90)'
    
    ! Catch invalid idx.
    if (lastVertexIdx < 1) call fatalError(Here, 'vertexShelf size must be +ve. Is: '//numToChar(lastVertexIdx)//'.')
    
    ! Create a temporary shelf and copy all the elements up to lastVertexIdx from the 
    ! original shelf into the temporary one.
    allocate(tempShelf % shelf(lastVertexIdx))
    tempShelf % shelf = self % shelf(1:lastVertexIdx)
    
    ! kill the original shelf and reallocate memory.
    call self % kill()
    allocate(self % shelf(lastVertexIdx))
    
    ! Copy back and kill the temporary shelf.
    self % shelf = tempShelf % shelf
    call tempShelf % kill()
  
  end subroutine collapse

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