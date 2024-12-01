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
    procedure                                       :: getSize
    procedure                                       :: init
    procedure                                       :: kill
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
  
  !! Subroutine 'init'
  !!
  !! Basic description:
  !!   Initialises the shelf from the 'points' file.
  !!
  !! Arguments:
  !!   folderPath [in] -> Path to the folder containing the mesh files.
  !!   nVertices [in]  -> Number of vertices in the mesh.
  !!
  subroutine init(self, folderPath, nVertices)
    class(vertexShelf), intent(inout) :: self
    character(*), intent(in)          :: folderPath
    integer(shortInt), intent(in)     :: nVertices
    integer(shortInt)                 :: i, j, leftBracketIdx, rightBracketIdx
    integer(shortInt), parameter      :: unit = 10
    real(defReal), dimension(3)       :: coordinates
    logical(defBool)                  :: singleLine
    character(150)                    :: string ! Note: here the string is longer than usual to deal
                                                ! with cases when all vertices are written on a
                                                ! single line.

    ! Open the 'points' file and read it until a line containing the symbol ')' is encountered.
    call openToRead(unit, folderPath//'points')
    read(unit, "(a)") string
    do while (index(string(1:len_trim(string)), ")") == 0)
      read(unit, "(a)") string

    end do

    ! Initialise singleLine = .false. and check if all vertices are written on a single line.
    ! If yes, update singleLine and erase the leftmost bracket from the string.
    singleLine = .false.
    if (string(2:2) == '(') then
      singleLine = .true.
      string(1:1) = ''

    end if

    ! Loop over all vertices.
    do i = 1, nVertices
      if (singleLine) then
        ! Locate the leftmost and rightmost brackets. Read the coordinates between the
        ! brackets and remove them from the string.
        leftBracketIdx = index(string, '(')
        rightBracketIdx = index(string, ')')
        read(string(leftBracketIdx + 1:rightBracketIdx - 1), *) coordinates
        string(leftBracketIdx:leftBracketIdx) = ''
        string(rightBracketIdx:rightBracketIdx) = ''

      else
        ! Read the coordinates of the current vertex and move onto the next line.
        read(string(2:len_trim(string) - 1), *) coordinates
        read(unit, "(a)") string

      end if

      ! Set the index and coordinates of the current vertex. Apply offset in the process.
      call self % shelf(i) % setIdx(i)
      coordinates = coordinates + self % offset
      call self % shelf(i) % setCoordinates(coordinates)

      ! Update extremal coordinates in the shelf if necessary.
      do j = 1, 3
        if (coordinates(j) < self % extremalCoordinates(j)) self % extremalCoordinates(j) = coordinates(j)
        if (coordinates(j) > self % extremalCoordinates(j + 3)) self % extremalCoordinates(j + 3) = coordinates(j)

      end do

    end do

    ! Close the 'points' file.
    close(unit)
  end subroutine init
  
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
end module vertexShelf_class