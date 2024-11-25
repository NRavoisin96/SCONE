module faceShelf_class
  
  use numPrecision
  use genericProcedures,   only : openToRead
  use face_class,          only : face
  use triangleShelf_class, only : triangleShelf
  use vertex_class,        only : vertex
  
  implicit none
  private
  
  !!
  !! Storage space for faces of an OpenFOAM mesh.
  !!
  !! Private members:
  !!   shelf -> Array to store faces.
  !!
  type, public :: faceShelf
    private
    type(face), dimension(:), allocatable, public :: shelf
  contains
    procedure                                     :: getSize
    procedure                                     :: init
    procedure                                     :: kill
  end type

contains
  
  !! Subroutine 'init'
  !!
  !! Basic description:
  !!   Initialises the shelf from the 'faces' file.
  !!
  !! Detailed description:
  !!   Opens the 'faces' file and reads it until a line with the symbol ')' is encountered, which 
  !!   marks the beginning of the faces' data listing. In the 'faces' file, each line corresponds to
  !!   a single face. The subroutine then retrieves the number of vertices in each face and their
  !!   indices sequentially line-by-line.
  !!
  !! Arguments:
  !!   folderPath [in]     -> Path of the folder containing the mesh files.
  !!   nInternalFaces [in] -> Number of internal faces in the mesh.
  !!
  subroutine init(self, folderPath, nInternalFaces)
    class(faceShelf), intent(inout)              :: self
    character(*), intent(in)                     :: folderPath
    integer(shortInt), intent(in)                :: nInternalFaces
    integer(shortInt), parameter                 :: unit = 10
    integer(shortInt)                            :: i, j, nVertices
    integer(shortInt), dimension(:), allocatable :: vertexIdxs
    character(100)                               :: string

    ! Open the 'faces' data file and read it until a line containing the symbol ')' is encountered.
    call openToRead(unit, folderPath//'faces')
    read(unit, "(a)") string
    do while (index(string(1:len_trim(string)), ")") == 0)
      read(unit, "(a)") string

    end do
    
    ! Loop through all faces in the shelf.
    do i = 1, size(self % shelf)
      ! Set the face index. If i > nInternalFaces set this face as a boundary face.
      call self % shelf(i) % setIdx(i)
      if (i > nInternalFaces) call self % shelf(i) % setBoundaryFace()
      
      ! Skip lines if they are blank (can happen for faces with a large number of vertices).
      if (len_trim(string) == 0) read(unit, "(a)") string
      
      ! Check if all the vertices are listed on a single line and retrieve nVertices accordingly.
      if (index(string(1:len_trim(string)), "(") > 0) then
        read(string(1:index(string, "(") - 1), *) nVertices

      else
        read(string, *) nVertices
        read(unit, "(a)") string

      end if

      ! Allocate memory.
      allocate(vertexIdxs(nVertices))
      
      ! If all vertices are listed on a single line copy them directly.
      if (index(string(1:len_trim(string)), ")") > 0) then
        read(string(index(string, "(") + 1:len_trim(string) - 1), *) vertexIdxs

      else
        ! Skip one line and read each vertex index line-by-line.
        read(unit, "(a)") string
        do j = 1, nVertices
          read(string, *) vertexIdxs(j)
          read(unit, "(a)") string

        end do
        
        ! Skip one line.
        read(unit, "(a)") string

      end if
      
      ! Add one to the vertices indices since Fortran starts indexing at one rather than zero and
      ! add each vertex to the current face.
      vertexIdxs = vertexIdxs + 1
      do j = 1, size(vertexIdxs)
        call self % shelf(i) % addVertexToFace(vertexIdxs(j))

      end do
      
      ! Free memory and move onto the next line.
      deallocate(vertexIdxs)
      read(unit, "(a)") string

    end do
    
    ! Close the 'faces' file.
    close(unit)
  end subroutine init
  
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
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(faceShelf), intent(inout) :: self
    
    if (allocated(self % shelf)) deallocate(self % shelf)

  end subroutine kill
  
end module faceShelf_class