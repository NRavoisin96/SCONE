module elementShelf_class
  
  use numPrecision
  use genericProcedures, only : append, findCommon, hasDuplicates, openToRead, removeDuplicates
  use element_class,     only : element
  use face_class,        only : face
  use faceShelf_class,   only : faceShelf
  use vertexShelf_class, only : vertexShelf
  
  implicit none
  private
  
  !!
  !! Storage space for elements in a given OpenFOAM mesh.
  !!
  !! Private members:
  !!   shelf -> Array to store elements.
  !!
  type, public :: elementShelf
    private
    type(element), dimension(:), allocatable, public :: shelf
  contains
    procedure                                        :: init
    procedure                                        :: kill
  end type elementShelf

contains
  
  !! Subroutine 'init'
  !!
  !! Basic description:
  !!   Initialises the shelf from the 'owner' and 'neighbour' files.
  !!
  !! Arguments:
  !!   folderPath [in]     -> Path of the folder containing the OpenFOAM mesh files.
  !!   nFaces [in]         -> Number of faces in the mesh.
  !!   nInternalFaces [in] -> Number of internal faces in the mesh.
  !!
  subroutine init(self, folderPath, nElements, nFaces, nInternalFaces)
    class(elementShelf), intent(inout)           :: self
    character(*), intent(in)                     :: folderPath
    integer(shortInt), intent(in)                :: nElements, nFaces, nInternalFaces
    integer(shortInt)                            :: i, elementIdx
    integer(shortInt), parameter                 :: unit = 10
    integer(shortInt), dimension(:), allocatable :: elementIdxs
    character(100)                               :: string

    ! Retrieve the number of elements and set the index of each element.
    do i = 1, nElements
      call self % shelf(i) % setIdx(i)

    end do
    
    ! If there is only one element in the mesh simply add all the faces to this element and return.
    if (nElements == 1) then
      do i = 1, nFaces
        call self % shelf(1) % addFaceToElement(i)

      end do
      return

    end if

    ! Open the 'owner' file and read it until a line containing the symbol '(' is encountered.
    call openToRead(unit, folderPath//'owner')
    read(unit, "(a)") string
    do while (index(string(1:len_trim(string)), "(") == 0)
      read(unit, "(a)") string

    end do
    
    ! Skip one more line and loop over all faces.
    read(unit, "(a)") string
    do i = 1, nFaces
      ! Read the current element index and add the current face to this element. Note
      ! that we add one since Fortran starts indexing at one and not zero.
      read(string, *) elementIdx
      call self % shelf(elementIdx + 1) % addFaceToElement(i)
      ! Move onto the next line.
      read(unit, "(a)") string

    end do

    ! Close the 'owner' file.
    close(unit)
    
    ! If nInternalFaces = 0 we can return early here. Else we need to repeat the above procedure
    ! for the 'neighbour file.
    if (nInternalFaces == 0) return
    
    ! Open the 'neighbour' file and read it until a line containing the symbol '(' is encountered.
    call openToRead(unit, folderPath//'neighbour')
    read(unit, "(a)") string
    do while (index(string(1:len_trim(string)), "(") == 0)
      read(unit, "(a)") string

    end do
    
    ! Check if the current line contains the symbol ')'. If it does, then all element indices
    ! are written on a single line.
    if (index(string(1:len_trim(string)), ")") > 0) then
      ! Allocate the number of entries in the 'elementIndices' array to the number of internal
      ! faces and copy element indices into this array.
      allocate(elementIdxs(nInternalFaces))
      read(string(index(string(1:len_trim(string)), "(") + 1:&
      index(string(1:len_trim(string)), ")") - 1), *) elementIdxs
      
      ! Loop over all internal faces.
      do i = 1, nInternalFaces
        ! Update mesh connectivity information.
        call self % shelf(elementIdxs(i) + 1) % addFaceToElement(-i)

      end do
    ! If not element indices are spread over multiple lines with one element index per line.
    else
      ! Skip one more line and loop over all internal faces.
      read(unit, "(a)") string
      do i = 1, nInternalFaces
        ! Read the current element index and add the current internal face to this element.
        read(string, *) elementIdx
        ! Update mesh connectivity information.
        call self % shelf(elementIdx + 1) % addFaceToElement(-i)
        
        ! Move onto the next line.
        read(unit, "(a)") string

      end do

    end if
    
    ! Close the 'neighbour' file.
    close(unit)
  end subroutine init
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(elementShelf), intent(inout) :: self
    
    if (allocated(self % shelf)) deallocate(self % shelf)

  end subroutine kill
  
end module elementShelf_class