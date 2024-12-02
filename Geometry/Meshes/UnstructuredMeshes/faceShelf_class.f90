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
    procedure                                     :: kill
  end type

contains
  
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