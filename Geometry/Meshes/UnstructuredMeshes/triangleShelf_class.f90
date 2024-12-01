module triangleShelf_class
  
  use numPrecision
  use triangle_class, only : triangle
  
  implicit none
  private
  
  !!
  !! Storage space for triangles in a given OpenFOAM mesh.
  !!
  !! Private members:
  !!   shelf -> An array of triangles.
  !!
  type, public :: triangleShelf
    private
    type(triangle), dimension(:), allocatable, public :: shelf
  contains
    procedure                                         :: kill
  end type triangleShelf

contains
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(triangleShelf), intent(inout) :: self
    
    if (allocated(self % shelf)) deallocate(self % shelf)
  end subroutine kill
end module triangleShelf_class