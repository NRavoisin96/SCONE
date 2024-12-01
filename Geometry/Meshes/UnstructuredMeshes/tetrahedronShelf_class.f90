module tetrahedronShelf_class
  
  use numPrecision
  use tetrahedron_class, only : tetrahedron
  
  implicit none
  private
  
  !!
  !! Storage space for tetrahedra in a given OpenFOAM mesh.
  !!
  !! Private members:
  !!   shelf -> An array of tetrahedra.
  !!
  type, public :: tetrahedronShelf
    private
    type(tetrahedron), dimension(:), allocatable, public :: shelf
  contains
    procedure                                            :: kill
  end type tetrahedronShelf

contains
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(tetrahedronShelf), intent(inout) :: self
    
    if (allocated(self % shelf)) deallocate(self % shelf)
  end subroutine kill
end module tetrahedronShelf_class