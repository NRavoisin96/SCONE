module elementShelf_class
  
  use numPrecision
  use element_class, only : element
  
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
    procedure                                        :: kill
  end type elementShelf

contains
  
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