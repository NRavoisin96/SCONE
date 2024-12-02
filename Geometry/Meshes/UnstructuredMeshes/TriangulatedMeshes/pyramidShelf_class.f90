module pyramidShelf_class
  
  use numPrecision
  use pyramid_class, only : pyramid
  
  implicit none
  private
  
  !!
  !! Storage space for pyramids in a given OpenFOAM mesh.
  !!
  !! Private members:
  !!   shelf -> Array of pyramids.
  !!
  type, public :: pyramidShelf
    private
    type(pyramid), dimension(:), allocatable, public :: shelf
  contains
    procedure                                        :: kill
  end type pyramidShelf

contains
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(pyramidShelf), intent(inout) :: self
    
    if (allocated(self % shelf)) deallocate(self % shelf)

  end subroutine kill
  
end module pyramidShelf_class