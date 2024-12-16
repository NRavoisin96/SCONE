module edgeShelf_class
  
  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use edge_class,        only : edge
  
  implicit none
  private
  
  type, public                                    :: edgeShelf
    private
    type(edge), dimension(:), allocatable, public :: shelf
  contains
    procedure                                     :: collapse
    procedure                                     :: kill
  end type edgeShelf

contains

  !! Subroutine 'collapse'
  !!
  !! Basic description:
  !!   Reduces the size of the shelf to lastIdx.
  !!
  !! Arguments:
  !!   lastIdx [in] -> Index of the last edge in the shelf to be collapsed.
  !!
  !! Errors:
  !!   - fatalError if lastIdx < 1.
  !!   - fatalError if lastIdx > size(self % shelf).
  !!
  subroutine collapse(self, lastIdx)
    class(edgeShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: lastIdx
    type(edgeShelf)                 :: tempShelf
    character(*), parameter         :: Here = 'collapse (edgeShelf_class.f90)'
    
    ! Catch invalid idx.
    if (lastIdx < 1) call fatalError(Here, 'edgeShelf size must be +ve. Is: '//numToChar(lastIdx)//'.')
    if (lastIdx > size(self % shelf)) call fatalError(Here, 'New shelf size must be smaller than original size.')
    
    ! Create a temporary shelf and copy all the elements up to lastVertexIdx from the 
    ! original shelf into the temporary one.
    allocate(tempShelf % shelf(lastIdx))
    tempShelf % shelf = self % shelf(1:lastIdx)
    
    ! kill the original shelf and reallocate memory.
    call self % kill()
    allocate(self % shelf(lastIdx))
    
    ! Copy back and kill the temporary shelf.
    self % shelf = tempShelf % shelf
    call tempShelf % kill()

  end subroutine collapse

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an unitialised state.
  !!
  elemental subroutine kill(self)
    class(edgeShelf), intent(inout) :: self

    if (allocated(self % shelf)) deallocate(self % shelf)

  end subroutine kill
  
end module edgeShelf_class