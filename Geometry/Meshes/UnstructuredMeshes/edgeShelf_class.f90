module edgeShelf_class
  
  use numPrecision
  use genericProcedures, only : findCommon, hasDuplicates
  use edge_class,        only : edge
  use faceShelf_class,   only : faceShelf
  
  implicit none
  private
  
  type, public :: edgeShelf
    private
    type(edge), dimension(:), allocatable :: shelf
  contains
    procedure                             :: getIdxsAndVertices
  end type edgeShelf

contains
  
  !!
  !!
  !!
  elemental subroutine getIdxsAndVertices(self, faces)
    class(edgeShelf), intent(inout)              :: self
    type(faceShelf), intent(in)                  :: faces
    integer(shortInt)                            :: nFaces, i, j, idx
    integer(shortInt), dimension(:), allocatable :: verticesFace1, verticesFace2, commonVertices
    type(edge)                                   :: edgeToAdd
    real(defReal), dimension(3)                  :: normalFace1, normalFace2
    
    idx = 1
    
    nFaces = size(faces % shelf)
    
    do i = 1, nFaces - 1
      verticesFace1 = faces % shelf(i) % getVertices()
      do j = i + 1, nFaces
        verticesFace2 = faces % shelf(j) % getVertices()
        if (hasDuplicates([verticesFace1, verticesFace2])) then
          commonVertices = findCommon(verticesFace1, verticesFace2)
        
          if (size(commonVertices) == 2) then
!            edgeToAdd % idx = idx
!            edgeToAdd % startVertex = minval(commonVertices)
!            edgeToAdd % endVertex = maxval(commonVertices)
        
            if (allocated(self % shelf)) then
          
              self % shelf = [self % shelf, edgeToAdd]
          
            else
          
              allocate(self % shelf(1))
              self % shelf(1) = edgeToAdd
          
            end if
          
            idx = idx + 1
        
          end if
      
        end if
      
      end do
    
    end do
  
  end subroutine getIdxsAndVertices
end module edgeShelf_class