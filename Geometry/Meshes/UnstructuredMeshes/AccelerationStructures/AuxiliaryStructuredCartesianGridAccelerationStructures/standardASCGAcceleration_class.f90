module standardASCGAcceleration_class

  use ASCGAcceleration_inter, only : ASCGAcceleration
  use CartesianCell_class,    only : CartesianCell
  use CartesianGrid_class,    only : CartesianGrid
  use edgeShelf_class,        only : edgeShelf
  use element_inter,          only : inclusionTestResult
  use elementShelf_class,     only : elementShelf
  use faceShelf_class,        only : faceShelf
  use numPrecision
  use vertexShelf_class,      only : vertexShelf
  use universalVariables,     only : INSIDE_ELEMENT

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(ASCGAcceleration) :: standardASCGAcceleration
    private
  contains
    procedure :: findHostElementIdx
  end type standardASCGAcceleration

contains
  !!
  !!
  !!
  subroutine findHostElementIdx(self, r, edges, elements, faces, vertices, elementIdx)
    class(standardASCGAcceleration), intent(in)  :: self
    real(defReal), dimension(3), intent(in)      :: r
    type(edgeShelf), intent(in)                  :: edges
    type(elementShelf), intent(in)               :: elements
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    integer(shortInt), intent(inout)             :: elementIdx
    integer(shortInt)                            :: i
    integer(shortInt), dimension(:), allocatable :: intersectedFaceIdxs, potentialElementIdxs
    type(CartesianCell), pointer                 :: terminalCellPtr
    type(inclusionTestResult)                    :: elementInclusionResults

    ! Get pointer to deepest grid.
    terminalCellPtr => self % searchGrids(r)

    ! Check if terminal cell is fully inside an element and return immediately if so.
    elementIdx = terminalCellPtr % getElementIdx()
    if(0 < elementIdx .or. terminalCellPtr % isOutside()) return

    intersectedFaceIdxs = terminalCellPtr % getIntersectedFaceIdxs()
    potentialElementIdxs = faces % getFaceElementIdxs(intersectedFaceIdxs)

    ! For multi-layered Patch-Search, check if terminal cell only intersects with a single face. In this case, perform an 
    ! element inclusion test on the elements sharing this face and return.
    if(1 < self % getDepth()) then
      if(size(intersectedFaceIdxs) == 1) then
        if(ZERO < dot_product(faces % getFaceCentroid(intersectedFaceIdxs(1)) - r, &
           faces % getFaceNormal(intersectedFaceIdxs(1)))) then
          elementIdx = minval(potentialElementIdxs)

        elseif(.not. faces % getFaceIsBoundary(intersectedFaceIdxs(1))) then
          elementIdx = maxval(potentialElementIdxs)

        end if
        return

      end if

    end if

    ! Else, search all potential elements.
    do i = 1, size(potentialElementIdxs)
      elementInclusionResults = elements % isPointInside(potentialElementIdxs(i), r, faces)
      if(elementInclusionResults % status == INSIDE_ELEMENT) then
        elementIdx = potentialElementIdxs(i)
        return

      end if

    end do

  end subroutine findHostElementIdx

end module standardASCGAcceleration_class