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
  subroutine findHostElementIdx(self, u, edges, elements, faces, vertices, elementIdx, r)
    class(standardASCGAcceleration), intent(in)  :: self
    real(defReal), dimension(3), intent(in)      :: u
    type(edgeShelf), intent(in)                  :: edges
    type(elementShelf), intent(in)               :: elements
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    integer(shortInt), intent(inout)             :: elementIdx
    real(defReal), dimension(3), intent(inout)   :: r
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
    if(1 < self % getDepth()) then
      if(size(intersectedFaceIdxs) == 1) then
        call faces % testFaceHalfSpace(intersectedFaceIdxs(1), r, elementIdx)
        return

      end if

    end if

    ! Else, search all potential elements.
    potentialElementIdxs = faces % getFaceElementIdxs(intersectedFaceIdxs)
    do i = 1, size(potentialElementIdxs)
      elementInclusionResults = elements % isPointInside(potentialElementIdxs(i), r, faces)
      if(elementInclusionResults % status == INSIDE_ELEMENT) then
        elementIdx = potentialElementIdxs(i)
        return

      end if

    end do

  end subroutine findHostElementIdx

end module standardASCGAcceleration_class