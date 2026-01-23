module terminalCartesianCell_class

  use CartesianCell_inter,        only : CartesianCell
  use cartesianGenericProcedures, only : testIntervalIntersection
  use edgeShelf_class,            only : edgeShelf
  use elementShelf_class,         only : elementShelf
  use face_inter,                 only : faceSATData
  use faceShelf_class,            only : faceShelf
  use numPrecision
  use vertexShelf_class,          only : vertexShelf

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(CartesianCell) :: terminalCartesianCell
    private
    integer(shortInt)               :: edgeIdx = 0, elementIdx = 0, gridIdx = 0, vertexIdx = 0
    integer(shortInt), dimension(2) :: intersectedFaceIdxs = 0
  contains
    procedure :: getEdgeIdx
    procedure :: getElementIdx
    procedure :: getGridIdx
    procedure :: getVertexIdx
    procedure :: isOutside
    procedure :: isUnprocessed
    procedure :: kill
    procedure :: setGridIdx
    procedure :: testElementInclusion
    procedure :: testFaceIntersection
  end type terminalCartesianCell

contains
  !!
  !!
  !!
  elemental function getEdgeIdx(self) result(edgeIdx)
    class(terminalCartesianCell), intent(in) :: self
    integer(shortInt)                        :: edgeIdx

    edgeIdx = self % edgeIdx

  end function getEdgeIdx

  !!
  !!
  !!
  elemental function getElementIdx(self) result(elementIdx)
    class(terminalCartesianCell), intent(in) :: self
    integer(shortInt)                        :: elementIdx

    elementIdx = self % elementIdx

  end function getElementIdx

  !!
  !!
  !!
  elemental function getGridIdx(self) result(gridIdx)
    class(terminalCartesianCell), intent(in) :: self
    integer(shortInt)                        :: gridIdx

    gridIdx = self % gridIdx

  end function getGridIdx

  !!
  !!
  !!
  elemental function getVertexIdx(self) result(vertexIdx)
    class(terminalCartesianCell), intent(in) :: self
    integer(shortInt)                        :: vertexIdx

    vertexIdx = self % vertexIdx

  end function getVertexIdx

  !!
  !!
  !!
  elemental function isOutside(self) result(isIt)
    class(terminalCartesianCell), intent(in) :: self
    logical(defBool)                         :: isIt

    isIt = self % edgeIdx == 0 .and. self % elementIdx == 0 .and. self % vertexIdx == 0

  end function isOutside

  !!
  !!
  !!
  elemental function isUnprocessed(self) result(isIt)
    class(terminalCartesianCell), intent(in) :: self
    logical(defBool)                         :: isIt

    isIt = self % elementIdx == 0 .and. all(self % intersectedFaceIdxs == 0)

  end function isUnprocessed

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(terminalCartesianCell), intent(inout) :: self

    ! Local.
    self % edgeIdx = 0
    self % elementIdx = 0
    self % gridIdx = 0
    self % vertexIdx = 0
    self % intersectedFaceIdxs = 0

  end subroutine kill

  !!
  !!
  !!
  elemental subroutine setGridIdx(self, gridIdx)
    class(terminalCartesianCell), intent(in) :: self
    integer(shortInt), intent(in)            :: gridIdx

    self % gridIdx = gridIdx

  end subroutine setGridIdx

  !!
  !!
  !!
  subroutine testElementInclusion(self, elementIdx, elementFaceIdxs, centroid, elements, faces)
    class(terminalCartesianCell), intent(inout) :: self
    integer(shortInt), intent(in)               :: elementIdx
    integer(shortInt), dimension(:), intent(in) :: elementFaceIdxs
    real(defReal), dimension(3), intent(in)     :: centroid
    type(elementShelf), intent(in)              :: elements
    type(faceShelf), intent(in)                 :: faces
    integer(shortInt)                           :: faceIdx, i

    ! Test centroid inclusion within element.
    do i = 1, size(elementFaceIdxs)
      faceIdx = elementFaceIdxs(i)
      if(dot_product(faces % getFaceCentroid(abs(faceIdx)) - centroid, faces % getFaceNormal(faceIdx)) < ZERO) return

    end do
    self % elementIdx = elementIdx

  end subroutine testElementInclusion

  !!
  !!
  !!
  subroutine testFaceIntersection(self, faceIdx, extraDistance, spacing, targetDistance, centroid, edges, cache, faces, vertices)
    class(terminalCartesianCell), intent(inout)  :: self
    integer(shortInt), intent(in)                :: faceIdx
    real(defReal), intent(in)                    :: extraDistance, spacing, targetDistance
    real(defReal), dimension(3), intent(in)      :: centroid
    type(edgeShelf), intent(in)                  :: edges
    type(faceSATData), intent(in)                :: cache
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    integer(shortInt)                            :: commonEdgeIdx, i, j
    integer(shortInt), dimension(2)              :: commonEdgeVertexIdxs
    integer(shortInt), dimension(:), allocatable :: faceEdgeIdxs
    real(defReal)                                :: distanceToVertex1Squared, distanceToVertex2Squared, projection
    real(defReal), dimension(3)                  :: commonEdgeUnitVector, commonEdgeVertex1Coords, projections, radii, temp

    ! Return immediately if mapping has already been assigned for this cell.
    if(0 < self % edgeIdx .and. 0 < self % vertexIdx) return

    ! Test along face normal.
    projection = dot_product(centroid, cache % faceNormal)
    if(.not. testIntervalIntersection(projection - extraDistance, projection + extraDistance, -cache % faceConstant, &
                                      -cache % faceConstant)) return

    ! Now test along edge axes using cached intervals.
    do i = 1, size(cache % edgeVectors, 2)
      ! Compute centroid projections and radii along each axis.
      projections(1) = cache % edgeVectors(3, i) * centroid(2) - cache % edgeVectors(2, i) * centroid(3)
      projections(2) = cache % edgeVectors(1, i) * centroid(3) - cache % edgeVectors(3, i) * centroid(1)
      projections(3) = cache % edgeVectors(2, i) * centroid(1) - cache % edgeVectors(1, i) * centroid(2)
      radii(1) = HALF * (abs(cache % edgeVectors(2, i)) + abs(cache % edgeVectors(3, i))) * spacing
      radii(2) = HALF * (abs(cache % edgeVectors(1, i)) + abs(cache % edgeVectors(3, i))) * spacing
      radii(3) = HALF * (abs(cache % edgeVectors(1, i)) + abs(cache % edgeVectors(2, i))) * spacing

      do j = 1, 3
        if(.not. testIntervalIntersection(projections(j) - radii(j), projections(j) + radii(j), &
                                          cache % edgeAxesIntervals(j, i, 1), cache % edgeAxesIntervals(j, i, 2))) return

      end do

    end do

    ! If reached here, the cell intersects the face so append the index to the list of intersected faces.
    if(self % intersectedFaceIdxs(1) == 0) then
      self % intersectedFaceIdxs(1) = faceIdx
      
      ! Assign edge mapping to first edge in the face.
      faceEdgeIdxs = faces % getFaceEdgeIdxs(faceIdx)
      self % edgeIdx = faceEdgeIdxs(1)

    else
      self % intersectedFaceIdxs(2) = faceIdx

      ! Check if the cell intersects the common edge between the two faces.
      commonEdgeIdx = faces % findCommonEdgeIdx(self % intersectedFaceIdxs(1), self % intersectedFaceIdxs(2))

      ! if there is no common edge, exit the subroutine early (other combinations of faces to be tried later)
      if (commonEdgeIdx == 0) return
      self % edgeIdx = 0

      ! calculate centre of intersection between the circumscribed ball and the plane parallel to the polygon
      commonEdgeVertexIdxs = edges % getEdgeVertexIdxs(commonEdgeIdx)
      commonEdgeUnitVector = edges % getEdgeUnitVector(commonEdgeIdx)
      commonEdgeVertex1Coords = vertices % getVertexCoordinates(commonEdgeVertexIdxs(1))

      ! Project centroid onto common edge then compute squared distance to vertices 1 and 2.
      projection = dot_product(centroid - commonEdgeVertex1Coords, commonEdgeUnitVector)
      distanceToVertex1Squared = projection * projection
      temp = commonEdgeVertex1Coords + projection * commonEdgeUnitVector - &
             vertices % getVertexCoordinates(commonEdgeVertexIdxs(2))
      distanceToVertex2Squared = dot_product(temp, temp)

      ! make comparison between the two distances (dotProduct1/2)
      self % vertexIdx = merge(commonEdgeVertexIdxs(1), commonEdgeVertexIdxs(2), &
                               distanceToVertex1Squared < distanceToVertex2Squared)
      if(targetDistance < distanceToVertex1Squared .and. targetDistance < distanceToVertex2Squared) &
      self % edgeIdx = commonEdgeIdx

    end if

  end subroutine testFaceIntersection

end module terminalCartesianCell_class