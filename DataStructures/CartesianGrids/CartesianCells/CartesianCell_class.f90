module CartesianCell_class

  use cartesianGenericProcedures, only : testIntervalIntersection
  use edgeShelf_class,            only : edgeShelf
  use elementShelf_class,         only : elementShelf
  use errors_mod,                 only : fatalError
  use face_inter,                 only : faceSATData
  use faceShelf_class,            only : faceShelf
  use genericProcedures,          only : append
  use numPrecision
  use vertexShelf_class,          only : vertexShelf

  implicit none
  private

  !!
  !!
  !!
  type, public :: CartesianCell
    private
    integer(shortInt)                            :: edgeIdx = 0, elementIdx = 0, subGridIdx = 0, vertexIdx = 0
    integer(shortInt), dimension(:), allocatable :: intersectedFaceIdxs
  contains
    procedure :: getEdgeIdx
    procedure :: getElementIdx
    procedure :: getIntersectedFaceIdxs
    procedure :: getSubGridIdx
    procedure :: getVertexIdx
    procedure :: isOutside
    procedure :: isSimple
    procedure :: isUnprocessed
    procedure :: kill
    procedure :: setSubGridIdx
    procedure :: testElementInclusion
    procedure :: testFaceIntersection
  end type CartesianCell

contains
  !!
  !!
  !!
  elemental function getEdgeIdx(self) result(edgeIdx)
    class(CartesianCell), intent(in) :: self
    integer(shortInt)                :: edgeIdx

    edgeIdx = self % edgeIdx

  end function getEdgeIdx

  !!
  !!
  !!
  elemental function getElementIdx(self) result(elementIdx)
    class(CartesianCell), intent(in) :: self
    integer(shortInt)                :: elementIdx

    elementIdx = self % elementIdx

  end function getElementIdx

  !!
  !!
  !!
  pure function getIntersectedFaceIdxs(self) result(intersectedFaceIdxs)
    class(CartesianCell), intent(in)             :: self
    integer(shortInt), dimension(:), allocatable :: intersectedFaceIdxs

    if(allocated(self % intersectedFaceIdxs)) then
      intersectedFaceIdxs = self % intersectedFaceIdxs

    else
      allocate(intersectedFaceIdxs(0))

    end if

  end function

  !!
  !!
  !!
  elemental function getSubGridIdx(self) result(subGridIdx)
    class(CartesianCell), intent(in) :: self
    integer(shortInt)                :: subGridIdx

    subGridIdx = self % subGridIdx

  end function getSubGridIdx

  !!
  !!
  !!
  elemental function getVertexIdx(self) result(vertexIdx)
    class(CartesianCell), intent(in) :: self
    integer(shortInt)                :: vertexIdx

    vertexIdx = self % vertexIdx

  end function getVertexIdx

  !!
  !!
  !!
  elemental function isOutside(self) result(isIt)
    class(CartesianCell), intent(in) :: self
    logical(defBool)                 :: isIt

    isIt = self % edgeIdx == 0 .and. self % elementIdx == 0 .and. self % vertexIdx == 0

  end function isOutside

  !!
  !!
  !!
  elemental function isSimple(self) result(isIt)
    class(CartesianCell), intent(in) :: self
    logical(defBool)                 :: isIt

    isIt = 0 < self % elementIdx .or. self % isOutside() .or. size(self % getIntersectedFaceIdxs()) < 2

  end function isSimple

  !!
  !!
  !!
  elemental function isUnprocessed(self) result(isIt)
    class(CartesianCell), intent(in) :: self
    logical(defBool)                 :: isIt

    isIt = self % elementIdx == 0 .and. .not. allocated(self % intersectedFaceIdxs)

  end function isUnprocessed

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(CartesianCell), intent(inout) :: self

    ! Local.
    self % edgeIdx = 0
    self % elementIdx = 0
    self % subGridIdx = 0
    self % vertexIdx = 0
    if(allocated(self % intersectedFaceIdxs)) deallocate(self % intersectedFaceIdxs)

  end subroutine kill

  !!
  !!
  !!
  elemental subroutine setSubGridIdx(self, subGridIdx)
    class(CartesianCell), intent(inout) :: self
    integer(shortInt), intent(in)       :: subGridIdx

    self % subGridIdx = subGridIdx

  end subroutine setSubGridIdx

  !!
  !!
  !!
  pure subroutine testElementInclusion(self, elementIdx, elementFaceIdxs, centroid, elements, faces)
    class(CartesianCell), intent(inout)         :: self
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
  subroutine testFaceIntersection(self, faceIdx, extraDistance, spacing, targetDistance, centroid, edges, cache, faces, &
                                       vertices)
    class(CartesianCell), intent(inout)          :: self
    integer(shortInt), intent(in)                :: faceIdx
    real(defReal), intent(in)                    :: extraDistance, spacing, targetDistance
    real(defReal), dimension(3), intent(in)      :: centroid
    type(edgeShelf), intent(in)                  :: edges
    type(faceSATData), intent(in)                :: cache
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    integer(shortInt)                            :: commonEdgeIdx, i, j, nIntersectedFaces
    integer(shortInt), dimension(2)              :: commonEdgeVertexIdxs
    integer(shortInt), dimension(:), allocatable :: faceEdgeIdxs, tempArray
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
    call append(self % intersectedFaceIdxs, faceIdx)
    nIntersectedFaces = size(self % intersectedFaceIdxs)
    if(nIntersectedFaces == 1) then
      ! Assign edge mapping to first edge in the face.
      faceEdgeIdxs = faces % getFaceEdgeIdxs(faceIdx)
      self % edgeIdx = faceEdgeIdxs(1)

    else
      ! Loop through all pairs of faces and try to find a common edge which intersects the cell.
      outerLoop: do i = 1, nIntersectedFaces - 1
        do j = i + 1, nIntersectedFaces
          commonEdgeIdx = faces % findCommonEdgeIdx(self % intersectedFaceIdxs(i), self % intersectedFaceIdxs(j))
          if(0 < commonEdgeIdx) then
            self % edgeIdx = 0
            ! Retrieve vertex indices, unit vector and first vertex coordinates for the common edge.
            commonEdgeVertexIdxs = edges % getEdgeVertexIdxs(commonEdgeIdx)
            commonEdgeUnitVector = edges % getEdgeUnitVector(commonEdgeIdx)
            commonEdgeVertex1Coords = vertices % getVertexCoordinates(commonEdgeVertexIdxs(1))

            ! Project centroid onto common edge then compute squared distance to vertices 1 and 2.
            projection = dot_product(centroid - commonEdgeVertex1Coords, commonEdgeUnitVector)
            distanceToVertex1Squared = projection * projection
            temp = commonEdgeVertex1Coords + projection * commonEdgeUnitVector - &
                   vertices % getVertexCoordinates(commonEdgeVertexIdxs(2))
            distanceToVertex2Squared = dot_product(temp, temp)

            ! Assign vertexIdx for this cell.
            self % vertexIdx = merge(commonEdgeVertexIdxs(1), commonEdgeVertexIdxs(2), &
                                     distanceToVertex1Squared < distanceToVertex2Squared)

            ! Check if cell is far enough from both edge vertices.
            if(targetDistance < distanceToVertex1Squared .and. targetDistance < distanceToVertex2Squared) then
              self % edgeIdx = commonEdgeIdx
              exit outerLoop

            end if

          end if

        end do

      end do outerLoop

    end if

  end subroutine testFaceIntersection

end module CartesianCell_class