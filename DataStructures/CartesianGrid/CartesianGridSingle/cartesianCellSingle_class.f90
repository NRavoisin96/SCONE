module cartesianCellSingle_class
  
  use cartesianGenericProcedures, only : testIntervalIntersection
  use cartesianInitProcedures
  use edgeShelf_class,            only : edgeShelf
  use elementShelf_class,         only : elementShelf
  use face_inter,                 only : faceSATData
  use faceShelf_class,            only : faceShelf
  use genericProcedures,          only : append
  use numPrecision
  use universalVariables,         only : INSIDE_ELEMENT
  use vertexShelf_class,          only : vertexShelf

  implicit none
  private
  
  !!
  !!
  type, public :: cartesianCellSingle
    private
    integer(shortInt)               :: cellToEdgeIdx = 0, cellToElementIdx = 0, cellToVertexIdx = 0
    integer(shortInt), dimension(2) :: intersectedFaceIdxs = 0
  contains
    procedure :: cellTestEdgeIntersection
    procedure :: cellTestPolyhedronInclusion
    procedure :: cellTestPolyhedronInclusion_new
    procedure :: cellTestFaceIntersection
    procedure :: cellTestFaceIntersection_new
    procedure :: cellConstructMapSingleFace
    procedure :: setIsOutsideMesh
    procedure :: cellFinitePrecision
    procedure :: getCellIsOutside
    procedure :: getCellToEdgeIdx
    procedure :: getCellToElementIdx
    procedure :: getCellToVertexIdx
    procedure :: getFaceIdxs
    procedure :: kill
  end type cartesianCellSingle

contains

  !!
  !!
  !! (needs to be changed) (Inefficiency due to refactoring original code) (mapping matrices in cartesianCell_class)
  !! (can be moved and stored as a attribute in cartesianGrid. Then call testEdgeIntersection directly from cartesianGrid_class)
  !! (For now this is fine becase only 4% of initialisation time is increased by this - in face no increase in init time observed later.)
  !! (Tried this before and after specifying "elemental" and "pure". Keeping the original architecture (with separte class for cells))
  !! (is faster than without for initialisation, and negligible difference for in-cycle)
  pure subroutine cellTestEdgeIntersection(self, vertices, edges, edgeIdx, circumscribedBallRadius, targetDistance, centroid, &
                                           currEdgeVector, currVertexIdxs, a)
    class(cartesianCellSingle), intent(inout)   :: self
    class(vertexShelf), intent(in)              :: vertices
    class(edgeShelf), intent(in)                :: edges
    integer(shortInt), intent(in)               :: edgeIdx
    real(defReal), intent(in)                   :: circumscribedBallRadius, targetDistance, a
    real(defReal), dimension(3), intent(in)     :: centroid, currEdgeVector
    integer(shortInt), dimension(2), intent(in) :: currVertexIdxs

    call testEdgeIntersection(vertices, edges, edgeIdx, circumscribedBallRadius, targetDistance, centroid, &
                              currEdgeVector, currVertexIdxs, a, self % cellToVertexIdx, self % cellToEdgeIdx)

  end subroutine cellTestEdgeIntersection

  !!
  !!
  !!
  subroutine cellTestPolyhedronInclusion(self, faces, currElementFaceIdxs, centroid, &
                                         faceNormalSigns, elementIdx)
    class(cartesianCellSingle), intent(inout)           :: self
    class(faceShelf), intent(in)                        :: faces
    integer(shortInt), dimension(:), intent(in)         :: currElementFaceIdxs
    real(defReal), dimension(3), intent(in)             :: centroid
    real(defReal), dimension(:,:), intent(in)           :: faceNormalSigns
    integer(shortInt), intent(in)                       :: elementIdx

    call testPolyhedronInclusion(faces, currElementFaceIdxs, centroid, faceNormalSigns, elementIdx, self % cellToElementIdx)

  end subroutine cellTestPolyhedronInclusion

  !!
  !!
  !!
  subroutine cellTestPolyhedronInclusion_new(self, elementIdx, elementFaceIdxs, centroid, elements, faces)
    class(cartesianCellSingle), intent(inout)   :: self
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
    self % cellToElementIdx = elementIdx

  end subroutine cellTestPolyhedronInclusion_new

  !!
  !!
  !!
  subroutine cellTestFaceIntersection(self, vertices, edges, faces, currVertexIdxs, extraDistance, currFaceNormal, centroid, &
                                      cellSpacing, faceIdx, currFaceEdgeIdxs, targetDistance)
    class(cartesianCellSingle), intent(inout)   :: self
    class(vertexShelf), intent(in)              :: vertices
    class(edgeShelf), intent(in)                :: edges
    class(faceShelf), intent(in)                :: faces
    integer(shortInt), dimension(:), intent(in) :: currVertexIdxs, currFaceEdgeIdxs
    real(defReal), dimension(3), intent(in)     :: currFaceNormal, centroid
    real(defReal), intent(in)                   :: extraDistance, cellSpacing, targetDistance
    integer(shortInt), intent(in)               :: faceIdx

    call testFaceIntersection(vertices, edges, faces, currVertexIdxs, extraDistance, currFaceNormal, centroid, cellSpacing, &
                              faceIdx, currFaceEdgeIdxs, targetDistance, self % cellToElementIdx, self % cellToVertexIdx, &
                              self % cellToEdgeIdx, self % intersectedFaceIdxs)

  end subroutine cellTestFaceIntersection

  !!
  !!
  !!
  pure subroutine cellTestFaceIntersection_new(self, faceIdx, cellSpacing, extraDistance, targetDistance, centroid, edges, &
                                               cache, faces, vertices)
    class(cartesianCellSingle), intent(inout)    :: self
    integer(shortInt), intent(in)                :: faceIdx
    real(defReal), intent(in)                    :: cellSpacing, extraDistance, targetDistance
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
    if(0 < self % cellToEdgeIdx .and. 0 < self % cellToVertexIdx) return

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
      radii(1) = HALF * (abs(cache % edgeVectors(2, i)) + abs(cache % edgeVectors(3, i))) * cellSpacing
      radii(2) = HALF * (abs(cache % edgeVectors(1, i)) + abs(cache % edgeVectors(3, i))) * cellSpacing
      radii(3) = HALF * (abs(cache % edgeVectors(1, i)) + abs(cache % edgeVectors(2, i))) * cellSpacing

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
      self % cellToEdgeIdx = faceEdgeIdxs(1)

    else
      self % intersectedFaceIdxs(2) = faceIdx

      ! Check if the cell intersects the common edge between the two faces.
      commonEdgeIdx = faces % findCommonEdgeIdx(self % intersectedFaceIdxs(1), self % intersectedFaceIdxs(2))

      ! if there is no common edge, exit the subroutine early (other combinations of faces to be tried later)
      if (commonEdgeIdx == 0) return
      self % cellToEdgeIdx = 0

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
      self % cellToVertexIdx = merge(commonEdgeVertexIdxs(1), commonEdgeVertexIdxs(2), &
                                     distanceToVertex1Squared < distanceToVertex2Squared)
      if(targetDistance < distanceToVertex1Squared .and. targetDistance < distanceToVertex2Squared) &
      self % cellToEdgeIdx = commonEdgeIdx

    end if

  end subroutine cellTestFaceIntersection_new

  !!
  !!
  !!
  pure subroutine cellconstructMapSingleFace(self, faces)
    class(cartesianCellSingle), intent(inout)    :: self
    class(faceShelf), intent(in)                 :: faces
    integer(shortInt), dimension(:), allocatable :: currFaceEdgeIdxs!, currFaceVertexIdxs
    integer(shortInt)                            :: temp

    call constructMapSingleFace(faces, self % intersectedFaceIdxs, self % cellToEdgeIdx)

  end subroutine cellConstructMapSingleFace

  !!
  !!
  !!
  elemental subroutine setIsOutsideMesh(self)!, centroid)
    class(cartesianCellSingle), intent(inout)           :: self
    !real(defReal), dimension(3), intent(in)             :: centroid  !!!

    ! if a cell intersects with neither any edge nor face, and it is not contained in a single polyhedron,
    ! then, this cell lies outside the computational domain for the unstructured mesh
    if (all([self % cellToElementIdx, self % cellToVertexIdx, self % cellToEdgeIdx] == 0)) self % cellToElementIdx = -1

  end subroutine setIsOutsideMesh

  !!
  !!
  !!
  subroutine cellFinitePrecision(self, faces, elements, centroid, elementIdx)
    class(cartesianCellSingle), intent(inout)           :: self
    class(faceShelf), intent(in)                        :: faces
    class(elementShelf), intent(in)                     :: elements
    real(defReal), dimension(3), intent(in)             :: centroid
    integer(shortInt), intent(in)                       :: elementIdx

    ! If the current Cartesian cell does not intersect with any faces nor included in a single mesh element,
    ! Test if centroid lies inside any mesh element. If yes, finite precision error messed it up. Hence, 
    ! update chi mapping info to that mesh element. If not, this cell lies within another element or is outside mesh domain.
    if(all([self % cellToElementIdx, self % cellToVertexIdx, self % cellToEdgeIdx] == 0)) &
    call coverFinitePrecision(faces, elements, centroid, elementIdx, self % cellToElementIdx)

  end subroutine cellFinitePrecision

  !!
  !!
  !!
  elemental function getCellIsOutside(self) result(isOutside)
    class(cartesianCellSingle), intent(in) :: self
    logical(defBool)                       :: isOutside

    isOutside = self % cellToEdgeIdx == 0 .and. self % cellToElementIdx == 0 .and. self % cellToVertexIdx == 0

  end function getCellIsOutside

  !!
  !!
  !!
  elemental function getCellToEdgeIdx(self) result(cellToEdgeIdx)
    class(cartesianCellSingle), intent(in) :: self
    integer(shortInt)                      :: cellToEdgeIdx

    cellToEdgeIdx = self % cellToEdgeIdx

  end function getCellToEdgeIdx

  !!
  !!
  !!
  elemental function getCellToElementIdx(self) result(cellToElementIdx)
    class(cartesianCellSingle), intent(in) :: self
    integer(shortInt)                      :: cellToElementIdx

    cellToElementIdx = self % cellToElementIdx

  end function getCellToElementIdx

  !!
  !!
  !!
  elemental function getCellToVertexIdx(self) result(cellToVertexIdx)
    class(cartesianCellSingle), intent(in) :: self
    integer(shortInt)                      :: cellToVertexIdx

    cellToVertexIdx = self % cellToVertexIdx

  end function getCellToVertexIdx

  !!
  !!
  !!
  pure function getFaceIdxs(self) result(intersectedFaceIdxs)
    class(cartesianCellSingle), intent(in) :: self
    integer(shortInt), dimension(2)        :: intersectedFaceIdxs

    intersectedFaceIdxs = self % intersectedFaceIdxs

  end function

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(cartesianCellSingle), intent(inout) :: self

    self % cellToEdgeIdx = 0
    self % cellToElementIdx = 0
    self % cellToVertexIdx = 0
    self % intersectedFaceIdxs = 0

  end subroutine kill

end module CartesianCellSingle_class