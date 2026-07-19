module CartesianCell_class

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
    procedure :: addIntersectedFaceIdx
    procedure :: getEdgeIdx
    procedure :: getElementIdx
    procedure :: getFirstIntersectedFaceIdx
    procedure :: getIntersectedFaceIdxs
    procedure :: getIntersectedFacesNumber
    procedure :: getStorageSize
    procedure :: getSubGridIdx
    procedure :: getVertexIdx
    procedure :: intersectsOnlyOneFace
    procedure :: isOutside
    procedure :: isSimple
    procedure :: isUnprocessed
    procedure :: kill
    procedure :: map
    procedure :: setElementIdx
    procedure :: setSubGridIdx
  end type CartesianCell

contains
  !!
  !!
  !!
  elemental subroutine addIntersectedFaceIdx(self, intersectedFaceIdx)
    class(CartesianCell), intent(inout) :: self
    integer(shortInt), intent(in)       :: intersectedFaceIdx

    call append(self % intersectedFaceIdxs, intersectedFaceIdx)

  end subroutine addIntersectedFaceIdx

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
  elemental function getFirstIntersectedFaceIdx(self) result(firstIntersectedFaceIdx)
    class(CartesianCell), intent(in) :: self
    integer(shortInt)                :: firstIntersectedFaceIdx

    firstIntersectedFaceIdx = self % intersectedFaceIdxs(1)

  end function getFirstIntersectedFaceIdx

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
  elemental function getIntersectedFacesNumber(self) result(nIntersectedFaces)
    class(CartesianCell), intent(in) :: self
    integer(shortInt)                :: nIntersectedFaces

    if(allocated(self % intersectedFaceIdxs)) then
      nIntersectedFaces = size(self % intersectedFaceIdxs)

    else
      nIntersectedFaces = 0

    end if

  end function getIntersectedFacesNumber

  !!
  !!
  !!
  elemental function getStorageSize(self) result(storageSize)
    class(CartesianCell), intent(in) :: self
    integer(longInt)                 :: storageSize

    storageSize = storage_size(self) / 8
    if(allocated(self % intersectedFaceIdxs)) storageSize = storageSize + 4 * size(self % intersectedFaceIdxs)

  end function getStorageSize

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
  elemental function intersectsOnlyOneFace(self) result(doesIt)
    class(CartesianCell), intent(in) :: self
    logical(defBool)                 :: doesIt

    doesIt = .false.
    if(allocated(self % intersectedFaceIdxs)) then
      doesIt = size(self % intersectedFaceIdxs) == 1

    end if

  end function intersectsOnlyOneFace

  !!
  !!
  !!
  elemental function isOutside(self) result(isIt)
    class(CartesianCell), intent(in) :: self
    logical(defBool)                 :: isIt

    isIt = self % elementIdx == 0 .and. .not. allocated(self % intersectedFaceIdxs)

  end function isOutside

  !!
  !!
  !!
  elemental function isSimple(self, singleFaceShortcut) result(isIt)
    class(CartesianCell), intent(in) :: self
    logical(defBool), intent(in)     :: singleFaceShortcut
    logical(defBool)                 :: isIt

    isIt = (0 < self % elementIdx) .or. (self % isOutside()) .or. &
           (self % getIntersectedFacesNumber() < 2 .and. singleFaceShortcut)

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
  pure subroutine map(self, targetDistance, centroid, edges, faces, vertices)
    class(CartesianCell), intent(inout)          :: self
    real(defReal), intent(in)                    :: targetDistance
    real(defReal), dimension(3), intent(in)      :: centroid
    type(edgeShelf), intent(in)                  :: edges
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    integer(shortInt)                            :: commonEdgeIdx, i, j, nIntersectedFaces
    integer(shortInt), dimension(2)              :: commonEdgeVertexIdxs
    integer(shortInt), dimension(:), allocatable :: faceEdgeIdxs
    real(defReal)                                :: distanceToVertex1Squared, distanceToVertex2Squared, projection
    real(defReal), dimension(3)                  :: commonEdgeUnitVector, commonEdgeVertex1Coords, temp

    ! Return immediately if mapping has already been assigned for this cell.
    if(.not. allocated(self % intersectedFaceIdxs) .or. (0 < self % edgeIdx .and. 0 < self % vertexIdx)) return
    nIntersectedFaces = size(self % intersectedFaceIdxs)
    if(nIntersectedFaces == 1) then
      ! Assign edge mapping to first edge in the face.
      faceEdgeIdxs = faces % getFaceEdgeIdxs(self % intersectedFaceIdxs(1))
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

  end subroutine map

  !!
  !!
  !!
  elemental subroutine setElementIdx(self, elementIdx)
    class(CartesianCell), intent(inout) :: self
    integer(shortInt), intent(in)       :: elementIdx

    self % elementIdx = elementIdx

  end subroutine setElementIdx

  !!
  !!
  !!
  elemental subroutine setSubGridIdx(self, subGridIdx)
    class(CartesianCell), intent(inout) :: self
    integer(shortInt), intent(in)       :: subGridIdx

    self % subGridIdx = subGridIdx

  end subroutine setSubGridIdx

end module CartesianCell_class