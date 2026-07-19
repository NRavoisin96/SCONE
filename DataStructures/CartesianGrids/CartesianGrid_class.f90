module CartesianGrid_class

  use CartesianCell_class,        only : CartesianCell
  use dictionary_class,           only : dictionary
  use edgeShelf_class,            only : edgeShelf
  use elementShelf_class,         only : elementShelf
  use face_inter,                 only : faceSATData
  use faceShelf_class,            only : faceShelf
  use genericProcedures,          only : crossProduct, findCommon
  use numPrecision
  use universalVariables,         only : INF
  use vertexShelf_class,          only : vertexShelf

  implicit none
  private

  !!
  !!
  !!
  type, public :: CartesianGrid
    private
    integer(shortInt), dimension(3)                      :: nCells = 0
    real(defReal)                                        :: spacing = ZERO, inverseSpacing = ZERO
    real(defReal), dimension(6)                          :: bounds = ZERO
    type(CartesianCell), dimension(:, :, :), allocatable :: cells
  contains
    procedure :: constructCellIdxs
    procedure :: findHostCellIdxs
    procedure :: getBounds
    procedure :: getCellEdgeIdx
    procedure :: getCellElementIdx
    procedure :: getCellIntersectedFaceIdxs
    procedure :: getCellPtr
    procedure :: getCellsNumber
    procedure :: getCellSubGridIdx
    procedure :: getCellVertexIdx
    procedure :: getInverseSpacing
    procedure :: getSpacing
    procedure :: getStorageSize
    procedure :: init
    procedure :: isCellOutside
    procedure :: isCellSimple
    procedure :: isCellUnprocessed
    procedure :: isOutsideBounds
    procedure :: kill
    procedure :: map
    procedure :: mapCell
    procedure :: setBounds
    procedure :: setCellsNumber
    procedure :: setCellSubGridIdx
  end type CartesianGrid

contains
  !!
  !!
  !!
  pure function constructCellIdxs(self, vertexIdxs, vertices) result(cellIdxs)
    class(CartesianGrid), intent(in)            :: self
    integer(shortInt), dimension(:), intent(in) :: vertexIdxs
    type(vertexShelf), intent(in)               :: vertices
    integer(shortInt)                           :: i
    integer(shortInt), dimension(6)             :: cellIdxs
    real(defreal), dimension(3)                 :: xyz_min, xyz_max, vertexCoords

    ! initialise xyz_min and xyz_max using the first vertex
    vertexCoords = vertices % getVertexCoordinates(vertexIdxs(1))
    xyz_max = vertexCoords
    xyz_min = vertexCoords

    ! find xyz_min and xyz_max 
    do i = 2, size(vertexIdxs)
      vertexCoords = vertices % getVertexCoordinates(vertexIdxs(i))
      xyz_min = min(xyz_min, vertexCoords)
      xyz_max = max(xyz_max, vertexCoords)

    end do

    ! Find cell indices.
    cellIdxs(1:3) = ceiling((xyz_min - self % bounds(1:3)) * self % inverseSpacing)
    cellIdxs(4:6) = ceiling((xyz_max - self % bounds(1:3)) * self % inverseSpacing)

  end function constructCellIdxs

  !!
  !!
  !!
  pure function findHostCellIdxs(self, r) result(cellIdxs)
    class(CartesianGrid), intent(in)        :: self
    real(defReal), dimension(3), intent(in) :: r
    integer(shortInt), dimension(3)         :: cellIdxs

    cellIdxs = ceiling((r - self % bounds(1:3)) * self % inverseSpacing)

  end function findHostCellIdxs

  !!
  !!
  !!
  pure function getBounds(self) result(bounds)
    class(CartesianGrid), intent(in) :: self
    real(defReal), dimension(6)      :: bounds

    bounds = self % bounds

  end function getBounds

!!
  !!
  !!
  elemental function getCellEdgeIdx(self, xIdx, yIdx, zIdx) result(edgeIdx)
    class(CartesianGrid), intent(in) :: self
    integer(shortInt), intent(in)    :: xIdx, yIdx, zIdx
    integer(shortInt)                :: edgeIdx

    edgeIdx = self % cells(xIdx, yIdx, zIdx) % getEdgeIdx()

  end function getCellEdgeIdx

  !!
  !!
  !!
  elemental function getCellElementIdx(self, xIdx, yIdx, zIdx) result(elementIdx)
    class(CartesianGrid), intent(in) :: self
    integer(shortInt), intent(in)    :: xIdx, yIdx, zIdx
    integer(shortInt)                :: elementIdx

    elementIdx = self % cells(xIdx, yIdx, zIdx) % getElementIdx()

  end function getCellElementIdx

  !!
  !!
  !!
  pure function getCellIntersectedFaceIdxs(self, xIdx, yIdx, zIdx) result(intersectedFaceIdxs)
    class(CartesianGrid), intent(in)             :: self
    integer(shortInt), intent(in)                :: xIdx, yIdx, zIdx
    integer(shortInt), dimension(:), allocatable :: intersectedFaceIdxs

    intersectedFaceIdxs = self % cells(xIdx, yIdx, zIdx) % getIntersectedFaceIdxs()

  end function getCellIntersectedFaceIdxs

  !!
  !!
  !!
  function getCellPtr(self, xIdx, yIdx, zIdx) result(cellPtr)
    class(CartesianGrid), target, intent(in) :: self
    integer(shortInt), intent(in)            :: xIdx, yIdx, zIdx
    type(CartesianCell), pointer             :: cellPtr

    cellPtr => self % cells(xIdx, yIdx, zIdx)

  end function getCellPtr

  !!
  !!
  !!
  pure function getCellsNumber(self) result(nCells)
    class(CartesianGrid), intent(in) :: self
    integer(shortInt), dimension(3)  :: nCells

    nCells = self % nCells

  end function getCellsNumber

  !!
  !!
  !!
  elemental function getCellSubGridIdx(self, xIdx, yIdx, zIdx) result(subGridIdx)
    class(CartesianGrid), intent(in) :: self
    integer(shortInt), intent(in)    :: xIdx, yIdx, zIdx
    integer(shortInt)                :: subGridIdx

    subGridIdx = self % cells(xIdx, yIdx, zIdx) % getSubGridIdx()

  end function getCellSubGridIdx

  !!
  !!
  !!
  elemental function getCellVertexIdx(self, xIdx, yIdx, zIdx) result(vertexIdx)
    class(CartesianGrid), intent(in) :: self
    integer(shortInt), intent(in)    :: xIdx, yIdx, zIdx
    integer(shortInt)                :: vertexIdx

    vertexIdx = self % cells(xIdx, yIdx, zIdx) % getVertexIdx()

  end function getCellVertexIdx

  !!
  !!
  !!
  elemental function getInverseSpacing(self) result(inverseSpacing)
    class(CartesianGrid), intent(in) :: self
    real(defReal)                    :: inverseSpacing

    inverseSpacing = self % inverseSpacing

  end function getInverseSpacing

  !!
  !!
  !!
  elemental function getSpacing(self) result(spacing)
    class(CartesianGrid), intent(in) :: self
    real(defReal)                    :: spacing

    spacing = self % spacing

  end function getSpacing

  !!
  !!
  !!
  elemental function getStorageSize(self) result(storageSize)
    class(CartesianGrid), intent(in) :: self
    integer(longInt)                 :: storageSize
    integer(shortInt)                :: i, j, k

    storageSize = storage_size(self) / 8
    if(allocated(self % cells)) then
      do k = 1, size(self % cells, 3)
        do j = 1, size(self % cells, 2)
          do i = 1, size(self % cells, 1)
            storageSize = storageSize + self % cells(i, j, k) % getStorageSize()

          end do

        end do

      end do

    end if

  end function getStorageSize

  !!
  !!
  !!
  subroutine init(self, spacing, bounds)
    class(CartesianGrid), intent(inout)     :: self
    real(defReal), intent(in)               :: spacing
    real(defReal), dimension(6), intent(in) :: bounds

    ! Set spacing and bounds then compute number of cells in each dimension.
    self % spacing = spacing
    self % inverseSpacing = ONE / self % spacing
    self % bounds = bounds
    self % nCells = nint((self % bounds(4:6) - self % bounds(1:3)) * self % inverseSpacing)

    ! Allocate memory.
    allocate(self % cells(self % nCells(1), self % nCells(2), self % nCells(3)))

  end subroutine init

  !!
  !!
  !!
  pure function isCellOutside(self, xIdx, yIdx, zIdx) result(isIt)
    class(CartesianGrid), intent(in) :: self
    integer(shortInt), intent(in)    :: xIdx, yIdx, zIdx
    logical(defBool)                 :: isIt

    isIt = self % cells(xIdx, yIdx, zIdx) % isOutside()

  end function isCellOutside

  !!
  !!
  !!
  elemental function isCellSimple(self, xIdx, yIdx, zIdx, singleFaceShortcut) result(isIt)
    class(CartesianGrid), intent(in) :: self
    integer(shortInt), intent(in)    :: xIdx, yIdx, zIdx
    logical(defBool), intent(in)     :: singleFaceShortcut
    logical(defBool)                 :: isIt

    isIt = self % cells(xIdx, yIdx, zIdx) % isSimple(singleFaceShortcut)

  end function isCellSimple

  !!
  !!
  !!
  elemental function isCellUnprocessed(self, xIdx, yIdx, zIdx) result(isIt)
    class(CartesianGrid), intent(in) :: self
    integer(shortInt), intent(in)    :: xIdx, yIdx, zIdx
    logical(defBool)                 :: isIt

    isIt = self % cells(xIdx, yIdx, zIdx) % isUnprocessed()

  end function isCellUnprocessed

  !!
  !!
  !!
  pure function isOutsideBounds(self, r) result(isIt)
    class(CartesianGrid), intent(in)        :: self
    real(defReal), dimension(3), intent(in) :: r
    logical(defBool)                        :: isIt

    isIt = any(r < self % bounds(1:3)) .or. any(self % bounds(4:6) < r)

  end function isOutsideBounds

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(CartesianGrid), intent(inout) :: self
    integer(shortInt)                   :: i, j, k

    ! Local.
    self % nCells = 0
    self % spacing = ZERO
    self % inverseSpacing = ZERO
    self % bounds = ZERO

    if(allocated(self % cells)) then
      do k = 1, size(self % cells, 3)
        do j = 1, size(self % cells, 2)
          do i = 1, size(self % cells, 1)
            call self % cells(i, j, k) % kill()

          end do

        end do

      end do
      deallocate(self % cells)

    end if

  end subroutine kill

  !!
  !!
  !!
  subroutine map(self, elementIdxs, faceIdxs, isFinestLayer, mapCells, naiveInitialisation, targetDistance, edges, &
                 elements, faces, vertices)
    class(CartesianGrid), intent(inout)          :: self
    integer(shortInt), dimension(:), intent(in)  :: elementIdxs, faceIdxs
    logical(defBool), intent(in)                 :: isFinestLayer, mapCells, naiveInitialisation
    real(defReal), intent(in)                    :: targetDistance
    type(edgeShelf), intent(in)                  :: edges
    type(elementShelf), intent(in)               :: elements
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    integer(shortInt)                            :: c, i, j, k, l
    integer(shortInt), dimension(:), allocatable :: cellIdxs, elementFaceIdxs, elementVertexIdxs, faceVertexIdxs
    logical(defBool)                             :: isInside
    real(defReal), dimension(3)                  :: centroid, corner
    real(defReal), dimension(3, 8), parameter    :: CORNERS = reshape([-ONE,-ONE,-ONE, ONE,-ONE,-ONE, -ONE, ONE,-ONE, &
                                                                       ONE, ONE,-ONE, -ONE,-ONE, ONE, ONE,-ONE, ONE, &
                                                                       -ONE, ONE, ONE, ONE, ONE, ONE], [3, 8])

    ! Loop over all faces.
    do i = 1, size(faceIdxs)
      ! Generate the indices of the cells contained in the current face's AABB.
      faceVertexIdxs = faces % getFaceVertexIdxs(faceIdxs(i))
      cellIdxs = self % constructCellIdxs(faceVertexIdxs, vertices)

      ! Loop over all cells.
      do l = max(1, cellIdxs(3)), min(self % nCells(3), cellIdxs(6))
        centroid(3) = self % bounds(3) + self % spacing * (l - HALF)
        do k = max(1, cellIdxs(2)), min(self % nCells(2), cellIdxs(5))
          centroid(2) = self % bounds(2) + self % spacing * (k - HALF)
          do j = max(1, cellIdxs(1)), min(self % nCells(1), cellIdxs(4))
            centroid(1) = self % bounds(1) + self % spacing * (j - HALF)
            ! If initialising without any optimisations, re-compute the SAT cache of the face before testing for an intersection.
            if(naiveInitialisation) then
              ! Test current cell for intersection with the current face.
              if(faces % intersectsFace_naive(faceIdxs(i), self % spacing, centroid, edges, vertices)) &
                call self % cells(j, k, l) % addIntersectedFaceIdx(faceIdxs(i)) 

            else
              ! Test current cell for intersection with the current face.
              if(faces % intersectsFace(faceIdxs(i), self % spacing, centroid)) &
                call self % cells(j, k, l) % addIntersectedFaceIdx(faceIdxs(i)) 

            end if

            ! If we are at the finest layer, map the cell.
            if(mapCells .and. isFinestLayer) call self % cells(j, k, l) % map(targetDistance, centroid, edges, faces, vertices)

          end do

        end do

      end do

    end do

    ! Now loop over all elements.
    do i = 1, size(elementIdxs)
      ! Generate the indices of the cells contained in the current element's AABB.
      elementFaceIdxs = elements % getElementFaceIdxs(elementIdxs(i))
      elementVertexIdxs = elements % getElementVertexIdxs(elementIdxs(i))
      cellIdxs = self % constructCellIdxs(elementVertexIdxs, vertices)

      ! Loop over all cells.
      do l = max(1, cellIdxs(3)), min(self % nCells(3), cellIdxs(6))
        centroid(3) = self % bounds(3) + self % spacing * (l - HALF)
        do k = max(1, cellIdxs(2)), min(self % nCells(2), cellIdxs(5))
          centroid(2) = self % bounds(2) + self % spacing * (k - HALF)
          do j = max(1, cellIdxs(1)), min(self % nCells(1), cellIdxs(4))
            if(naiveInitialisation) then
              ! Compute centroid unconditionally and initialise isInside = .true.
              centroid(1) = self % bounds(1) + self % spacing * (j - HALF)
              isInside = .true.

              ! Check if each corner is inside and abort check if not.
              do c = 1, 8
                corner = centroid + HALF * self % spacing * CORNERS(:, c)
                if(.not. elements % isPointInsideElementNoBoundaryCheck(elementIdxs(i), corner, faces)) then
                  isInside = .false.
                  exit

                end if

              end do
              ! If still inside then the cell is contained within the element.
              if(isInside) call self % cells(j, k, l) % setElementIdx(elementIdxs(i))

            elseif(self % cells(j, k, l) % isUnprocessed()) then
              centroid(1) = self % bounds(1) + self % spacing * (j - HALF)
              if(elements % isPointInsideElementNoBoundaryCheck(elementIdxs(i), centroid, faces)) &
              call self % cells(j, k, l) % setElementIdx(elementIdxs(i))

            end if

          end do

        end do

      end do

    end do

  end subroutine map

  !!
  !!
  !!
  pure subroutine mapCell(self, xIdx, yIdx, zIdx, targetDistance, centroid, edges, faces, vertices)
    class(CartesianGrid), intent(inout)     :: self
    integer(shortInt), intent(in)           :: xIdx, yIdx, zIdx
    real(defReal), intent(in)               :: targetDistance
    real(defReal), dimension(3), intent(in) :: centroid
    type(edgeShelf), intent(in)             :: edges
    type(faceShelf), intent(in)             :: faces
    type(vertexShelf), intent(in)           :: vertices

    call self % cells(xIdx, yIdx, zIdx) % map(targetDistance, centroid, edges, faces, vertices)

  end subroutine mapCell

  !!
  !!
  !!
  pure subroutine setCellsNumber(self, nCells)
    class(CartesianGrid), intent(inout)         :: self
    integer(shortInt), dimension(3), intent(in) :: nCells

    self % nCells = nCells

  end subroutine setCellsNumber

  !!
  !!
  !!
  elemental subroutine setCellSubGridIdx(self, xIdx, yIdx, zIdx, subGridIdx)
    class(CartesianGrid), intent(inout) :: self
    integer(shortInt), intent(in)       :: xIdx, yIdx, zIdx, subGridIdx

    call self % cells(xIdx, yIdx, zIdx) % setSubGridIdx(subGridIdx)

  end subroutine setCellSubGridIdx

  !!
  !!
  !!
  pure subroutine setBounds(self, bounds)
    class(CartesianGrid), intent(inout)     :: self
    real(defReal), dimension(6), intent(in) :: bounds

    self % bounds = bounds

  end subroutine setBounds

end module CartesianGrid_class