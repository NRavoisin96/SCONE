module CartesianGrid_class

  use CartesianCell_class,        only : CartesianCell
  use cartesianGenericProcedures, only : calculateAvgEdgeLength, constructAABB, findMinDihedralAngle, findMinFaceAngle
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
    procedure :: findHostCellIdxs
    procedure :: getBounds
    procedure :: getCellsNumber
    procedure :: getCellSubGridIdx
    procedure :: getCellEdgeIdx
    procedure :: getCellElementIdx
    procedure :: getCellIntersectedFaceIdxs
    procedure :: getCellVertexIdx
    procedure :: getInverseSpacing
    procedure :: getSpacing
    procedure :: init
    procedure :: isCellOutside
    procedure :: isCellSimple
    procedure :: isCellUnprocessed
    procedure :: isOutsideBounds
    procedure :: kill
    procedure :: setBounds
    procedure :: setCellsNumber
    procedure :: setCellSubGridIdx
    procedure :: testCellElementInclusion
    procedure :: testCellFaceIntersection
  end type CartesianGrid

contains
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
  elemental function isCellSimple(self, xIdx, yIdx, zIdx) result(isIt)
    class(CartesianGrid), intent(in) :: self
    integer(shortInt), intent(in)    :: xIdx, yIdx, zIdx
    logical(defBool)                 :: isIt

    isIt = self % cells(xIdx, yIdx, zIdx) % isSimple()

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

  !!
  !!
  !!
  pure subroutine testCellElementInclusion(self, elementIdx, xIdx, yIdx, zIdx, elementFaceIdxs, centroid, elements, faces)
    class(CartesianGrid), intent(inout)         :: self
    integer(shortInt), intent(in)               :: elementIdx, xIdx, yIdx, zIdx
    integer(shortInt), dimension(:), intent(in) :: elementFaceIdxs
    real(defReal), dimension(3), intent(in)     :: centroid
    type(elementShelf), intent(in)              :: elements
    type(faceShelf), intent(in)                 :: faces

    call self % cells(xIdx, yIdx, zIdx) % testElementInclusion(elementIdx, elementFaceIdxs, centroid, elements, faces)

  end subroutine testCellElementInclusion

  !!
  !!
  !!
  subroutine testCellFaceIntersection(self, faceIdx, xIdx, yIdx, zIdx, extraDistance, targetDistance, centroid, edges, &
                                           cache, faces, vertices)
    class(CartesianGrid), intent(inout)     :: self
    integer(shortInt), intent(in)           :: faceIdx, xIdx, yIdx, zIdx
    real(defReal), intent(in)               :: extraDistance, targetDistance
    real(defReal), dimension(3), intent(in) :: centroid
    type(edgeShelf), intent(in)             :: edges
    type(faceSATData), intent(in)           :: cache
    type(faceShelf), intent(in)             :: faces
    type(vertexShelf), intent(in)           :: vertices

    call self % cells(xIdx, yIdx, zIdx) % testFaceIntersection(faceIdx, extraDistance, self % spacing, targetDistance, &
                                                               centroid, edges, cache, faces, vertices)

  end subroutine testCellFaceIntersection

end module CartesianGrid_class