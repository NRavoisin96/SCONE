module ASCGAcceleration_inter

  use accelerationStructure_inter, only : accelerationStructure, getStorageSize_super => getStorageSize
  use CartesianCell_class,         only : CartesianCell
  use CartesianGrid_class,         only : CartesianGrid
  use dictionary_class,            only : dictionary
  use edgeShelf_class,             only : edgeShelf
  use elementShelf_class,          only : elementShelf
  use errors_mod,                  only : fatalError
  use face_inter,                  only : faceSATData
  use faceShelf_class,             only : faceShelf
  use genericProcedures,           only : computePseudoAngle, crossProduct, findCommon, numToChar
  use numPrecision
  use patchSearchStatistics_mod
  use universalVariables,          only : INF
  use vertexShelf_class,           only : vertexShelf

  implicit none
  private

  ! Public procedures.
  public :: init

  !!
  !!
  !!
  type, abstract, public, extends(accelerationStructure) :: ASCGAcceleration
    private
    integer(shortInt)                              :: depth = 0, nSubGrids = 0
    logical(defBool)                               :: mapCells = .false., naiveInitialisation = .false., singleFaceShortcut = .true.
    real(defReal)                                  :: minimumAngle = ZERO, minimumEdgeLength = ZERO, targetDistance = ZERO, &
                                                      wStar = ZERO
    real(defReal), dimension(6)                    :: meshBounds = ZERO
    real(defReal), dimension(:), allocatable       :: spacings
    type(CartesianGrid)                            :: rootGrid
    type(CartesianGrid), dimension(:), allocatable :: subGrids
  contains
    procedure :: addSubGrid
    procedure :: computeCosineMaximumDihedralAngle
    procedure :: computeCosineMaximumFaceAngle
    procedure :: computeGeometricParameters
    procedure :: getDepth
    procedure :: getSingleFaceShortcut
    procedure :: getStorageSize
    procedure :: getWStar
    procedure :: init
    procedure :: kill
    procedure :: refineGrid
    procedure :: searchGrids
    procedure :: setMapCells
  end type ASCGAcceleration

contains
  !!
  !!
  !!
  elemental subroutine addSubGrid(self)
    class(ASCGAcceleration), intent(inout)         :: self
    integer(shortInt)                              :: nSubGrids
    type(CartesianGrid), dimension(:), allocatable :: temp

    ! Check if shelf has already been allocated and allocate it if so. Use 100 as an initial guess.
    if(.not. allocated(self % subGrids)) allocate(self % subGrids(100))
    self % nSubGrids = self % nSubGrids + 1

    ! Compute the number of subgrids in the shelf and expand it if needed.
    nSubGrids = size(self % subGrids)
    if(nSubGrids < self % nSubGrids) then
      allocate(temp(2 * nSubGrids))
      temp(1:nSubGrids) = self % subGrids
      call move_alloc(temp, self % subGrids)

    end if

  end subroutine addSubGrid

  !!
  !!
  !!
  pure function computeCosineMaximumDihedralAngle(self, edges, elements, faces) result(cosineMaximumDihedralAngle)
    class(ASCGAcceleration), intent(in)          :: self
    type(edgeShelf), intent(in)                  :: edges
    type(elementShelf), intent(in)               :: elements
    type(faceShelf), intent(in)                  :: faces
    integer(shortInt)                            :: i, j, k
    integer(shortInt), dimension(2)              :: commonFaceIdxs, candidateElementIdxs, signArray
    integer(shortInt), dimension(:), allocatable :: elementEdgeIdxs, elementFaceIdxs
    real(defReal)                                :: cosineMaximumDihedralAngle

    ! Initialise cosineMaximumDihedralAngle = -ONE then loop over all elements.
    cosineMaximumDihedralAngle = -ONE
    do i = 1, elements % getSize()
      elementFaceIdxs = abs(elements % getElementFaceIdxs(i))
      elementEdgeIdxs = elements % getElementEdgeIdxs(i)

      do j = 1, size(elementEdgeIdxs)
        ! Find indices of faces in the current element which share the current edge.
        commonFaceIdxs = findCommon(elementFaceIdxs, edges % getEdgeFaceIdxs(elementEdgeIdxs(j)))

        do k = 1, 2
          candidateElementIdxs = faces % getFaceElementIdxs(commonFaceIdxs(k))
          signArray(k) = merge(1, -1, any(i < candidateElementIdxs))

        end do
        ! Update cosineMaximumDihedralAngle.
        cosineMaximumDihedralAngle = max(cosineMaximumDihedralAngle, -dot_product(faces % getFaceNormal(commonFaceIdxs(1)), &
                                         faces % getFaceNormal(commonFaceIdxs(2))) * product(signArray))
        
      end do

    end do

  end function computeCosineMaximumDihedralAngle

  !!
  !!
  !!
  pure function computeCosineMaximumFaceAngle(self, edges, faces) result(cosineMaximumFaceAngle)
    class(ASCGAcceleration), intent(in)          :: self
    type(edgeShelf), intent(in)                  :: edges
    type(faceShelf), intent(in)                  :: faces
    integer(shortInt)                            :: i, j, k, l, nEdgesInFace, m
    integer(shortInt), dimension(2)              :: edge1VertexIdxs, edge2VertexIdxs
    integer(shortInt), dimension(:), allocatable :: faceEdgeIdxs
    real(defReal)                                :: cosineMaximumFaceAngle

    ! Initialise cosineMaximumFaceAngle = -ONE then loop over all faces.
    cosineMaximumFaceAngle = -ONE
    do i = 1, faces % getSize()
      faceEdgeIdxs = faces % getFaceEdgeIdxs(i)
      nEdgesInFace = size(faceEdgeIdxs)

      do j = 1, nEdgesInFace - 1
        edge1VertexIdxs = edges % getEdgeVertexIdxs(faceEdgeIdxs(j))

        do k = j + 1, nEdgesInFace
          edge2VertexIdxs = edges % getEdgeVertexIdxs(faceEdgeIdxs(k))

          do l = 1, 2
            do m = 1, 2
              ! Update maxCosValue. Correct the direction of unit vector of each edge.
              if (edge1VertexIdxs(l) == edge2VertexIdxs(m)) &
              cosineMaximumFaceAngle = max(cosineMaximumFaceAngle, &
                                           dot_product(edges % getEdgeUnitvector(faceEdgeIdxs(j)), &
                                                       edges % getEdgeUnitvector(faceEdgeIdxs(k))) * merge(ONE, -ONE, l == m))

            end do

          end do

        end do

      end do

    end do

  end function computeCosineMaximumFaceAngle

  !!
  !!
  !!
  pure subroutine computeGeometricParameters(self, edges, elements, faces, vertices)
    class(ASCGAcceleration), intent(inout) :: self
    type(edgeShelf), intent(in)            :: edges
    type(elementShelf), intent(in)         :: elements
    type(faceShelf), intent(in)            :: faces
    type(vertexShelf), intent(in)          :: vertices
    integer(shortInt)                      :: i

    ! Compute minimum edge length.
    self % minimumEdgeLength = INF
    do i = 1, edges % getSize()
      self % minimumEdgeLength = min(self % minimumEdgeLength, edges % getEdgeLength(i))

    end do

    ! Compute minimum angle.
    self % minimumAngle = acos(max(self % computeCosineMaximumDihedralAngle(edges, elements, faces), &
                                   self % computeCosineMaximumFaceAngle(edges, faces)))

  end subroutine computeGeometricParameters

  !!
  !!
  !!
  elemental function getDepth(self) result(depth)
    class(ASCGAcceleration), intent(in) :: self
    integer(shortInt)                   :: depth

    depth = self % depth

  end function getDepth

  !!
  !!
  !!
  elemental function getSingleFaceShortcut(self) result(singleFaceShortcut)
    class(ASCGAcceleration), intent(in) :: self
    logical(defBool)                    :: singleFaceShortcut

    singleFaceShortcut = self % singleFaceShortcut

  end function getSingleFaceShortcut

  !!
  !!
  !!
  elemental function getStorageSize(self) result(storageSize)
    class(ASCGAcceleration), intent(in) :: self
    integer(longInt)                    :: storageSize
    integer(shortInt)                   :: i

    storageSize = getStorageSize_super(self)
    if(allocated(self % spacings)) storageSize = storageSize + 8 * size(self % spacings)
    storageSize = storageSize + self % rootGrid % getStorageSize()
    if(allocated(self % subGrids)) then
      do i = 1, size(self % subGrids)
        storageSize = storageSize + self % subGrids(i) % getStorageSize()

      end do

    end if

  end function getStorageSize

  !!
  !!
  !!
  elemental function getWStar(self) result(wStar)
    class(ASCGAcceleration), intent(in) :: self
    real(defReal)                       :: wStar

    wStar = self % wStar

  end function getWStar

  !!
  !!
  !!
  subroutine init(self, dict, vertices, edges, faces, elements)
    class(ASCGAcceleration), intent(inout)         :: self
    type(dictionary), intent(in)                   :: dict
    type(vertexShelf), intent(in)                  :: vertices
    type(edgeShelf), intent(inout)                 :: edges
    type(faceShelf), intent(inout)                 :: faces
    type(elementShelf), intent(in)                 :: elements
    integer(shortInt)                              :: i, j, k, l
    integer(shortInt), dimension(3)                :: nCells
    real(defReal)                                  :: alphaGeneral, coarsestLayerSpacing, ratio, sineAlpha, sineHalfAlpha, &
                                                      cellVolume, spacing
    real(defReal), dimension(3)                    :: extraDistances
    real(defReal), dimension(6)                    :: gridBounds
    type(CartesianCell), pointer                   :: cellPtr
    type(CartesianGrid), dimension(:), allocatable :: temp
    character(*), parameter                        :: HERE = 'init (ASCGAcceleration_inter.f90)'

    ! Retrieve whether to use shortcut for single face intersections. Default to .true.
    call dict % getOrDefault(self % singleFaceShortcut, 'singleFaceShortcut', .true.)

    ! Retrieve whether to initialise structure without optimisations. Default to .false.
    call dict % getOrDefault(self % naiveInitialisation, 'naiveInitialisation', .false.)

    ! Retrieve number of layers from dictionary and allocate memory.
    call dict % getOrDefault(self % depth, 'depth', 1)
    if(self % depth < 1) call fatalError(HERE, 'Depth must be at least 1. Is: '//numToChar(self % depth)//'.')
    allocate(self % spacings(self % depth))

    ! Compute geometric parameters and angular sectors for edges.
    call self % computeGeometricParameters(edges, elements, faces, vertices)

    ! Set face constants and compute caches.
    if(.not. self % naiveInitialisation) then
      do i = 1, faces % getSize()
        call faces % computeFaceSATData(i, edges, vertices)

      end do

    end if

    ! Pre-compute sin(alpha) and sin(alpha / 2) then compute wStar and targetDistance.
    alphaGeneral = min(self % minimumAngle, THIRD * PI)
    sineAlpha = sin(alphaGeneral)
    sineHalfAlpha = sin(HALF * alphaGeneral)
    self % wStar = self % minimumEdgeLength * min(HALF, sin(self % minimumAngle))
    self % targetDistance = (self % wStar / (ONE + sineAlpha)) ** 2
    
    ! Now compute spacing in the finest layer and its inverse.
    self % spacings(self % depth) = TWO * self % wStar * (ONE - epsilon(ONE)) * sineAlpha * sineHalfAlpha / &
                                    (sqrt(THREE) * (ONE + sineAlpha) * (ONE + sineHalfAlpha))

    ! Set mesh bounds and compute bounds for the root grid.
    self % meshBounds = vertices % getExtremalCoordinates()
    if(1 < self % depth) then
      ! Estimate spacing in the coarsest layer. Use 100 subdivisions along the dimension of maximum extent in the mesh
      ! to avoid memory explosion.
      coarsestLayerSpacing = maxval(self % meshBounds(4:6) - self % meshBounds(1:3)) / 100

      ! Compute spacing ratio between the coarsest and finest layers. Use a minimum ratio of TWO.
      ratio = max(real(nint((coarsestLayerSpacing / self % spacings(self % depth)) ** (ONE / (self % depth - 1))), defReal), TWO)

      ! Now propagate spacings from the finest layer up to the coarsest layer.
      do i = self % depth - 1, 1, -1
        self % spacings(i) = self % spacings(i + 1) * ratio

      end do

    end if
    extraDistances = HALF * (self % spacings(1) - mod(self % meshBounds(4:6) - self % meshBounds(1:3), self % spacings(1)))
    gridBounds(1:3) = self % meshBounds(1:3) - extraDistances
    gridBounds(4:6) = self % meshBounds(4:6) + extraDistances

    ! Compute total volume here.
    totalVolume = (gridBounds(4) - gridBounds(1)) * (gridBounds(5) - gridBounds(2)) * (gridBounds(6) - gridBounds(3))

    ! Initialise root grid and map it.
    call self % rootGrid % init(self % spacings(1), gridBounds)
    call self % rootGrid % map([(i, i = 1, elements % getSize())], [(i, i = 1, faces % getSize())], self % depth == 1, &
                               self % mapCells, self % naiveInitialisation, self % targetDistance, edges, elements, &
                               faces, vertices)

    ! Refine root grid if needed.
    if(1 < self % depth) then
      call self % refineGrid(1, 0, edges, elements, faces, vertices)
      if(self % nSubGrids < size(self % subGrids)) then
        allocate(temp(self % nSubGrids))
        temp = self % subGrids(1:self % nSubGrids)
        call move_alloc(temp, self % subGrids)

      end if

    end if

    ! Now print volumes of each cell category.
    nCells = self % rootGrid % getCellsNumber()
    spacing = self % rootGrid % getSpacing()
    cellVolume = spacing * spacing * spacing
    do k = 1, nCells(3)
      do j = 1, nCells(2)
        do i = 1, nCells(1)
          cellPtr => self % rootGrid % getCellPtr(i, j, k)
          if(0 < cellPtr % getSubGridIdx()) cycle
          if(cellPtr % isOutside()) then
            outsideVolume = outsideVolume + cellVolume

          elseif(0 < cellPtr % getElementIdx()) then
            elementMappingVolume = elementMappingVolume + cellVolume

          elseif(cellPtr % intersectsOnlyOneFace() .and. 1 < self % depth) then
            singleFaceVolume = singleFaceVolume + cellVolume

          elseif(0 < cellPtr % getEdgeIdx()) then
            edgeMappingVolume = edgeMappingVolume + cellVolume

          else
            vertexMappingVolume = vertexMappingVolume + cellVolume

          end if

        end do

      end do

    end do

    if(allocated(self % subGrids)) then
      ! Loop over all sub-grids.
      do l = 1, size(self % subGrids)
        nCells = self % subGrids(l) % getCellsNumber()
        spacing = self % subGrids(l) % getSpacing()
        cellVolume = spacing * spacing * spacing

        do k = 1, nCells(3)
          do j = 1, nCells(2)
            do i = 1, nCells(1)
              cellPtr => self % subGrids(l) % getCellPtr(i, j, k)
              if(0 < cellPtr % getSubGridIdx()) cycle
              if(cellPtr % isOutside()) then
                outsideVolume = outsideVolume + cellVolume

              elseif(0 < cellPtr % getElementIdx()) then
                elementMappingVolume = elementMappingVolume + cellVolume

              elseif(cellPtr % intersectsOnlyOneFace() .and. 1 < self % depth) then
                singleFaceVolume = singleFaceVolume + cellVolume

              elseif(0 < cellPtr % getEdgeIdx()) then
                edgeMappingVolume = edgeMappingVolume + cellVolume

              else
                vertexMappingVolume = vertexMappingVolume + cellVolume

              end if

            end do

          end do

        end do

      end do

    end if

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(ASCGAcceleration), intent(inout) :: self
    integer(shortInt)                      :: i

    ! Local.
    self % depth = 0
    self % nSubGrids = 0
    self % mapCells = .false.
    self % naiveInitialisation = .false.
    self % singleFaceShortcut = .true.
    self % minimumAngle = ZERO
    self % minimumEdgeLength = ZERO
    self % targetDistance = ZERO
    self % wStar = ZERO
    if(allocated(self % spacings)) deallocate(self % spacings)
    call self % rootGrid % kill()
    if(allocated(self % subGrids)) then
      do i = 1, size(self % subGrids)
        call self % subGrids(i) % kill()

      end do
      deallocate(self % subGrids)

    end if

  end subroutine kill

  !!
  !!
  !!
  recursive subroutine refineGrid(self, depth, subGridIdx, edges, elements, faces, vertices)
    class(ASCGAcceleration), target, intent(inout) :: self
    integer(shortInt), intent(in)                  :: depth, subGridIdx
    type(edgeShelf), intent(in)                    :: edges
    type(elementShelf), intent(in)                 :: elements
    type(faceShelf), intent(inout)                 :: faces
    type(vertexShelf), intent(in)                  :: vertices
    integer(shortInt)                              :: i, j, k, newDepth
    integer(shortInt), dimension(3)                :: nCells
    integer(shortInt), dimension(:), allocatable   :: cellIntersectedFaceIdxs
    real(defReal)                                  :: gridSpacing
    real(defReal), dimension(6)                    :: cellBounds, gridBounds
    type(CartesianGrid), pointer                   :: currentGridPtr

    if(self % depth <= depth) return

    ! Get pointer to appropriate grid.
    if(subGridIdx == 0) then
      currentGridPtr => self % rootGrid

    else
      currentGridPtr => self % subGrids(subGridIdx)

    end if

    ! Loop over all cells in the current grid.
    nCells = currentGridPtr % getCellsNumber()
    gridSpacing = currentGridPtr % getSpacing()
    gridBounds = currentGridPtr % getBounds()
    do k = 1, nCells(3)
      cellBounds(6) = gridBounds(3) + gridSpacing * k
      cellBounds(3) = cellBounds(6) - gridSpacing
      do j = 1, nCells(2)
        cellBounds(5) = gridBounds(2) + gridSpacing * j
        cellBounds(2) = cellBounds(5) - gridSpacing
        do i = 1, nCells(1)
          ! Check if the current cell needs to be refined and request index for next available subgrid.
          if(currentGridPtr % isCellSimple(i, j, k, self % singleFaceShortcut)) cycle
          call self % addSubGrid()

          ! Re-acquire pointer since memory reallocation may have corrupted pointers.
          if(subGridIdx == 0) then
            currentGridPtr => self % rootGrid

          else
            currentGridPtr => self % subGrids(subGridIdx)

          end if

          ! Set the current grid cell's subgrid index.
          call currentGridPtr % setCellSubGridIdx(i, j, k, self % nSubGrids)

          ! Initialise current sub grid.
          newDepth = depth + 1
          cellBounds(4) = gridBounds(1) + gridSpacing * i
          cellBounds(1) = cellBounds(4) - gridSpacing
          call self % subGrids(self % nSubGrids) % init(self % spacings(newDepth), cellBounds)

          ! Get indices of faces intersected by the current cell, map current subgrid then refine if deepest layer has not been 
          ! reached yet.
          cellIntersectedFaceIdxs = currentGridPtr % getCellIntersectedFaceIdxs(i, j, k)
          call self % subGrids(self % nSubGrids) % map(faces % getFaceElementIdxs(cellIntersectedFaceIdxs), &
                                                       cellIntersectedFaceIdxs, self % depth == newDepth, self % mapCells, &
                                                       self % naiveInitialisation, self % targetDistance, edges, elements, &
                                                       faces, vertices)
          call self % refineGrid(newDepth, self % nSubGrids, edges, elements, faces, vertices)

          ! Now re-acquire pointer.
          if(subGridIdx == 0) then
            currentGridPtr => self % rootGrid

          else
            currentGridPtr => self % subGrids(subGridIdx)

          end if

        end do

      end do

    end do

  end subroutine refineGrid

  !!
  !!
  !!
  function searchGrids(self, r) result(deepestCellPtr)
    class(ASCGAcceleration), target, intent(in) :: self
    real(defReal), dimension(3), intent(in)     :: r
    integer(shortInt)                           :: subGridIdx
    integer(shortInt), dimension(3)             :: hostCellIdxs
    type(CartesianCell), pointer                :: deepestCellPtr
    type(CartesianGrid), pointer                :: currentGridPtr

    ! Start at the root grid then search until we hit a terminal cell.
    currentGridPtr => self % rootGrid
    do
      hostCellIdxs = currentGridPtr % findHostCellIdxs(r)
      subGridIdx = currentGridPtr % getCellSubGridIdx(hostCellIdxs(1), hostCellIdxs(2), hostCellIdxs(3))
      if(subGridIdx == 0) exit
      currentGridPtr => self % subGrids(subGridIdx)
      
    end do
    deepestCellPtr => currentGridPtr % getCellPtr(hostCellIdxs(1), hostCellIdxs(2), hostCellIdxs(3))

  end function searchGrids

  !!
  !!
  !!
  elemental subroutine setMapCells(self, mapCells)
    class(ASCGAcceleration), intent(inout) :: self
    logical(defBool), intent(in)           :: mapCells

    self % mapCells = mapCells

  end subroutine setMapCells

end module ASCGAcceleration_inter