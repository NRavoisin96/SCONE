module patchSearchAcceleration_class

  use accelerationStructure_inter, only : accelerationStructure
  use CartesianGrid_class,         only : CartesianGrid
  use coord_class,                 only : coord
  use dictionary_class,            only : dictionary
  use edgeShelf_class,             only : edgeShelf
  use elementShelf_class,          only : elementShelf
  use errors_mod,                  only : fatalError
  use face_inter,                  only : faceSATData
  use faceShelf_class,             only : faceShelf
  use genericProcedures,           only : computePseudoAngle, crossProduct, findCommon, numToChar
  use numPrecision
  use universalVariables,          only : INF
  use vertexShelf_class,           only : vertexShelf

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(accelerationStructure) :: patchSearchAcceleration
    private
    integer(shortInt)                              :: nLayers = 0, nSubGrids = 0
    real(defReal)                                  :: minimumAngle = ZERO, minimumEdgeLength = ZERO, targetDistance = ZERO, &
                                                      wStar = ZERO
    real(defReal), dimension(6)                    :: meshBounds
    real(defReal), dimension(:), allocatable       :: spacings
    type(CartesianGrid)                            :: rootGrid
    type(CartesianGrid), dimension(:), allocatable :: subGrids
    type(faceSATData), dimension(:), allocatable   :: faceCaches
  contains
    procedure :: addSubGrid
    procedure :: computeCosineMaximumDihedralAngle
    procedure :: computeCosineMaximumFaceAngle
    procedure :: computeEdgeAngularSectors
    procedure :: computeGeometricParameters
    procedure :: findHostElement
    procedure :: init
    procedure :: kill
    procedure :: mapGrid
    procedure :: refineGrid
    procedure :: searchGrids
  end type patchSearchAcceleration

contains
  !!
  !!
  !!
  elemental subroutine addSubGrid(self)
    class(patchSearchAcceleration), intent(inout)  :: self
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
    class(patchSearchAcceleration), intent(in)   :: self
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
    class(patchSearchAcceleration), intent(in)   :: self
    type(edgeShelf), intent(in)                  :: edges
    type(faceShelf), intent(in)                  :: faces
    integer(shortInt)                            :: i, j, k, l , m
    integer(shortInt), dimension(2)              :: edge1VertexIdxs, edge2VertexIdxs
    integer(shortInt), dimension(:), allocatable :: faceEdgeIdxs
    real(defReal)                                :: cosineMaximumFaceAngle

    ! Initialise cosineMaximumFaceAngle = -ONE then loop over all faces.
    cosineMaximumFaceAngle = -ONE
    do i = 1, faces % getSize()
      faceEdgeIdxs = faces % getFaceEdgeIdxs(i)

      do j = 1, size(faceEdgeIdxs) - 1
        edge1VertexIdxs = edges % getEdgeVertexIdxs(faceEdgeIdxs(j))

        do k = j + 1, size(faceEdgeIdxs)
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
  pure subroutine computeEdgeAngularSectors(self, elements, faces, vertices, edges)
    class(patchSearchAcceleration), intent(in)   :: self
    type(elementShelf), intent(in)               :: elements
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    type(edgeShelf), intent(inout)               :: edges
    integer(shortInt)                            :: i, idx, j, k, nAngularSectors, nElements
    integer(shortInt), dimension(2)              :: commonEdgeVertexIdxs, edgeVertexIdxs
    integer(shortInt), dimension(:), allocatable :: commonEdgeIdxs, edgeElementIdxs, sharingEdgeIdxs, finalIdxsArray
    real(defReal), dimension(2)                  :: temp
    real(defReal), dimension(3)                  :: commonEdgeUnitVector, edgeUnitVector, localBasis1, localBasis2
    real(defReal), dimension(:, :), allocatable  :: anglesArray, finalAnglesArray

    ! Loop through all edges.
    do i = 1, edges % getSize()
      !---------------------------------------------------------------------------------------------------------------
      ! construct 2D local coordinates system (localBasis1,localBasis2) on the plane whose normal is given as the 
      ! current edge's unit vector and contains the second vertex of the edge.
      !---------------------------------------------------------------------------------------------------------------
      edgeUnitVector = edges % getEdgeUnitVector(i)
      
      ! construct localBasis1
      if(abs(edgeUnitVector(1)) <= abs(edgeUnitVector(2)) .and. abs(edgeUnitVector(1)) <= abs(edgeUnitVector(3))) then
        localBasis1 = [ZERO, edgeUnitVector(3), -edgeUnitVector(2)]

      elseif(abs(edgeUnitVector(2)) <= abs(edgeUnitVector(3))) then 
        localBasis1 = [-edgeUnitVector(3), ZERO, edgeUnitVector(1)]

      else
        localBasis1 = [edgeUnitVector(2), -edgeUnitVector(1), ZERO]

      end if
      localBasis1 = localBasis1 / norm2(localBasis1)

      ! construct localBasis2 (cross product gives the normalised vector)
      localBasis2 = crossProduct(edgeUnitVector, localBasis1)

      ! store localBasis1 and localBasis2 to each associated edge
      call edges % setEdgeLocalBasis1(i, localBasis1)
      call edges % setEdgeLocalBasis2(i, localBasis2)

      !---------------------------------------------------------------------------------------------------------------
      ! construct (unsorted) arrays for angles and associated elementIdxs
      !---------------------------------------------------------------------------------------------------------------
      ! retrieve relevant information
      edgeElementIdxs = edges % getEdgeElementIdxs(i)
      edgeVertexIdxs = edges % getEdgeVertexIdxs(i)
      nElements = size(edgeElementIdxs)

      ! Retrieve all edges connected to the second vertex of the current edge.
      sharingEdgeIdxs = vertices % getVertexEdgeIdxs(edgeVertexIdxs(2))

      ! initialise arrays for angle and face index
      if(allocated(anglesArray)) deallocate(anglesArray)
      allocate(anglesArray(nElements, 2))

      ! Loop through all the elements containing the current edge.
      nAngularSectors = 0
      do j = 1, nElements
        ! Retrieve edges in the current element then find common edges with those sharing the second vertex of the current edge.
        commonEdgeIdxs = findCommon(sharingEdgeIdxs, elements % getElementEdgeIdxs(edgeElementIdxs(j)))

        ! Reset idx = 0 then loop through all common edges.
        idx = 0
        do k = 1, size(commonEdgeIdxs)
          ! One of these edges is the original edge itself so skip it.
          if(i == commonEdgeIdxs(k)) cycle
          idx = idx + 1

          ! Get the vertices in the current edge then compute pseudo-angle.
          commonEdgeVertexIdxs = edges % getEdgeVertexIdxs(commonEdgeIdxs(k))
          commonEdgeUnitVector = edges % getEdgeUnitVector(commonEdgeIdxs(k))
          if(commonEdgeVertexIdxs(1) /= edgeVertexIdxs(2)) commonEdgeUnitVector = -commonEdgeUnitVector
          temp(idx) = computePseudoAngle(commonEdgeUnitVector, localBasis1, localBasis2)

        end do
        anglesArray(j, :) = [minval(temp), maxval(temp)]
        nAngularSectors = nAngularSectors + merge(1, 2, abs(anglesArray(j, 2) - anglesArray(j, 1)) <= TWO)

      end do

      ! Now split angular intervals which are outside the intervals [-2, 0] or [0, 2].
      if(allocated(finalAnglesArray)) deallocate(finalAnglesArray)
      if(allocated(finalIdxsArray)) deallocate(finalIdxsArray)
      allocate(finalAnglesArray(nAngularSectors, 2), finalIdxsArray(nAngularSectors))
      idx = 0
      do j = 1, nElements
        if(abs(anglesArray(j, 2) - anglesArray(j, 1)) <= TWO) then
          idx = idx + 1
          finalAnglesArray(idx, :) = anglesArray(j, :)
          finalIdxsArray(idx) = edgeElementIdxs(j)

        else
          idx = idx + 1
          finalAnglesArray(idx, :) = [-TWO, anglesArray(j, 1)]
          finalIdxsArray(idx) = edgeElementIdxs(j)

          idx = idx + 1
          finalAnglesArray(idx, :) = [anglesArray(j, 2), TWO]
          finalIdxsArray(idx) = edgeElementIdxs(j)

        end if

      end do

      ! pass and set elementIdxsArray and anglesArray to each corresponding edge
      call edges % setEdgeAnglesArray(i, finalAnglesArray)
      call edges % setEdgeElementIdxsArray(i, finalIdxsArray)

    end do

  end subroutine computeEdgeAngularSectors

  !!
  !!
  !!
  pure subroutine computeGeometricParameters(self, elements, faces, vertices, edges)
    class(patchSearchAcceleration), intent(inout) :: self
    type(elementShelf), intent(in)                :: elements
    type(faceShelf), intent(in)                   :: faces
    type(vertexShelf), intent(in)                 :: vertices
    type(edgeShelf), intent(inout)                :: edges
    integer(shortInt)                             :: i
    integer(shortInt), dimension(2)               :: edgeVertexIdxs
    real(defReal)                                 :: edgeLength
    real(defReal), dimension(3)                   :: edgeVector

    ! Compute minimum edge length.
    self % minimumEdgeLength = INF
    do i = 1, edges % getSize()
      edgeVertexIdxs = edges % getEdgeVertexIdxs(i)
      edgeVector = vertices % getVertexCoordinates(edgeVertexIdxs(2)) - vertices % getVertexCoordinates(edgeVertexIdxs(1))
      edgeLength = norm2(edgeVector)

      call edges % setEdgeUnitVector(i, edgeVector / edgeLength)
      call edges % setEdgeLength(i, edgeLength)
      self % minimumEdgeLength = min(self % minimumEdgeLength, edgeLength)

    end do

    ! Compute minimum angle.
    self % minimumAngle = acos(max(self % computeCosineMaximumDihedralAngle(edges, elements, faces), &
                                   self % computeCosineMaximumFaceAngle(edges, faces)))

  end subroutine computeGeometricParameters

  !!
  !!
  !!
  subroutine findHostElement(self, vertices, edges, faces, elements, coords)
    class(patchSearchAcceleration), intent(in)   :: self
    class(vertexShelf), intent(in)               :: vertices
    class(edgeShelf), intent(in)                 :: edges
    type(faceShelf), intent(in)                  :: faces
    type(elementShelf), intent(in)               :: elements
    type(coord), intent(inout)                   :: coords
    integer(shortInt)                            :: elementIdx, n
    real(defReal), dimension(3)                  :: r

    ! retrieve the coordinates of neutron
    r = coords % getPositionToNudge()

    ! Return immediately if the particle is outside the mesh bounds.
    if(any(r < self % meshBounds(1:3)) .or. any(self % meshBounds(4:6) < r)) return

    ! Search subgrids or not depending on whether multiple layers have been defined.
    elementIdx = 0
    n = 0
    call self % searchGrids(r, edges, elements, faces, vertices, elementIdx, n)

    ! Update coords if a valid element has been found.
    if(0 < elementIdx) then
      call coords % setElementIdx(elementIdx)
      call coords % setParentElementIdx(elements % getElementParentIdx(elementIdx))
      call coords % setLocalId(elements % getElementLocalId(elementIdx))

    end if

  end subroutine findHostElement

  !!
  !!
  !!
  subroutine init(self, dict, vertices, edges, faces, elements)
    class(patchSearchAcceleration), intent(inout)  :: self
    type(dictionary), intent(in)                   :: dict
    type(vertexShelf), intent(in)                  :: vertices
    type(edgeShelf), intent(inout)                 :: edges
    type(faceShelf), intent(inout)                 :: faces
    type(elementShelf), intent(in)                 :: elements
    integer(shortInt)                              :: i, nFaces
    real(defReal)                                  :: coarsestLayerSpacing, ratio, sineAlpha, sineHalfAlpha
    real(defReal), dimension(3)                    :: extraDistances
    real(defReal), dimension(6)                    :: gridBounds
    type(CartesianGrid), dimension(:), allocatable :: temp
    character(*), parameter                        :: HERE = 'init (patchSearchAcceleration_class.f90)'

    ! Retrieve number of layers from dictionary and allocate memory.
    call dict % getOrDefault(self % nLayers, 'nLayers', 1)
    if(self % nLayers < 1) call fatalError(HERE, 'Number of layers must be at least 1. Is: '//numToChar(self % nLayers)//'.')
    allocate(self % spacings(self % nLayers))

    ! Compute geometric parameters and angular sectors for edges.
    call self % computeGeometricParameters(elements, faces, vertices, edges)
    call self % computeEdgeAngularSectors(elements, faces, vertices, edges)

    ! Set face constants and compute caches.
    nFaces = faces % getSize()
    allocate(self % faceCaches(nFaces))
    do i = 1, nFaces
      call faces % setFaceConst(i, -dot_product(faces % getFaceNormal(i), faces % getFaceCentroid(i)))
      self % faceCaches(i) = faces % computeFaceSATData(i, edges, vertices)

    end do

    ! Pre-compute sin(alpha) and sin(alpha / 2) then compute wStar and targetDistance.
    sineAlpha = sin(self % minimumAngle)
    sineHalfAlpha = sin(HALF * self % minimumAngle)
    self % wStar = self % minimumEdgeLength * min(HALF, sineAlpha) * (ONE - epsilon(ONE))
    self % targetDistance = (self % wStar / (ONE + sin(self % minimumAngle))) ** 2
    
    ! Now compute spacing in the finest layer and its inverse.
    self % spacings(self % nLayers) = &
    TWO * self % wStar * sineAlpha * sineHalfAlpha / (sqrt(THREE) * (ONE + sineAlpha) * (ONE + sineHalfAlpha))

    ! Set mesh bounds and compute bounds for the root grid.
    self % meshBounds = vertices % getExtremalCoordinates()
    if(1 < self % nLayers) then
      ! Estimate spacing in the coarsest layer. Use 100 subdivisions along the dimension of maximum extent in the mesh
      ! to avoid memory explosion.
      coarsestLayerSpacing = maxval(self % meshBounds(4:6) - self % meshBounds(1:3)) / 100

      ! Compute spacing ratio between the coarsest and finest layers. Use a minimum ratio of TWO.
      ratio = max(real(nint((coarsestLayerSpacing / self % spacings(self % nLayers)) ** (ONE / (self % nLayers - 1))), defReal), &
                  TWO)

      ! Now propagate spacings from the finest layer up to the coarsest layer.
      do i = self % nLayers - 1, 1, -1
        self % spacings(i) = self % spacings(i + 1) * ratio

      end do

    end if
    extraDistances = HALF * (self % spacings(1) - mod(self % meshBounds(4:6) - self % meshBounds(1:3), self % spacings(1)))
    gridBounds(1:3) = self % meshBounds(1:3) - extraDistances
    gridBounds(4:6) = self % meshBounds(4:6) + extraDistances

    ! Initialise root grid and map it.
    call self % rootGrid % init(self % spacings(1), gridBounds)
    call self % mapGrid([(i, i = 1, elements % getSize())], [(i, i = 1, faces % getSize())], self % nLayers == 1, edges, &
                        elements, faces, vertices, self % rootGrid)

    ! Refine root grid if needed.
    if(1 < self % nLayers) then
      call self % refineGrid(1, 0, edges, elements, faces, vertices)
      if(self % nSubGrids < size(self % subGrids)) then
        allocate(temp(self % nSubGrids))
        temp = self % subGrids(1:self % nSubGrids)
        call move_alloc(temp, self % subGrids)

      end if

    end if

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(patchSearchAcceleration), intent(inout) :: self
    integer(shortInt)                             :: i

    ! Local.
    self % nLayers = 0
    self % nSubGrids = 0
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
    if(allocated(self % faceCaches)) deallocate(self % faceCaches)

  end subroutine kill

  !!
  !!
  !!
  subroutine mapGrid(self, elementIdxs, faceIdxs, isFinestLayer, edges, elements, faces, vertices, grid)
    class(patchSearchAcceleration), intent(in)   :: self
    integer(shortInt), dimension(:), intent(in)  :: elementIdxs, faceIdxs
    logical(defBool), intent(in)                 :: isFinestLayer
    type(edgeShelf), intent(in)                  :: edges
    type(elementShelf), intent(in)               :: elements
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    type(CartesianGrid), intent(inout)           :: grid
    integer(shortInt)                            :: i, j, k, l
    integer(shortInt), dimension(3)              :: nCells
    integer(shortInt), dimension(6)              :: cellIdxs
    integer(shortInt), dimension(:), allocatable :: elementFaceIdxs, elementVertexIdxs, faceVertexIdxs
    real(defReal)                                :: extraDistance, spacing
    real(defReal), dimension(3)                  :: centroid
    real(defReal), dimension(6)                  :: bounds

    ! Retrieve grid spacing and bounds.
    nCells = grid % getCellsNumber()
    spacing = grid % getSpacing()
    bounds = grid % getBounds()

    ! Loop over all faces.
    do i = 1, size(faceIdxs)
      extraDistance = HALF * sum(abs(self % faceCaches(faceIdxs(i)) % faceNormal)) * spacing

      ! Generate the indices of the cells contained in the current face's AABB.
      faceVertexIdxs = faces % getFaceVertexIdxs(faceIdxs(i))
      cellIdxs = grid % constructCellIdxs(faceVertexIdxs, vertices)

      ! Loop over all cells.
      do l = max(1, cellIdxs(3)), min(nCells(3), cellIdxs(6))
        centroid(3) = bounds(3) + spacing * (l - HALF)
        do k = max(1, cellIdxs(2)), min(nCells(2), cellIdxs(5))
          centroid(2) = bounds(2) + spacing * (k - HALF)
          do j = max(1, cellIdxs(1)), min(nCells(1), cellIdxs(4))
            centroid(1) = bounds(1) + spacing * (j - HALF)
            ! Now test current cell for intersection with the current face.
            call grid % testCellFaceIntersection(faceIdxs(i), j, k, l, extraDistance, centroid, self % faceCaches(faceIdxs(i)))

            ! If we are at the finest layer, map the cell.
            if(isFinestLayer) call grid % mapCell(j, k, l, self % targetDistance, centroid, edges, faces, vertices)

          end do

        end do

      end do

    end do

    ! Now loop over all elements.
    do i = 1, size(elementIdxs)
      ! Generate the indices of the cells contained in the current element's AABB.
      elementFaceIdxs = elements % getElementFaceIdxs(elementIdxs(i))
      elementVertexIdxs = elements % getElementVertexIdxs(elementIdxs(i))
      cellIdxs = grid % constructCellIdxs(elementVertexIdxs, vertices)

      ! Loop over all cells.
      do l = max(1, cellIdxs(3)), min(nCells(3), cellIdxs(6))
        centroid(3) = bounds(3) + spacing * (l - HALF)
        do k = max(1, cellIdxs(2)), min(nCells(2), cellIdxs(5))
          centroid(2) = bounds(2) + spacing * (k - HALF)
          do j = max(1, cellIdxs(1)), min(nCells(1), cellIdxs(4))
            if(grid % isCellUnprocessed(j, k, l)) then
              centroid(1) = bounds(1) + spacing * (j - HALF)
              call grid % testCellElementInclusion(elementIdxs(i), j, k, l, elementFaceIdxs, centroid, elements, faces)

            end if

          end do

        end do

      end do

    end do

  end subroutine mapGrid

  !!
  !!
  !!
  recursive subroutine refineGrid(self, depth, subGridIdx, edges, elements, faces, vertices)
    class(patchSearchAcceleration), target, intent(inout) :: self
    integer(shortInt), intent(in)                         :: depth, subGridIdx
    type(edgeShelf), intent(in)                           :: edges
    type(elementShelf), intent(in)                        :: elements
    type(faceShelf), intent(in)                           :: faces
    type(vertexShelf), intent(in)                         :: vertices
    integer(shortInt)                                     :: i, j, k, newDepth
    integer(shortInt), dimension(3)                       :: nCells
    integer(shortInt), dimension(:), allocatable          :: cellIntersectedFaceIdxs
    real(defReal)                                         :: gridSpacing
    real(defReal), dimension(6)                           :: cellBounds, gridBounds
    type(CartesianGrid), pointer                          :: currentGridPtr

    if(self % nLayers <= depth) return

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
          if(currentGridPtr % isCellSimple(i, j, k)) cycle
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
          call self % mapGrid(faces % getFaceElementIdxs(cellIntersectedFaceIdxs), cellIntersectedFaceIdxs, &
                              self % nLayers == newDepth, edges, elements, faces, vertices, self % subGrids(self % nSubGrids))
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
  recursive subroutine searchGrids(self, r, edges, elements, faces, vertices, elementIdx, n)
    class(patchSearchAcceleration), target, intent(in) :: self
    real(defReal), dimension(3), intent(in)            :: r
    type(edgeShelf), intent(in)                        :: edges
    type(elementShelf), intent(in)                     :: elements
    type(faceShelf), intent(in)                        :: faces
    type(vertexShelf), intent(in)                      :: vertices
    integer(shortInt), intent(inout)                   :: elementIdx, n
    integer(shortInt)                                  :: edgeIdx, i, subGridIdx
    integer(shortInt), dimension(2)                    :: edgeVertexIdxs
    integer(shortInt), dimension(3)                    :: hostCellIdxs
    integer(shortInt), dimension(:), allocatable       :: elementIdxsArray, intersectedFaceIdxs, potentialElementIdxs
    real(defReal)                                      :: thetaHat
    real(defReal), dimension(3)                        :: displacementVector, vertexCoords
    real(defReal), dimension(:, :), allocatable        :: angularSectorsArray
    type(CartesianGrid), pointer                       :: currentGridPtr

    ! Start at the root grid then search until we hit a terminal cell.
    currentGridPtr => self % rootGrid
    do
      hostCellIdxs = currentGridPtr % findHostCellIdxs(r)
      subGridIdx = currentGridPtr % getCellSubGridIdx(hostCellIdxs(1), hostCellIdxs(2), hostCellIdxs(3))
      if(subGridIdx == 0) exit
      currentGridPtr => self % subGrids(subGridIdx)
      
    end do

    ! Check if terminal cell is fully outside and return immediately if so.
    if(currentGridPtr % isCellOutside(hostCellIdxs(1), hostCellIdxs(2), hostCellIdxs(3))) return

    ! Check if terminal cell is fully inside an element and return immediately if so.
    elementIdx = currentGridPtr % getCellElementIdx(hostCellIdxs(1), hostCellIdxs(2), hostCellIdxs(3))
    if(0 < elementIdx) return

    ! For multi-layered Patch-Search, check if terminal cell only intersects with a single face. In this case, perform an 
    ! element inclusion test on the elements sharing this face and return.
    if(1 < self % nLayers) then
      intersectedFaceIdxs = currentGridPtr % getCellIntersectedFaceIdxs(hostCellIdxs(1), hostCellIdxs(2), hostCellIdxs(3))
      if(size(intersectedFaceIdxs) == 1) then
        potentialElementIdxs = faces % getFaceElementIdxs(intersectedFaceIdxs(1))
        if(ZERO < dot_product(faces % getFaceCentroid(intersectedFaceIdxs(1)) - r, &
           faces % getFaceNormal(intersectedFaceIdxs(1)))) then
          elementIdx = minval(potentialElementIdxs)

        elseif(.not. faces % getFaceIsBoundary(intersectedFaceIdxs(1))) then
          elementIdx = maxval(potentialElementIdxs)

        end if
        return

      end if

    end if

    ! Else, begin Patch-Search procedure.
    edgeIdx = currentGridPtr % getCellEdgeIdx(hostCellIdxs(1), hostCellIdxs(2), hostCellIdxs(3))
    if(edgeIdx == 0) then
      vertexCoords = vertices % getVertexCoordinates(currentGridPtr % getCellVertexIdx(hostCellIdxs(1), hostCellIdxs(2), &
                                                                                       hostCellIdxs(3)))
      displacementVector = r - vertexCoords
      call self % searchGrids(vertexCoords + self % wStar * displacementVector / norm2(displacementVector), edges, elements, &
                              faces, vertices, elementIdx, n)

    else
      ! Compute pseudo-angle.
      edgeVertexIdxs = edges % getEdgeVertexIdxs(edgeIdx) 
      thetaHat = computePseudoAngle(r - vertices % getVertexCoordinates(edgeVertexIdxs(2)), edges % getEdgeLocalBasis1(edgeIdx), &
                                    edges % getEdgeLocalBasis2(edgeIdx))

      ! Search and return.
      elementIdxsArray = edges % getEdgeElementIdxsArray(edgeIdx)
      angularSectorsArray = edges % getEdgeAnglesArray(edgeIdx)
      do i = 1, size(angularSectorsArray, 1)
        if(angularSectorsArray(i, 1) <= thetaHat .and. thetaHat <= angularSectorsArray(i, 2)) then
          elementIdx = elementIdxsArray(i)
          return

        end if

      end do

    end if

  end subroutine searchGrids

end module patchSearchAcceleration_class