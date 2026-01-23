module uniformCartesianGrid_class

  use CartesianCell_inter,         only : CartesianCell
  use cartesianGenericProcedures,  only : calculateAvgEdgeLength, constructAABB
  use CartesianGrid_inter,         only : CartesianGrid, init_super => init
  use dictionary_class,            only : dictionary
  use edgeShelf_class,             only : edgeShelf
  use elementShelf_class,          only : elementShelf
  use errors_mod,                  only : fatalError
  use face_inter,                  only : faceSATData
  use faceShelf_class,             only : faceShelf
  use genericProcedures,           only : crossProduct, findCommon
  use numPrecision
  use terminalCartesianCell_class, only : terminalCartesianCell
  use vertexShelf_class,           only : vertexShelf

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(CartesianGrid) :: uniformCartesianGrid
    private
  contains
    procedure :: allocateCell
    procedure :: constructMapping
    procedure :: findHostElementIdx
    procedure :: init
    procedure :: sortAngles
  end type uniformCartesianGrid

contains
  !!
  !!
  !!
  subroutine allocateCell(self, cell)
    class(uniformCartesianGrid), intent(in)        :: self
    class(CartesianCell), allocatable, intent(out) :: cell

    ! Allocate.
    allocate(terminalCartesianCell :: cell)

  end subroutine allocateCell

  !!
  !!
  !!
  subroutine constructMapping(self, edges, elements, faces, vertices)
    class(uniformCartesianGrid), intent(inout)   :: self
    type(edgeShelf), intent(inout)               :: edges
    type(elementShelf), intent(in)               :: elements
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    integer(shortInt)                            :: i, j, k, l, debugEdgeIdx
    integer(shortInt), dimension(6)              :: cellIdxs
    integer(shortInt), dimension(:), allocatable :: elementFaceIdxs, elementVertexIdxs, faceVertexIdxs
    real(defReal)                                :: extraDistance, targetDistance, reciprocalSpacing, spacing
    real(defReal), dimension(2)                  :: tempIdxs
    real(defReal), dimension(3)                  :: centroid, minimumGridBounds, realIdxs
    type(faceSATData)                            :: cache

    ! Pre-compute constants.
    minimumGridBounds = self % getMinimumGridBounds()
    reciprocalSpacing = self % getReciprocalSpacing()
    spacing = self % getSpacing()
    targetDistance = (self % getWStar() / (ONE + sin(self % getMinimumAngle()))) ** 2

    ! Loop over all faces in the mesh.
    do i = 1, faces % getSize()
      ! Generate SAT data cache for current face.
      cache = faces % computeFaceSATData(i, edges, vertices)
      extraDistance = HALF * sum(abs(cache % faceNormal)) * spacing

      ! Generate the indices of the cells contained in the current face's AABB.
      faceVertexIdxs = faces % getFaceVertexIdxs(i)
      cellIdxs = constructAABB(vertices, faceVertexIdxs, minimumGridBounds, reciprocalSpacing)

      ! Loop over all cells.
      do l = cellIdxs(3), cellIdxs(6)
        centroid(3) = minimumGridBounds(3) + spacing * (l - HALF)
        do k = cellIdxs(2), cellIdxs(5)
          centroid(2) = minimumGridBounds(2) + spacing * (k - HALF)
          do j = cellIdxs(1), cellIdxs(4)
            centroid(1) = minimumGridBounds(1) + spacing * (j - HALF)
            ! Now test current cell for intersection with the current face.
            call self % testCellFaceIntersection(i, j, k, l, extraDistance, targetDistance, centroid, edges, cache, faces, &
                                                 vertices)

          end do

        end do

      end do

    end do

    ! Now loop over all elements in the mesh.
    do i = 1, elements % getSize()
      ! Generate the indices of the cells contained in the current face's AABB.
      elementFaceIdxs = elements % getElementFaceIdxs(i)
      elementVertexIdxs = elements % getElementVertexIdxs(i)
      cellIdxs = constructAABB(vertices, elementVertexIdxs, minimumGridBounds, reciprocalSpacing)

      ! Loop over all cells.
      do l = cellIdxs(3), cellIdxs(6)
        tempIdxs(2) = l - HALF
        do k = cellIdxs(2), cellIdxs(5)
          tempIdxs(1) = k - HALF
          do j = cellIdxs(1), cellIdxs(4)
            if(self % isCellUnprocessed(j, k, l)) then
              centroid(1) = minimumGridBounds(1) + spacing * (j - HALF)
              centroid(2) = minimumGridBounds(2) + spacing * tempIdxs(1)
              centroid(3) = minimumGridBounds(3) + spacing * tempIdxs(2)
              call self % testCellElementInclusion(i, j, k, l, elementFaceIdxs, centroid, elements, faces)

            end if

          end do

        end do

      end do

    end do

  end subroutine constructMapping

  !!
  !!
  !!
  recursive subroutine findHostElementIdx(self, reciprocalSpacing, minimumBounds, r, edges, vertices, edgeIdx, hostElementIdx)
    class(uniformCartesianGrid), intent(in)      :: self
    real(defReal), intent(in)                    :: reciprocalSpacing
    real(defReal), dimension(3), intent(in)      :: minimumBounds, r
    type(edgeShelf), intent(in)                  :: edges
    type(vertexShelf), intent(in)                :: vertices
    integer(shortInt), intent(inout)             :: edgeIdx, hostElementIdx
    integer(shortInt)                            :: cellToElementIdx, i
    integer(shortInt), dimension(2)              :: currEdgeVertexIdxs
    integer(shortInt), dimension(3)              :: cellIdxs
    integer(shortInt), dimension(:), allocatable :: elementIdxsArray
    real(defReal)                                :: thetaHat, xLocalCoord, yLocalCoord
    real(defReal), dimension(3)                  :: displacementVector, rLocalCoords, vertexCoords
    real(defReal), dimension(:, :), allocatable  :: angularSectorsArray

    ! Find cartesian cell indices. In case the current cell lies outside the computational domain for the unstructured mesh, 
    ! return.
    cellIdxs = ceiling((r - minimumBounds) * reciprocalSpacing)
    if(self % isCellOutside(cellIdxs(1), cellIdxs(2), cellIdxs(3))) return

    ! retrieve element index from chi mapping.
    cellToElementIdx = self % getCellToElementIdx(cellIdxs(1), cellIdxs(2), cellIdxs(3))

    ! if element index is valid (the current cell, characterised by "cellIdxs", is fully contained within that element)
    if(0 < cellToElementIdx) then
      hostElementIdx = cellToElementIdx
      return

    end if
      
    ! otherwise, the current cell intersects with either face(s) or edge(s). Start patch searching.
    edgeIdx = self % getCellToEdgeIdx(cellIdxs(1), cellIdxs(2), cellIdxs(3))
    if (edgeIdx == 0) then
      ! push the coordinates away from the current vertex (= phi)
      ! (needs to be changed) (possible improvement/acceleration for the rest of the subroutine below?)
      vertexCoords = vertices % getVertexCoordinates(self % getCellToVertexIdx(cellIdxs(1), cellIdxs(2), cellIdxs(3)))
      displacementVector = r - vertexCoords
      call self % findHostElementIdx(reciprocalSpacing, minimumBounds, &
                                     vertexCoords + self % getWStar() * displacementVector / norm2(displacementVector), &
                                     edges, vertices, edgeIdx, hostElementIdx)
      if(0 < hostElementIdx) return

    end if
  
    ! calculate pseudo angle
    currEdgeVertexIdxs = edges % getEdgeVertexIdxs(edgeIdx) 
    rLocalCoords = r - vertices % getVertexCoordinates(currEdgeVertexIdxs(2))
    xLocalCoord = dot_product(rLocalCoords, edges % getEdgeLocalBasis1(edgeIdx))
    yLocalCoord = dot_product(rLocalCoords, edges % getEdgeLocalBasis2(edgeIdx))
    thetaHat = sign(ONE - xLocalCoord / (abs(xLocalCoord) + abs(yLocalCoord)), yLocalCoord)

    ! perform binary search and return index pointer
    elementIdxsArray = edges % getEdgeElementIdxsArray(edgeIdx)
    angularSectorsArray = edges % getEdgeAnglesArray(edgeIdx)

    do i = 1, size(angularSectorsArray, 1)
      if(angularSectorsArray(i, 1) <= thetaHat .and. thetaHat <= angularSectorsArray(i, 2)) then
        hostElementIdx = elementIdxsArray(i)
        return

      end if

    end do

  end subroutine findHostElementIdx

  !!
  !!
  !!
  subroutine init(self, dict, elements, vertices, edges, faces)
    class(uniformCartesianGrid), intent(inout) :: self
    type(dictionary), intent(in)               :: dict
    type(elementShelf), intent(in)             :: elements
    type(vertexShelf), intent(in)              :: vertices
    type(edgeShelf), intent(inout)             :: edges
    type(faceShelf), intent(inout)             :: faces
    integer(shortInt)                          :: i, j, k
    integer(shortInt), dimension(3)            :: cellsNumber
    real(defReal)                              :: halfSpacingLessExtraRoom, reciprocalSpacing, spacing
    real(defReal), dimension(3)                :: gridMaximumBounds, gridMinimumBounds, meshMaximumBounds, meshMinimumBounds

    ! Initialise superclass.
    call init_super(self, dict, elements, vertices, edges, faces)

    reciprocalSpacing = self % getReciprocalSpacing()
    spacing = self % getSpacing()
    meshMaximumBounds = self % getMaximumMeshBounds()
    meshMinimumBounds = self % getMinimumMeshBounds()

    ! Compute number of cells in each dimension.
    do i = 1, 3
      halfSpacingLessExtraRoom = HALF * (spacing - mod(meshMaximumBounds(i) - meshMinimumBounds(i), spacing))
      gridMinimumBounds(i) = meshMinimumBounds(i) - halfSpacingLessExtraRoom
      gridMaximumBounds(i) = meshMaximumBounds(i) + halfSpacingLessExtraRoom
      cellsNumber(i) = nint((gridMaximumBounds(i) - gridMinimumBounds(i)) * reciprocalSpacing)
        
    end do
    call self % setMaximumGridBounds(gridMaximumBounds)
    call self % setMinimumGridBounds(gridMinimumBounds)
    call self % setCellsNumber(cellsNumber)

    ! print cartesian grid parameters and mesh quality
    print*, "----------------------------------------------------"
    print*, "/\/\ Cartesian grid parameters and mesh quality /\/\"
    print*, "Minimum angle            : ", self % getMinimumAngle()
    print*, "Minimum edge length      : ", self % getMinimumEdgeLength()
    print*, "No. of vertices          : ", vertices % getSize()
    print*, "No. of edges             : ", edges % getSize()
    print*, "No. of faces             : ", faces % getSize()
    print*, "No. of elements          : ", elements % getSize()
    print*, "Grid spacing             : ", spacing
    print*, "Grid size in x           : ", cellsNumber(1)
    print*, "Grid size in y           : ", cellsNumber(2)
    print*, "Grid size in z           : ", cellsNumber(3)
    print*, "Grid lower bounds in xyz : ", gridMinimumBounds
    print*, "Grid upper bounds in xyz : ", gridMaximumBounds
    print*, "----------------------------------------------------"

    ! Allocate memory then allocate all cells.
    call self % allocateGrid(cellsNumber)

    !-----------------------------------------------------------------------------------------
    !initialise for patch search
    !-----------------------------------------------------------------------------------------
    call self % constructMapping(edges, elements, faces, vertices)
    call self % sortAngles(edges, elements, faces, vertices)

    print*, "average edge length", calculateAvgEdgeLength(edges)

  end subroutine init

  !!
  !!
  !!
  subroutine sortAngles(self, edges, elements, faces, vertices)
    class(uniformCartesianGrid), intent(inout)   :: self
    type(edgeShelf), intent(inout)               :: edges
    type(elementShelf), intent(in)               :: elements
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    integer(shortInt)                            :: i, idx, j, k, l, pointerIdx, outer2LoopSize, nAngularSectors, nElements
    real(defReal), dimension(2)                  :: temp
    real(defReal), dimension(3)                  :: currEdgeUnitVector, localBasis1, localBasis2, currUnitVector
    integer(shortInt), dimension(:), allocatable :: currFaceEdgeIdxs, currEdgeElementIdxs, currEdgeFaceIdxs, faceIdxsArray, &
                                                    elementIdxsArray, sharingEdgeIdxs, currElementEdgeIdxs, commonEdgeIdxs, &
                                                    elementEdgeVertexIdxs, idxsArray, finalIdxsArray
    integer(shortInt), dimension(2)              :: currEdgeVertexIdxs, currVertexIdxs, face1ElementIdxs, face2ElementIdxs
    real(defReal)                                :: x, y, thetaHat
    real(defReal), dimension(:, :), allocatable  :: anglesArray, finalAnglesArray

    ! Loop through all edges.
    do i = 1, edges % getSize()
      !---------------------------------------------------------------------------------------------------------------
      ! construct 2D local coordinates system (localBasis1,localBasis2) on the plane whose normal is given as the 
      ! current edge's unit vector and contains the second vertex of the edge.
      !---------------------------------------------------------------------------------------------------------------
      currEdgeUnitVector = edges % getEdgeUnitVector(i)
      
      ! construct localBasis1
      if(abs(currEdgeUnitVector(1)) <= abs(currEdgeUnitVector(2)) .and. &
         abs(currEdgeUnitVector(1)) <= abs(currEdgeUnitVector(3))) then
        localBasis1 = [ZERO, currEdgeUnitVector(3), -currEdgeUnitVector(2)]

      elseif(abs(currEdgeUnitVector(2)) <= abs(currEdgeUnitVector(3))) then 
        localBasis1 = [-currEdgeUnitVector(3), ZERO, currEdgeUnitVector(1)]

      else
        localBasis1 = [currEdgeUnitVector(2), -currEdgeUnitVector(1), ZERO]

      end if
      localBasis1 = localBasis1 / norm2(localBasis1)

      ! construct localBasis2 (cross product gives the normalised vector)
      localBasis2 = crossProduct(currEdgeUnitVector, localBasis1)

      ! store localBasis1 and localBasis2 to each associated edge
      call edges % setEdgeLocalBasis1(i, localBasis1)
      call edges % setEdgeLocalBasis2(i, localBasis2)

      !---------------------------------------------------------------------------------------------------------------
      ! construct (unsorted) arrays for angles and associated elementIdxs
      !---------------------------------------------------------------------------------------------------------------
      ! retrieve relevant information
      currEdgeElementIdxs = edges % getEdgeElementIdxs(i)
      currEdgeFaceIdxs = edges % getEdgeFaceIdxs(i)
      currEdgeVertexIdxs = edges % getEdgeVertexIdxs(i)
      nElements = size(currEdgeElementIdxs)

      ! Retrieve all edges connected to the second vertex of the current edge.
      sharingEdgeIdxs = vertices % getVertexEdgeIdxs(currEdgeVertexIdxs(2))

      ! initialise arrays for angle and face index
      if(allocated(anglesArray)) deallocate(anglesArray)
      if(allocated(faceIdxsArray)) deallocate(faceIdxsArray)
      if(allocated(elementIdxsArray)) deallocate(elementIdxsArray)
      allocate(anglesArray(nElements, 2), faceIdxsArray(size(currEdgeFaceIdxs)), elementIdxsArray(nElements))

      ! Loop through all the elements containing the current edge.
      nAngularSectors = 0
      do j = 1, nElements
        ! Retrieve edges in the current element then find common edges with those sharing the second vertex of the current edge.
        currElementEdgeIdxs = elements % getElementEdgeIdxs(currEdgeElementIdxs(j))
        commonEdgeIdxs = findCommon(sharingEdgeIdxs, currElementEdgeIdxs)

        ! Reset idx = 0 then loop through all common edges.
        idx = 0
        do k = 1, size(commonEdgeIdxs)
          ! One of these edges is the original edge itself so skip it.
          if(i == commonEdgeIdxs(k)) cycle
          idx = idx + 1

          ! Get the vertices in the current edge then compute pseudo-angle.
          elementEdgeVertexIdxs = edges % getEdgeVertexIdxs(commonEdgeIdxs(k))
          currUnitVector = edges % getEdgeUnitVector(commonEdgeIdxs(k))
          if(elementEdgeVertexIdxs(1) /= currEdgeVertexIdxs(2)) currUnitVector = -currUnitVector
          x = dot_product(currUnitVector, localBasis1)
          y = dot_product(currUnitVector, localBasis2)
          temp(idx) = sign(ONE - x / (abs(x) + abs(y)), y)

        end do
        anglesArray(j, :) = [minval(temp), maxval(temp)]
        nAngularSectors = nAngularSectors + merge(1, 2, abs(anglesArray(j, 2) - anglesArray(j, 1)) <= TWO)

      end do

      ! Now split angular intervals which are outside the intervals [-2, 0] or [0, 2].
      idxsArray = [(j, j = 1, nAngularSectors)]
      if(allocated(finalAnglesArray)) deallocate(finalAnglesArray)
      if(allocated(finalIdxsArray)) deallocate(finalIdxsArray)
      allocate(finalAnglesArray(nAngularSectors, 2), finalIdxsArray(nAngularSectors))
      idx = 0
      do j = 1, nElements
        if(abs(anglesArray(j, 2) - anglesArray(j, 1)) <= TWO) then
          idx = idx + 1
          finalAnglesArray(idx, :) = anglesArray(j, :)
          finalIdxsArray(idx) = currEdgeElementIdxs(j)

        else
          idx = idx + 1
          finalAnglesArray(idx, :) = [-TWO, anglesArray(j, 1)]
          finalIdxsArray(idx) = currEdgeElementIdxs(j)

          idx = idx + 1
          finalAnglesArray(idx, :) = [anglesArray(j, 2), TWO]
          finalIdxsArray(idx) = currEdgeElementIdxs(j)

        end if

      end do

      ! pass and set elementIdxsArray and anglesArray to each corresponding edge
      call edges % setEdgeAnglesArray(i, finalAnglesArray)
      call edges % setEdgeElementIdxsArray(i, finalIdxsArray)

    end do

  end subroutine sortAngles

end module uniformCartesianGrid_class