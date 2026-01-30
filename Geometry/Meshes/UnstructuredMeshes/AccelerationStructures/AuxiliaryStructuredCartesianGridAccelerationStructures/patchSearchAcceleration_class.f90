module patchSearchAcceleration_class

  use accelerationStructure_inter, only : accelerationStructure
  use ASCGAcceleration_inter,      only : ASCGAcceleration, init_super => init
  use CartesianCell_class,         only : CartesianCell
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
  type, public, extends(ASCGAcceleration) :: patchSearchAcceleration
    private
  contains
    procedure :: computeEdgeAngularSectors
    procedure :: findHostElementIdx
    procedure :: init
  end type patchSearchAcceleration

contains
  !!
  !!
  !!
  subroutine computeEdgeAngularSectors(self, elements, faces, vertices, edges)
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
  recursive subroutine findHostElementIdx(self, r, edges, elements, faces, vertices, elementIdx)
    class(patchSearchAcceleration), intent(in)   :: self
    real(defReal), dimension(3), intent(in)      :: r
    type(edgeShelf), intent(in)                  :: edges
    type(elementShelf), intent(in)               :: elements
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    integer(shortInt), intent(inout)             :: elementIdx
    integer(shortInt)                            :: edgeIdx, i
    integer(shortInt), dimension(2)              :: edgeVertexIdxs
    integer(shortInt), dimension(:), allocatable :: elementIdxsArray, intersectedFaceIdxs, potentialElementIdxs
    real(defReal)                                :: thetaHat
    real(defReal), dimension(3)                  :: displacementVector, vertexCoords
    real(defReal), dimension(:, :), allocatable  :: angularSectorsArray
    type(CartesianCell), pointer                 :: terminalCellPtr

    ! Start at the root grid then search until we hit a terminal cell.
    terminalCellPtr => self % searchGrids(r)

    ! Check if terminal cell is fully inside an element and return immediately if so.
    elementIdx = terminalCellPtr % getElementIdx()
    if(0 < elementIdx .or. terminalCellPtr % isOutside()) return

    ! For multi-layered Patch-Search, check if terminal cell only intersects with a single face. In this case, perform an 
    ! element inclusion test on the elements sharing this face and return.
    if(1 < self % getDepth()) then
      if(size(intersectedFaceIdxs) == 1) then
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
    edgeIdx = terminalCellPtr % getEdgeIdx()
    if(edgeIdx == 0) then
      vertexCoords = vertices % getVertexCoordinates(terminalCellPtr % getVertexIdx())
      displacementVector = r - vertexCoords
      call self % findHostElementIdx(vertexCoords + self % getWStar() * displacementVector / norm2(displacementVector), edges, &
                                     elements, faces, vertices, elementIdx)

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

  end subroutine findHostElementIdx

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

    call self % setMapCells(.true.)
    call init_super(self, dict, vertices, edges, faces, elements)
    call self % computeEdgeAngularSectors(elements, faces, vertices, edges)

  end subroutine init

end module patchSearchAcceleration_class