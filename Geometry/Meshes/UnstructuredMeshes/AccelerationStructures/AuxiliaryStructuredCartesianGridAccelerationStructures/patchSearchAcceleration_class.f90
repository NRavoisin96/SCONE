module patchSearchAcceleration_class

  use accelerationStructure_inter, only : accelerationStructure
  use ASCGAcceleration_inter,      only : ASCGAcceleration, init_super => init
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
  subroutine computeEdgeAngularSectors(self, elements, faces, edges)
    class(patchSearchAcceleration), intent(in)   :: self
    type(elementShelf), intent(in)               :: elements
    type(faceShelf), intent(in)                  :: faces
    type(edgeShelf), intent(inout)               :: edges
    integer(shortInt)                            :: i, idx, j, k, l, nAngularSectors, nElements
    integer(shortInt), dimension(2)              :: edgeVertexIdxs, faceEdgeVertexIdxs
    integer(shortInt), dimension(:), allocatable :: edgeElementIdxs, edgeFaceIdxs, elementFaceIdxs, &
                                                    faceEdgeIdxs, finalIdxsArray
    real(defReal), dimension(2)                  :: temp
    real(defReal), dimension(3)                  :: boundingEdgeUnitVector, edgeUnitVector, localBasis1, localBasis2
    real(defReal), dimension(:, :), allocatable  :: anglesArray, finalAnglesArray
    character(*), parameter                      :: HERE = 'computeEdgeAngularSectors (patchSearchAcceleration_class.f90)'

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
      edgeFaceIdxs = edges % getEdgeFaceIdxs(i)
      edgeVertexIdxs = edges % getEdgeVertexIdxs(i)
      nElements = size(edgeElementIdxs)

      ! initialise arrays for angle and face index
      if(allocated(anglesArray)) deallocate(anglesArray)
      allocate(anglesArray(nElements, 2))

      ! Loop through all the elements containing the current edge.
      nAngularSectors = 0
      do j = 1, nElements
        ! Retrieve the indices of the faces in the current element.
        elementFaceIdxs = elements % getElementFaceIdxs(edgeElementIdxs(j))

        ! Reset idx = 0 and temp = ZERO then loop through all faces.
        idx = 0
        temp = ZERO
        do k = 1, size(elementFaceIdxs)
          ! Skip this face if it does not contain the current edge.
          if(.not. any(edgeFaceIdxs == elementFaceIdxs(k))) cycle
          idx = idx + 1

          ! Retrieve the edge in the current face that contains the second vertex of the current edge.
          faceEdgeIdxs = faces % getFaceEdgeIdxs(elementFaceIdxs(k))

          ! Loop over all edges in the face.
          do l = 1, size(faceEdgeIdxs)
            ! Skip this edge if it is the current edge.
            if(i == faceEdgeIdxs(l)) cycle

            ! Check if this edge contains the second vertex of the current edge.
            faceEdgeVertexIdxs = edges % getEdgeVertexIdxs(faceEdgeIdxs(l))
            if(any(faceEdgeVertexIdxs == edgeVertexIdxs(2))) then
              ! Compute pseudo-angle.
              boundingEdgeUnitVector = edges % getEdgeUnitVector(faceEdgeIdxs(l))
              if(faceEdgeVertexIdxs(1) /= edgeVertexIdxs(2)) boundingEdgeUnitVector = -boundingEdgeUnitVector
              temp(idx) = computePseudoAngle(boundingEdgeUnitVector, localBasis1, localBasis2)
              exit

            end if

          end do

        end do
        ! Call fatalError if the current element didn't contribute exactly two pseudo-angles.
        if(idx /= 2) call fatalError(HERE, 'Element '//numToChar(edgeElementIdxs(j))//' did not contribute two &
                                           &pseudo-angles for edge '//numToChar(i)//'.')

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
  recursive subroutine findHostElementIdx(self, u, edges, elements, faces, vertices, elementIdx, r)
    class(patchSearchAcceleration), intent(in) :: self
    real(defReal), dimension(3), intent(in)    :: u
    type(edgeShelf), intent(in)                :: edges
    type(elementShelf), intent(in)             :: elements
    type(faceShelf), intent(in)                :: faces
    type(vertexShelf), intent(in)              :: vertices
    integer(shortInt), intent(inout)           :: elementIdx
    real(defReal), dimension(3), intent(inout) :: r
    integer(shortInt)                          :: edgeIdx
    real(defReal), dimension(3)                :: displacementVector, rPrime, vertexCoords
    type(CartesianCell), pointer               :: terminalCellPtr

    terminalCellPtr => self % searchGrids(r)

    ! Check if terminal cell is fully inside an element and return immediately if so.
    elementIdx = terminalCellPtr % getElementIdx()
    if(0 < elementIdx) then
      return

    end if

    if(terminalCellPtr % isOutside()) then
      return

    end if

    ! For multi-layered Patch-Search, check if terminal cell only intersects with a single face. In this case, perform an 
    ! element inclusion test on the elements sharing this face and return.
    if(1 < self % getDepth() .and. self % getSingleFaceShortcut() .and. terminalCellPtr % intersectsOnlyOneFace()) then
      call faces % testFaceHalfSpace(terminalCellPtr % getFirstIntersectedFaceIdx(), r, elementIdx)
      return

    end if

    ! Else, begin Patch-Search procedure.
    edgeIdx = terminalCellPtr % getEdgeIdx()
    if(edgeIdx == 0) then
      vertexCoords = vertices % getVertexCoordinates(terminalCellPtr % getVertexIdx())
      displacementVector = r - vertexCoords
      rPrime = vertexCoords + self % getWStar() * displacementVector / norm2(displacementVector)
      call self % findHostElementIdx(u, edges, elements, faces, vertices, elementIdx, rPrime)

    else
      call edges % findElementIdxFromEdgeAngularSectorSearch(edgeIdx, r, vertices, elementIdx)

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
    call self % computeEdgeAngularSectors(elements, faces, edges)

  end subroutine init

end module patchSearchAcceleration_class