module patchMultiAcceleration_class

  use accelerationStructure_inter, only : accelerationStructure
  use coord_class,                 only : coord
  use elementShelf_class,          only : elementShelf
  use faceShelf_class,             only : faceShelf
  use numPrecision
  use vertexShelf_class,           only : vertexShelf
  use edgeShelf_class,             only : edgeShelf
  use cartesianGridCoarsest_class, only : cartesianGridCoarsest
  use cartesianGenericProcedures,  only : binarySearchAngle
  use genericProcedures,            only : fatalError

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(accelerationStructure) :: patchMultiAcceleration
    private
    type(cartesianGridCoarsest)                :: grid
  contains
    procedure                                  :: findHostElement
    procedure                                  :: init
    procedure                                  :: kill
  end type patchMultiAcceleration

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! bit-trick (not saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! !!
  ! !!
  ! !!
  ! subroutine findHostElement(self, vertices, edges, faces, elements, coords)
  !   class(patchMultiAcceleration), intent(in)    :: self
  !   class(vertexShelf), intent(in)               :: vertices
  !   class(edgeShelf), intent(in)                 :: edges
  !   type(faceShelf), intent(in)                  :: faces
  !   type(elementShelf), intent(in)               :: elements
  !   type(coord), intent(inout)                   :: coords
  !   integer(shortInt)                            :: potentialElementIdx, i, edgeIdx, vertexIdx, pointerIdx
  !   integer(shortInt), dimension(3)              :: baseIntegerCoord
  !   real(defReal), dimension(3)                  :: r, v_eCoord, rLocalCoord, dummyVector, phiCoord, gridBounds_min
  !   integer(shortInt), dimension(2)              :: currEdgeVertexIdxs
  !   real(defReal)                                :: thetaHat, xLocalCoord, yLocalCoord, gridSpacingReciprocal
  !   integer(shortInt), dimension(:), allocatable :: elementIdxsArray

  !   ! retrieve the coordinates of neutron
  !   ! (needs to be changed) (needs checking) (is it correct to use "getPositionToNudge" or other coordinates?)
  !   r = coords % getPositionToNudge()
  
  !   !!!!!
  !   ! print*, "-----------------------------------------------"
  !   ! print*, "coord", r
  !   !!!!!

  !   ! check if the position of neutron is inside the catesian grid bounds
  !   if (self % grid % getGridIsOutsideBounds(r)) return

  !   ! retrieve grid properties
  !   gridBounds_min = self % grid % getGridBounds_min()
  !   gridSpacingReciprocal = self % grid % getSpacingReciprocal() ! = grid spacing for the finest layer

  !   ! find cartesian cell indices
  !   ! do i = 1, 3
  !   !   baseIntegerCoord(i) = floor((r(i) - gridBounds_min(i))*(gridSpacingReciprocal))
  !   ! end do
  !   baseIntegerCoord(:) = floor((r(:) - gridBounds_min(:))*(gridSpacingReciprocal))

  !   ! retrieve element index from chi mapping.
  !   potentialElementIdx = self % grid % getGridChi(baseIntegerCoord)

  !   !!!!!
  !   !print*, "indices", baseIntegerCoord
  !   !print*, "PotentialElementIdx", potentialElementIdx
  !   !!!!!
    
  !   ! if element index is valid (the current cell, characterised by "cellIdxs", is fully contained within that element)
  !   if (potentialElementIdx > 0) then
  !     call coords % setElementIdx(potentialElementIdx)

  !   !!!!!
  !   !print*, "indices", baseIntegerCoord
  !   !print*, "PotentialElementIdx", potentialElementIdx
  !   !!!!!
      
  !     call coords % setParentElementIdx(elements % getElementParentIdx(potentialElementIdx))
  !     return

  !   ! in case the current cell lies outside the computational domain for the unstructured mesh, return.
  !   ! this is tested after testing if (chi > 0) because that is the most likely case in terms of the number of the cells
  !   elseif (potentialElementIdx == -1) then
  !     return
      
  !   ! otherwise, the current cell intersects with either face(s) or edge(s). Start patch searching.
  !   else
  !     edgeIdx = self % grid % getGridPhiCapital(baseIntegerCoord)

  !     !!!!!
  !     !print*, "edgeIdx", edgeIdx
  !     !!!!!

  !     if (edgeIdx == 0) then
  !       ! push the coordinates away from the current vertex (= phi)
  !       ! (needs to be changed) (possible improvement/acceleration for the rest of the subroutine below?)
  !       vertexIdx = self % grid % getGridPhi(baseIntegerCoord)
  !       phiCoord = vertices % getVertexCoordinates(vertexIdx)
  !       dummyVector = r - phiCoord
  !       r = phiCoord + (self % grid % getGridWStar())/(norm2(dummyVector))*(dummyVector)

  !       !!!!!
  !       !print*, "Phi", vertexIdx
  !       !print*, "pushed coord", r
  !       !!!!!

  !       ! find updated cartesian cell indices
  !       ! do i = 1, 3
  !       !   baseIntegerCoord(i) = floor((r(i) - gridBounds_min(i))*(gridSpacingReciprocal))
  !       ! end do
  !       baseIntegerCoord(:) = floor((r(:) - gridBounds_min(:))*(gridSpacingReciprocal))

  !       !!!!!
  !       ! print*, "baseCoordPUSHED", baseIntegerCoord
  !       ! print*, self % grid % getGridChi(baseIntegerCoord)
  !       ! print*, self % grid % getGridPhi(baseIntegerCoord)
  !       ! print*, self % grid % getGridPhiCapital(baseIntegerCoord)
  !       !!!!!

  !       ! get updated edge and element index
  !       !!!!!edgeIdx = self % grid % getGridPhiCapital(baseIntegerCoord)
        
  !       !!!!!
  !       !print*, edgeIdx
  !       !!!!!

  !       potentialElementIdx = self % grid % getGridChi(baseIntegerCoord)

  !       ! if pushed coordinate has direct mapping for element index, use that
  !       ! (needs to be changed) (possible acceleration for this and other parts of the subroutine)
  !       ! (needs checking) (is pushed position has direct mapping for element idx, is it guaranteed to lie inside. OW, ">=" not "/=")
  !       if (potentialElementIdx /= 0) then
  !         call coords % setElementIdx(potentialElementIdx)
  !         ! (needs to be changed) (temp:there is no internal subdivision)
  !         !call coords % setParentElementIdx(elements % getElementParentIdx(potentialElementIdx))
  !         call coords % setParentElementIdx(potentialElementIdx)
  !         return
  !       end if

  !       edgeIdx = self % grid % getGridPhiCapital(baseIntegerCoord)

  !     end if
    
  !     ! calculate pseudo angle
  !     currEdgeVertexIdxs = edges % getEdgeVertexIdxs(edgeIdx)
  !     v_eCoord = vertices % getVertexCoordinates(currEdgeVertexIdxs(2))
  !     rLocalCoord = r - v_eCoord
  !     xLocalCoord = dot_product(rLocalCoord, edges % getEdgeLocalBasis1(edgeIdx))
  !     yLocalCoord = dot_product(rLocalCoord, edges % getEdgeLocalBasis2(edgeIdx))
  !     thetaHat = SIGN(1 - (xLocalCoord / (abs(xLocalCoord) + abs(yLocalCoord))), yLocalCoord)

  !     ! perform binary search and return index pointer
  !     ! (needs checking) (index order and mechanics)
  !     !isBoundary = edges % getEdgeIsBoundary(edgeIdx) !!! not needed anymore
  !     pointerIdx = binarySearchAngle(edges % getEdgeAnglesArray(edgeIdx), thetaHat)
  !     elementIdxsArray = edges % getEdgeElementIdxsArray(edgeIdx)
  !     potentialElementIdx = elementIdxsArray(pointerIdx)

  !     ! if the neutron turns out to lie outside the mesh domain, return 
  !     if (potentialElementIdx == 0) return

  !     call coords % setElementIdx(potentialElementIdx)
  !     call coords % setParentElementIdx(elements % getElementParentIdx(potentialElementIdx))

  !   !!!!!
  !   ! print*, self % grid % getGridPhiCapital(baseIntegerCoord)
  !   ! print*, self % grid % getGridPhi(baseIntegerCoord)
  !   ! call fatalError("end","end")
  !   !!!!!

  !     return

  !   end if

  ! end subroutine findHostElement

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! bit-trick (saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !!
  !!
  !!
  subroutine findHostElement(self, vertices, edges, faces, elements, coords)
    class(patchMultiAcceleration), intent(in)      :: self
    class(vertexShelf), intent(in)                 :: vertices
    class(edgeShelf), intent(in)                   :: edges
    type(faceShelf), intent(in)                    :: faces
    type(elementShelf), intent(in)                 :: elements
    type(coord), intent(inout)                     :: coords
    integer(shortInt)                              :: potentialElementIdx, i, edgeIdx, vertexIdx, pointerIdx
    integer(shortInt), dimension(3)                :: baseIntegerCoord
    real(defReal), dimension(3)                    :: r, v_eCoord, rLocalCoord, dummyVector, phiCoord, gridBounds_min
    integer(shortInt), dimension(2)                :: currEdgeVertexIdxs
    real(defReal)                                  :: thetaHat, xLocalCoord, yLocalCoord, gridSpacingReciprocal
    integer(shortInt), dimension(:), allocatable   :: elementIdxsArray
    integer(shortInt), dimension(:,:), allocatable :: cellIdxsMat

    ! retrieve the coordinates of neutron
    ! (needs to be changed) (needs checking) (is it correct to use "getPositionToNudge" or other coordinates?)
    r = coords % getPositionToNudge()
  
    !!!!!
    ! print*, "-----------------------------------------------"
    ! print*, "coord", r
    !!!!!

    ! check if the position of neutron is inside the catesian grid bounds
    if (self % grid % getGridIsOutsideBounds(r)) return

    ! retrieve grid properties
    gridBounds_min = self % grid % getGridBounds_min()
    gridSpacingReciprocal = self % grid % getSpacingReciprocal() ! = grid spacing for the finest layer

    ! find cartesian cell indices
    ! do i = 1, 3
    !   baseIntegerCoord(i) = floor((r(i) - gridBounds_min(i))*(gridSpacingReciprocal))
    ! end do
    baseIntegerCoord(:) = floor((r(:) - gridBounds_min(:))*(gridSpacingReciprocal))

    ! retrieve element index from chi mapping.
    potentialElementIdx = self % grid % getGridChi(baseIntegerCoord, cellIdxsMat)

    !!!!!
    !print*, "indices", baseIntegerCoord
    !print*, "PotentialElementIdx", potentialElementIdx
    !!!!!
    
    ! if element index is valid (the current cell, characterised by "cellIdxs", is fully contained within that element)
    if (potentialElementIdx > 0) then
      call coords % setElementIdx(potentialElementIdx)

    !!!!!
    !print*, "indices", baseIntegerCoord
    !print*, "PotentialElementIdx", potentialElementIdx
    !!!!!
      
      call coords % setParentElementIdx(elements % getElementParentIdx(potentialElementIdx))
      return

    ! in case the current cell lies outside the computational domain for the unstructured mesh, return.
    ! this is tested after testing if (chi > 0) because that is the most likely case in terms of the number of the cells
    elseif (potentialElementIdx == -1) then
      return
      
    ! otherwise, the current cell intersects with either face(s) or edge(s). Start patch searching.
    else
      edgeIdx = self % grid % getGridPhiCapital(cellIdxsMat)

      !!!!!
      !print*, "edgeIdx", edgeIdx
      !!!!!

      if (edgeIdx == 0) then
        ! push the coordinates away from the current vertex (= phi)
        ! (needs to be changed) (possible improvement/acceleration for the rest of the subroutine below?)
        vertexIdx = self % grid % getGridPhi(cellIdxsMat)
        phiCoord = vertices % getVertexCoordinates(vertexIdx)
        dummyVector = r - phiCoord
        r = phiCoord + (self % grid % getGridWStar())/(norm2(dummyVector))*(dummyVector)

        !!!!!
        !print*, "Phi", vertexIdx
        !print*, "pushed coord", r
        !!!!!

        ! find updated cartesian cell indices
        ! do i = 1, 3
        !   baseIntegerCoord(i) = floor((r(i) - gridBounds_min(i))*(gridSpacingReciprocal))
        ! end do
        baseIntegerCoord(:) = floor((r(:) - gridBounds_min(:))*(gridSpacingReciprocal))

        !!!!!
        ! print*, "baseCoordPUSHED", baseIntegerCoord
        ! print*, self % grid % getGridChi(baseIntegerCoord)
        ! print*, self % grid % getGridPhi(baseIntegerCoord)
        ! print*, self % grid % getGridPhiCapital(baseIntegerCoord)
        !!!!!

        ! get updated edge and element index
        !!!!!edgeIdx = self % grid % getGridPhiCapital(baseIntegerCoord)
        
        !!!!!
        !print*, edgeIdx
        !!!!!

        potentialElementIdx = self % grid % getGridChi(baseIntegerCoord, cellIdxsMat)

        ! if pushed coordinate has direct mapping for element index, use that
        ! (needs to be changed) (possible acceleration for this and other parts of the subroutine)
        ! (needs checking) (is pushed position has direct mapping for element idx, is it guaranteed to lie inside. OW, ">=" not "/=")
        if (potentialElementIdx /= 0) then
          call coords % setElementIdx(potentialElementIdx)
          ! (needs to be changed) (temp:there is no internal subdivision)
          !call coords % setParentElementIdx(elements % getElementParentIdx(potentialElementIdx))
          call coords % setParentElementIdx(potentialElementIdx)
          return
        end if

        edgeIdx = self % grid % getGridPhiCapital(cellIdxsMat)

      end if
    
      ! calculate pseudo angle
      currEdgeVertexIdxs = edges % getEdgeVertexIdxs(edgeIdx)
      v_eCoord = vertices % getVertexCoordinates(currEdgeVertexIdxs(2))
      rLocalCoord = r - v_eCoord
      xLocalCoord = dot_product(rLocalCoord, edges % getEdgeLocalBasis1(edgeIdx))
      yLocalCoord = dot_product(rLocalCoord, edges % getEdgeLocalBasis2(edgeIdx))
      thetaHat = SIGN(1 - (xLocalCoord / (abs(xLocalCoord) + abs(yLocalCoord))), yLocalCoord)

      ! perform binary search and return index pointer
      ! (needs checking) (index order and mechanics)
      !isBoundary = edges % getEdgeIsBoundary(edgeIdx) !!! not needed anymore
      pointerIdx = binarySearchAngle(edges % getEdgeAnglesArray(edgeIdx), thetaHat)
      elementIdxsArray = edges % getEdgeElementIdxsArray(edgeIdx)
      potentialElementIdx = elementIdxsArray(pointerIdx)

      ! if the neutron turns out to lie outside the mesh domain, return 
      if (potentialElementIdx == 0) return

      call coords % setElementIdx(potentialElementIdx)
      call coords % setParentElementIdx(elements % getElementParentIdx(potentialElementIdx))

    !!!!!
    ! print*, self % grid % getGridPhiCapital(baseIntegerCoord)
    ! print*, self % grid % getGridPhi(baseIntegerCoord)
    ! call fatalError("end","end")
    !!!!!

      return

    end if

  end subroutine findHostElement

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Non-bit-trick (saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! !!
  ! !!
  ! !!
  ! subroutine findHostElement(self, vertices, edges, faces, elements, coords)
  !   class(patchMultiAcceleration), intent(in)      :: self
  !   class(vertexShelf), intent(in)                 :: vertices
  !   class(edgeShelf), intent(in)                   :: edges
  !   type(faceShelf), intent(in)                    :: faces
  !   type(elementShelf), intent(in)                 :: elements
  !   type(coord), intent(inout)                     :: coords
  !   integer(shortInt)                              :: potentialElementIdx, i, edgeIdx, vertexIdx, pointerIdx
  !   integer(shortInt), dimension(3)                :: baseIntegerCoord
  !   real(defReal), dimension(3)                    :: r, v_eCoord, rLocalCoord, dummyVector, phiCoord
  !   integer(shortInt), dimension(2)                :: currEdgeVertexIdxs
  !   real(defReal)                                  :: thetaHat, xLocalCoord, yLocalCoord
  !   integer(shortInt), dimension(:), allocatable   :: elementIdxsArray
  !   integer(shortInt), dimension(:,:), allocatable :: cellIdxsMat

  !   ! retrieve the coordinates of neutron
  !   ! (needs to be changed) (needs checking) (is it correct to use "getPositionToNudge" or other coordinates?)
  !   r = coords % getPositionToNudge()
  
  !   !!!!!
  !   ! print*, "-----------------------------------------------"
  !   ! print*, "coord", r
  !   !!!!!

  !   ! check if the position of neutron is inside the catesian grid bounds
  !   if (self % grid % getGridIsOutsideBounds(r)) return

  !   ! retrieve element index from chi mapping.
  !   potentialElementIdx = self % grid % getGridChi(r, cellIdxsMat)

  !   !!!!!
  !   !print*, "indices", baseIntegerCoord
  !   !print*, "PotentialElementIdx", potentialElementIdx
  !   !!!!!
    
  !   ! if element index is valid (the current cell, characterised by "cellIdxs", is fully contained within that element)
  !   if (potentialElementIdx > 0) then
  !     call coords % setElementIdx(potentialElementIdx)

  !   !!!!!
  !   !print*, "indices", baseIntegerCoord
  !   !print*, "PotentialElementIdx", potentialElementIdx
  !   !!!!!
      
  !     call coords % setParentElementIdx(elements % getElementParentIdx(potentialElementIdx))
  !     return

  !   ! in case the current cell lies outside the computational domain for the unstructured mesh, return.
  !   ! this is tested after testing if (chi > 0) because that is the most likely case in terms of the number of the cells
  !   elseif (potentialElementIdx == -1) then
  !     return
      
  !   ! otherwise, the current cell intersects with either face(s) or edge(s). Start patch searching.
  !   else
  !     edgeIdx = self % grid % getGridPhiCapital(cellIdxsMat)

  !     !!!!!
  !     !print*, "edgeIdx", edgeIdx
  !     !!!!!

  !     if (edgeIdx == 0) then
  !       ! push the coordinates away from the current vertex (= phi)
  !       ! (needs to be changed) (possible improvement/acceleration for the rest of the subroutine below?)
  !       vertexIdx = self % grid % getGridPhi(cellIdxsMat)
  !       phiCoord = vertices % getVertexCoordinates(vertexIdx)
  !       dummyVector = r - phiCoord
  !       r = phiCoord + (self % grid % getGridWStar())/(norm2(dummyVector))*(dummyVector)

  !       !!!!!
  !       !print*, "Phi", vertexIdx
  !       !print*, "pushed coord", r
  !       !!!!!

  !       !!!!!
  !       ! print*, "baseCoordPUSHED", baseIntegerCoord
  !       ! print*, self % grid % getGridChi(baseIntegerCoord)
  !       ! print*, self % grid % getGridPhi(baseIntegerCoord)
  !       ! print*, self % grid % getGridPhiCapital(baseIntegerCoord)
  !       !!!!!

  !       ! get updated edge and element index
  !       !!!!!edgeIdx = self % grid % getGridPhiCapital(baseIntegerCoord)
        
  !       !!!!!
  !       !print*, edgeIdx
  !       !!!!!

  !       potentialElementIdx = self % grid % getGridChi(r, cellIdxsMat)

  !       ! if pushed coordinate has direct mapping for element index, use that
  !       ! (needs to be changed) (possible acceleration for this and other parts of the subroutine)
  !       ! (needs checking) (is pushed position has direct mapping for element idx, is it guaranteed to lie inside. OW, ">=" not "/=")
  !       if (potentialElementIdx /= 0) then
  !         call coords % setElementIdx(potentialElementIdx)
  !         ! (needs to be changed) (temp:there is no internal subdivision)
  !         !call coords % setParentElementIdx(elements % getElementParentIdx(potentialElementIdx))
  !         call coords % setParentElementIdx(potentialElementIdx)
  !         return
  !       end if

  !       edgeIdx = self % grid % getGridPhiCapital(cellIdxsMat)

  !     end if
    
  !     ! calculate pseudo angle
  !     currEdgeVertexIdxs = edges % getEdgeVertexIdxs(edgeIdx)
  !     v_eCoord = vertices % getVertexCoordinates(currEdgeVertexIdxs(2))
  !     rLocalCoord = r - v_eCoord
  !     xLocalCoord = dot_product(rLocalCoord, edges % getEdgeLocalBasis1(edgeIdx))
  !     yLocalCoord = dot_product(rLocalCoord, edges % getEdgeLocalBasis2(edgeIdx))
  !     thetaHat = SIGN(1 - (xLocalCoord / (abs(xLocalCoord) + abs(yLocalCoord))), yLocalCoord)

  !     ! perform binary search and return index pointer
  !     ! (needs checking) (index order and mechanics)
  !     !isBoundary = edges % getEdgeIsBoundary(edgeIdx) !!! not needed anymore
  !     pointerIdx = binarySearchAngle(edges % getEdgeAnglesArray(edgeIdx), thetaHat)
  !     elementIdxsArray = edges % getEdgeElementIdxsArray(edgeIdx)
  !     potentialElementIdx = elementIdxsArray(pointerIdx)

  !     ! if the neutron turns out to lie outside the mesh domain, return 
  !     if (potentialElementIdx == 0) return

  !     call coords % setElementIdx(potentialElementIdx)
  !     call coords % setParentElementIdx(elements % getElementParentIdx(potentialElementIdx))

  !   !!!!!
  !   ! print*, self % grid % getGridPhiCapital(baseIntegerCoord)
  !   ! print*, self % grid % getGridPhi(baseIntegerCoord)
  !   ! call fatalError("end","end")
  !   !!!!!

  !     return

  !   end if

  ! end subroutine findHostElement

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  !!
  !!
  !!
  subroutine init(self, vertices, edges, faces, elements)
    class(patchMultiAcceleration), intent(inout)     :: self
    type(vertexShelf), intent(in)                    :: vertices
    type(faceShelf), intent(inout)                   :: faces
    type(elementShelf), intent(in)                   :: elements
    type(edgeShelf), intent(inout)                   :: edges

    ! Simply initialise the cartesian multi-layered grid
    call self % grid % init(vertices, edges, faces, elements)

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(patchMultiAcceleration), intent(inout) :: self

  end subroutine kill

end module patchMultiAcceleration_class