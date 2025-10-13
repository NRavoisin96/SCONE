module cartesianCellIntermediate2_class
  
  use numPrecision
  use universalVariables,              only : ZERO
  use vertexShelf_class,               only : vertexShelf
  use edgeShelf_class,                 only : edgeShelf
  use faceShelf_class,                 only : faceShelf
  use cartesianInitProcedures
  use cartesianGenericProcedures,      only : intBinarySearch, remove_elements  
  use cartesianGridSubLayer_inter,     only : cartesianGridSubLayer
  use cartesianGridFinest_class,       only : cartesianGridFinest
  use genericProcedures,               only : append, fatalError
  use ragged3dMatrix_class,            only : ragged3d
  use dynamic2dMatSet_class,           only : dynamic2dMatSet

  implicit none
  private
  
  !!
  !!
  type, public                                          :: cartesianCellIntermediate2
    private
    class(cartesianGridSubLayer), pointer               :: subGrid => null()
    integer(shortInt)                                   :: chi = 0
    ! (needs to be changed) (cell centre should be a property to avoid repetative calc.)
    ! (due to limited memory, this is calculated each time needed)
    ! (try to avoid adding properties tho due to memory)

  contains

    ! (needs to be changed) (If we mark cells who lie outside of mesh in upper layers,)
    ! (we can say a huge memory and a bit of initialisation time and in-cycle time.)
    ! (see OneNote for more detail: Ideas>To do)

    ! Build procedures.
    procedure                                    :: cellTestPolyhedronInclusion
    procedure                                    :: refineCell
    procedure                                    :: refineCell2
    procedure                                    :: constructNRefineCell
    ! Runtime procedures.
    procedure                                    :: getChi
    procedure                                    :: getPhi
    procedure                                    :: getPhiCapital
    ! Analysis procedures
    procedure                                    :: cellGetNumberOfCells
  end type cartesianCellIntermediate2

contains

  !!
  !!
  !!
  subroutine cellTestPolyhedronInclusion(self, faces, currElementFaceIdxs, centroid, &
                                              faceNormalSigns, elementIdx)
    class(cartesianCellIntermediate2), intent(inout)     :: self
    class(faceShelf), intent(in)                        :: faces
    integer(shortInt), dimension(:), intent(in)         :: currElementFaceIdxs
    real(defReal), dimension(3), intent(in)             :: centroid
    real(defReal), dimension(:,:), intent(in)           :: faceNormalSigns
    integer(shortInt), intent(in)                       :: elementIdx

    ! if current cell is found to be entirely contained within a polyhedron, there cannot not be other polyhedra
    if (self % chi /= 0) return

    ! Otherwise, continue testing and appending the array of candidate element indices 
    call testPolyhedronInclusion(faces, currElementFaceIdxs, centroid, faceNormalSigns, elementIdx, self % chi)

  end subroutine cellTestPolyhedronInclusion

  !!
  !!
  !!
  subroutine refineCell(self, vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                             currLayer, candidateElementIdxs, newGridBoundsMin, alpha, wStar)
    class(cartesianCellIntermediate2), intent(inout)     :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(inout)                     :: edges
    class(faceShelf), intent(inout)                     :: faces
    class(elementShelf), intent(in)                     :: elements
    real(defReal), dimension(:), intent(in)             :: spacing, spacingInv
    integer(shortInt), dimension(:,:), intent(in)       :: n_xyz
    integer(shortInt), intent(in)                       :: n_layers, currLayer
    integer(shortInt), dimension(:), intent(in)         :: candidateElementIdxs
    real(defReal), dimension(3), intent(in)             :: newGridBoundsMin
    real(defReal), intent(in)                           :: alpha, wStar

    !if the current cell is not entirely contained within a polyhedron, then refine the grid
    if (self % chi == 0) then

      ! (needs to be changed) (pointers cannot point to the same class)
      ! if (n_layers > currLayer+1) then 
      !   allocate(cartesianGridIntermediate3:: self % subGrid)
      ! else
        allocate(cartesianGridFinest:: self % subGrid)
      ! end if

      call self % subgrid % init(vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                                 currLayer + 1, candidateElementIdxs, newGridBoundsMin, alpha, wStar)

    end if 

  end subroutine refineCell

  !!
  !!
  !!
  subroutine refineCell2(self, vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                    currLayer, candidateElementIdxsOld, newGridBoundsMin, alpha, wStar, normalSignsMatOld)
    class(cartesianCellIntermediate2), intent(inout)     :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(inout)                     :: edges
    class(faceShelf), intent(inout)                     :: faces
    class(elementShelf), intent(in)                     :: elements
    real(defReal), dimension(:), intent(in)             :: spacing, spacingInv
    integer(shortInt), dimension(:,:), intent(in)       :: n_xyz
    integer(shortInt), intent(in)                       :: n_layers, currLayer
    integer(shortInt), dimension(:), intent(in)         :: candidateElementIdxsOld
    real(defReal), dimension(3), intent(in)             :: newGridBoundsMin
    real(defReal), intent(in)                           :: alpha, wStar
    type(dynamic2dMatSet), intent(in)                   :: normalSignsMatOld
    type(dynamic2dMatSet)                               :: normalSignsMat
    real(defReal), dimension(:,:), allocatable          :: faceNormalSigns
    integer(shortInt), dimension(:), allocatable        :: currElementFaceIdxs, removedFaceIdxsInArr, &
                                                           candidateElementIdxs, removedElementIdxsInArr
    integer(shortInt)                                   :: i
    real(defReal), dimension(3)                         :: centroid
    logical(defBool)                                    :: isOut

    normalSignsMat = normalSignsMatOld
    candidateElementIdxs = candidateElementIdxsOld
    if (allocated(removedElementIdxsInArr)) deallocate(removedElementIdxsInArr)

    ! Calculate centroid from newGridBoundsMin
    do i = 1, 3
      centroid(i) = newGridBoundsMin(i) + spacing(currLayer)*0.5
    end do

    allocate(faceNormalSigns(3, 2))
    allocate(currElementFaceIdxs(1))

    do i = 1, normalSignsMat % nslices()

      deallocate(faceNorMalSigns)
      deallocate(currElementFaceIdxs)
      if (allocated(removedFaceIdxsInArr)) deallocate(removedFaceIdxsInArr)

      call normalSignsMat % get_copy(i, faceNormalSigns, currElementFaceIdxs)

      call testPolyhedronInclusion3(faces, currElementFaceIdxs, centroid, faceNormalSigns, &
                            candidateElementIdxs(i), self % chi, removedFaceIdxsInArr, isOut)

      if (self % chi > 0) return

      if (isOut) then
        call append(removedElementIdxsInArr, i)
      else
        if (allocated(removedFaceIdxsInArr)) then
          call normalSignsMat % delete_columns(i, removedFaceIdxsInArr)
          !print*, "Yes"
          ! print*, size(removedFaceIdxsInArr)
          ! print*, "@@@", size(currElementFaceIdxs)
        else
          !print*, "Not"
        end if
      end if




    end do

    !@@ Update NormalSignsMat (remove and scale)
    if (allocated(removedElementIdxsInArr)) then
      call remove_elements(candidateElementIdxs, removedElementIdxsInArr)
      call normalSignsMat % deleteMany(removedElementIdxsInArr)
    end if
    call normalSignsMat % scale(spacingInv(currLayer)*spacing(currLayer+1))



    !if the current cell is not entirely contained within a polyhedron, then refine the grid
    if (self % chi == 0) then !@@ can be removed

      ! (needs to be changed) (pointers cannot point to the same class)
      ! if (n_layers > currLayer+1) then 
      !   allocate(cartesianGridIntermediate2:: self % subGrid)
      ! else
        allocate(cartesianGridFinest:: self % subGrid)
      ! end if

      call self % subgrid % init2(vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                        currLayer + 1, candidateElementIdxs, newGridBoundsMin, alpha, wStar, normalSignsMat)

    end if 

    call normalSignsMat % kill()

  end subroutine refineCell2

  !!
  !!
  !!
  subroutine constructNRefineCell(self, vertices, edges, faces, elements, spacing, spacingInv, &
                   n_xyz, n_layers, currLayer, intersectedFaceIdxsOld, newGridBoundsMin, alpha, wStar, &
                   extraDistanceArrOld, candidateElementIdxs, normalSignsMatOld, &
                   circumscribedBallRadius, targetDistance, targetDistanceSqr)
    class(cartesianCellIntermediate2), intent(inout)    :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(inout)                     :: edges
    class(faceShelf), intent(inout)                     :: faces
    class(elementShelf), intent(in)                     :: elements
    real(defReal), dimension(:), intent(in)             :: spacing, spacingInv
    integer(shortInt), dimension(:,:), intent(in)       :: n_xyz
    integer(shortInt), intent(in)                       :: n_layers, currLayer
    integer(shortInt), dimension(:), intent(in)         :: intersectedFaceIdxsOld, candidateElementIdxs
    real(defReal), dimension(3), intent(in)             :: newGridBoundsMin
    real(defReal), intent(in)                           :: alpha, wStar, circumscribedBallRadius, targetDistance, &
                                                           targetDistanceSqr
    real(defReal), dimension(:), intent(in)             :: extraDistanceArrOld
    type(ragged3d), intent(in)                          :: normalSignsMatOld
    type(ragged3d)                                      :: normalSignsMat
    real(defReal), dimension(3)                         :: centroid, currFaceNormal
    integer(shortInt)                                   :: i, currFaceIdx, currElementIdxInArr
    integer(shortInt), dimension(:), allocatable        :: currVertexIdxs, currFaceEdgeIdxs, &
                                                           intersectedFaceIdxs, removedFaceIdxs, &
                                                           currCandidateElementIdxs, testElementIdxs, &
                                                           currElementFaceIdxs 
    real(defReal)                                       :: extraDistance, currFaceConst
    real(defReal), dimension(:), allocatable            :: extraDistanceArr
    real(defReal), dimension(:,:), allocatable          :: faceNormalSigns

    ! Save a copy of intersectedFaceIdxs and extraDistanceArr because the array will be modified and passed to the sub-layer.
    if (allocated(intersectedFaceIdxs)) deallocate(intersectedFaceIdxs)
    if (allocated(extraDistanceArr)) deallocate(extraDistanceArr)
    call normalSignsMat % kill()
    intersectedFaceIdxs = intersectedFaceIdxsOld
    extraDistanceArr = extraDistanceArrOld
    normalSignsMat = normalSignsMatOld

    ! Calculate centroid from newGridBoundsMin
    do i = 1, 3
      centroid(i) = newGridBoundsMin(i) + spacing(currLayer)*0.5
    end do

    ! Face intersection tests agains all faces in the array "intersectedFaceIdxs"
    ! Looping index decreases by one each iteration because currFaceIdx can be removed.
    do i = size(intersectedFaceIdxs), 1, -1
      currFaceIdx = intersectedFaceIdxs(i)
      extraDistance = extraDistanceArr(i)

      ! Retrieve face-specific parameters
      if (allocated(currVertexIdxs)) deallocate(currVertexIdxs)
      if (allocated(currFaceEdgeIdxs)) deallocate(currFaceEdgeIdxs)
      if (allocated(removedFaceIdxs)) deallocate(removedFaceIdxs)
      currVertexIdxs = faces % getFaceVertexIdxs(currFaceIdx)
      currFaceNormal = faces % getFaceNormal(currFaceIdx)
      currFaceEdgeIdxs = faces % getFaceEdgeIdxs(currFaceIdx)

      call testFaceIntersectionNonCoarsest(vertices, edges, faces, currVertexIdxs, extraDistance, &
                    currFaceNormal, centroid, spacing(currLayer), currFaceIdx, currFaceEdgeIdxs, &
                    intersectedFaceIdxs, i, extraDistanceArr, removedFaceIdxs)

    end do

    ! If there turns out to be no face intersecting with the cell, perform polyhedron inclusion test
    if (.NOT. allocated(intersectedFaceIdxs)) then

      ! If there is a single face removed from the array "intersectedFaceIdxs", perform special (simplified) inclusion test
      if (size(removedFaceIdxs) == 1) then
        if (allocated(testElementIdxs)) deallocate(testElementIdxs)
        currFaceNormal = faces % getFaceNormal(removedFaceIdxs(1))
        currFaceConst = faces % getFaceConst(removedFaceIdxs(1))
        testElementIdxs = faces % getFaceElementIdxs(removedFaceIdxs(1))

        ! Test if the particle lies in the owner element of the face.
        ! If true, the particle lies in the non-owner element of the face
        ! Currently, the element indices of a give face is: [owner element, non-owner element]
        if (faceHalfSpaceTest(currFaceNormal, centroid, currFaceConst)) then

          ! If it is a bondary face, the only element attached is the owner of the face.
          ! Hence, if faceHalfSpaceTest tells that the neutron lies outside of the owner element, then we know its outside of the mesh.
          ! (There was no problem in this logic, but might not work with non-OpenFoam mesh data format)
          if (faces % getFaceIsBoundary(removedFaceIdxs(1))) then
            self % chi = -1
            return
          else
            self % chi = testElementIdxs(2)
            return
          end if

        else

          self % chi = testElementIdxs(1)
          return

        end if
        
      ! If there is more than one removed face indices
      else

        if (allocated(currCandidateElementIdxs)) deallocate(currCandidateElementIdxs)
        currCandidateElementIdxs = faces % getFaceElementIdxs(removedFaceIdxs)

        ! Loop over all candidate elements within which the neutron can possibly lie
        do i = 1, size(currCandidateElementIdxs)

          ! Consturct an array holding face indices of the current element
          if (allocated(currElementFaceIdxs)) deallocate(currElementFaceIdxs)
          currElementFaceIdxs = elements % getElementFaceIdxs(currCandidateElementIdxs(i))

          ! Construct a matrix holding faceNormal signs by retrieving pre-calculated data.
          if (allocated(faceNormalSigns)) deallocate(faceNormalSigns)
          currElementIdxInArr = intBinarySearch(candidateElementIdxs, currCandidateElementIdxs(i))
          faceNormalSigns = normalSignsMat % get(currElementIdxInArr)

          ! Perform polyhedron inclusion test.
          if (testPolyhedronInclusionnew(faces, currElementFaceIdxs, centroid, faceNormalSigns)) then
            ! If included, update chi mapping and return.
            self % chi = currCandidateElementIdxs(i)
            return
          end if

        end do

        ! If survived to this point, the particle does not lie within any single polyhedron
        ! nor the cell intersects with any faces. Hence, the particle lies outside of the mesh domain
        self % chi = -1
        return

      end if

    else 

      ! If there is still faces intersecting the cell, then we refine further
      ! (needs to be changed) (pointers cannot point to the same class)
      ! if (n_layers > currLayer+1) then 
      !   allocate(cartesianGridIntermediate2:: self % subGrid)
      ! else
        allocate(cartesianGridFinest:: self % subGrid)
      ! end if

      call normalSignsMat % scale(spacingInv(currLayer)*spacing(currLayer+1))
      extraDistanceArr(:) = extraDistanceArr(:)*(spacingInv(currLayer)*spacing(currLayer+1))

      call self % subgrid % initt(vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                                  currLayer+1, intersectedFaceIdxs, newGridBoundsMin, alpha, wStar, &
                                  extraDistanceArr, candidateElementIdxs, normalSignsMat, &
                                  circumscribedBallRadius, targetDistance, targetDistanceSqr)

    end if


  end subroutine constructNRefineCell

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! bit-trick (not saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! !!
  ! !!
  ! !!
  ! function getChi(self, baseIntegerCoord, shift, mask, currLayer) result(chi)
  !   class(cartesianCellIntermediate2), intent(in)        :: self
  !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
  !   integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: chi

  !   if (self % chi /= 0) then
  !     chi = self % chi
  !   else
  !     chi = self % subGrid % getGridChi(baseIntegerCoord, shift, mask, currLayer+1)
  !   end if

  ! end function getChi

  ! !!
  ! !!
  ! !!
  ! function getPhi(self, baseIntegerCoord, shift, mask, currLayer) result(phi)
  !   class(cartesianCellIntermediate2), intent(in)        :: self
  !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
  !   integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: phi

  !   phi = self % subGrid % getGridPhi(baseIntegerCoord, shift, mask, currLayer+1)

  ! end function getPhi

  ! !!
  ! !!
  ! !!
  ! function getPhiCapital(self, baseIntegerCoord, shift, mask, currLayer) result(phiCapital)
  !   class(cartesianCellIntermediate2), intent(in)        :: self
  !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
  !   integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: phiCapital

  !   phiCapital = self % subGrid % getGridPhiCapital(baseIntegerCoord, shift, mask, currLayer+1)

  ! end function getPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! bit-trick (saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !!
  !!
  !!
  function getChi(self, baseIntegerCoord, shift, mask, currLayer, cellIdxsMat) result(chi)
    class(cartesianCellIntermediate2), intent(in)        :: self
    integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
    integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
    integer(shortInt), dimension(:,:), intent(inout)    :: cellIdxsMat
    integer(shortInt), intent(in)                       :: currLayer
    integer(shortInt)                                   :: chi

    if (self % chi /= 0) then
      chi = self % chi
    else
      chi = self % subGrid % getGridChi(baseIntegerCoord, shift, mask, currLayer+1, cellIdxsMat)
    end if

  end function getChi

  !!
  !!
  !!
  function getPhi(self, cellIdxsMat, currLayer) result(phi)
    class(cartesianCellIntermediate2), intent(in)        :: self
    integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
    integer(shortInt), intent(in)                       :: currLayer
    integer(shortInt)                                   :: phi

    phi = self % subGrid % getGridPhi(cellIdxsMat, currLayer+1)

  end function getPhi

  !!
  !!
  !!
  function getPhiCapital(self, cellIdxsMat, currLayer) result(phiCapital)
    class(cartesianCellIntermediate2), intent(in)        :: self
    integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
    integer(shortInt), intent(in)                       :: currLayer
    integer(shortInt)                                   :: phiCapital

    phiCapital = self % subGrid % getGridPhiCapital(cellIdxsMat, currLayer+1)

  end function getPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Non bit-trick (saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! !!
  ! !!
  ! !!
  ! function getChi(self, r, gridBounds_min, spacingInv, currLayer, cellIdxsMat, nSub_xyz) result(chi)
  !   class(cartesianCellIntermediate2), intent(in)        :: self
  !   real(defReal), dimension(3), intent(in)             :: r, gridBounds_min                                    
  !   real(defReal), dimension(:), intent(in)             :: spacingInv
  !   integer(shortInt), dimension(:,:), intent(inout)    :: cellIdxsMat
  !   integer(shortInt), dimension(:,:), intent(in)       :: nSub_xyz
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: chi

  !   if (self % chi /= 0) then
  !     chi = self % chi
  !   else
  !     chi = self % subGrid % getGridChi(r, gridBounds_min, spacingInv, currLayer+1, cellIdxsMat, nSub_xyz)
  !   end if

  ! end function getChi

  ! !!
  ! !!
  ! !!
  ! function getPhi(self, cellIdxsMat, currLayer) result(phi)
  !   class(cartesianCellIntermediate2), intent(in)        :: self
  !   integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: phi

  !   phi = self % subGrid % getGridPhi(cellIdxsMat, currLayer+1)

  ! end function getPhi

  ! !!
  ! !!
  ! !!
  ! function getPhiCapital(self, cellIdxsMat, currLayer) result(phiCapital)
  !   class(cartesianCellIntermediate2), intent(in)        :: self
  !   integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: phiCapital

  !   phiCapital = self % subGrid % getGridPhiCapital(cellIdxsMat, currLayer+1)

  ! end function getPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  !!
  !!
  !!
  function cellGetNumberOfCells(self, localNxyz, n_layers, currLayer) result(output)
    class(cartesianCellIntermediate2), intent(in)        :: self
    integer(shortInt), dimension(:,:), intent(in)       :: localNxyz
    integer(shortInt), intent(in)                       :: n_layers, currLayer
    integer(shortInt), dimension(:), allocatable        :: output

    ! Allocate and initialise output array
    allocate(output(n_layers+1))
    output = 0

    if (self % chi /= 0) then
      output(currLayer) = 1
    else
      output = self % subgrid % getNumberOfCells(localNxyz, n_layers, currLayer + 1)
    end if
    
  end function cellGetNumberOfCells


end module cartesianCellIntermediate2_class