module cartesianCellFinest_class
  
  use numPrecision
  use universalVariables,           only : ZERO
  use vertexShelf_class,            only : vertexShelf
  use edgeShelf_class,              only : edgeShelf
  use faceShelf_class,              only : faceShelf
  use cartesianGenericProcedures,   only : sortByHighestFrequency
  use cartesianInitProcedures
  use cartesianGenericProcedures,   only : intBinarySearch
  use ragged3dMatrix_class,         only : ragged3d
  use genericProcedures,            only : append, fatalError

  implicit none
  private
  
  !!
  !!
  type, public                                          :: cartesianCellFinest
    private
    integer(shortInt), allocatable                      :: phi, phiCapital, chi !###
    integer(shortInt), dimension(:), allocatable        :: faceIdxs !###
    ! (needs to be changed) (cell centre should be a property to avoid repetative calc.)
    ! (due to limited memory, this is calculated each time needed)
    ! (try to avoid adding properties tho due to memory)

  contains
    ! Build procedures.
    procedure                                    :: cellTestEdgeIntersection
    procedure                                    :: cellTestPolyhedronInclusion
    procedure                                    :: cellTestPolyhedronInclusion2
    procedure                                    :: cellTestFaceIntersection
    procedure                                    :: cellTestFaceIntersection2
    procedure                                    :: cellConstructMapSingleFace
    procedure                                    :: setIsOutsideMesh
    procedure                                    :: cellFinitePrecision
    procedure                                    :: cellClearRedundancy
    procedure                                    :: cellAllocateAttributes
    procedure                                    :: constructNRefineCell
    ! Runtime procedures.
    procedure                                    :: getPhiCapital
    procedure                                    :: getPhi
    procedure                                    :: getChi
    ! Analysis procedures
    procedure                                    :: cellGetNumberOfCells

  end type cartesianCellFinest

contains

  !!
  !!
  !! (needs to be changed) (Inefficiency due to refactoring original code) (mapping matrices in cartesianCell_class)
  !! (can be moved and stored as a attribute in cartesianGrid. Then call testEdgeIntersection directly from cartesianGrid_class)
  !! (For now this is fine becase only 4% of initialisation time is increased by this - in face no increase in init time observed later.)
  !! (Tried this before and after specifying "elemental" and "pure". Keeping the original architecture (with separte class for cells))
  !! (is faster than without for initialisation, and negligible difference for in-cycle)
  subroutine cellTestEdgeIntersection(self, vertices, edges, edgeIdx, circumscribedBallRadius, &
                                      targetDistance, centroid, currEdgeVector, currVertexIdxs, a)
    class(cartesianCellFinest), intent(inout)           :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(in)                        :: edges
    integer(shortInt), intent(in)                       :: edgeIdx
    real(defReal), intent(in)                           :: circumscribedBallRadius, targetDistance, a
    real(defReal), dimension(3), intent(in)             :: centroid, currEdgeVector
    integer(shortInt), dimension(2), intent(in)         :: currVertexIdxs

          ! if (edgeIdx == 395) then
          !   print*, "circumscribedBallRadius", circumscribedBallRadius
          !   print*, "targetDistance", targetDistance
          !   print*, "EdgeVector", edges % getEdgeVector(edgeIdx)
          !   print*, "EdgeVertices", edges % getEdgeVertexIdxs(edgeIdx)
          !   print*, "the dot product", edges % getEdgeDotProductOfVector(edgeIdx)

          !   call fatalError("done", "done")


          ! end if

    call testEdgeIntersection(vertices, edges, edgeIdx, circumscribedBallRadius, targetDistance, centroid, &
                              currEdgeVector, currVertexIdxs, a, self % phi, self % phiCapital)

  end subroutine cellTestEdgeIntersection

  !!
  !!
  !!
  subroutine cellTestPolyhedronInclusion(self, faces, currElementFaceIdxs, centroid, &
                                         faceNormalSigns, elementIdx)
    class(cartesianCellFinest), intent(inout)           :: self
    class(faceShelf), intent(in)                        :: faces
    integer(shortInt), dimension(:), intent(in)         :: currElementFaceIdxs
    real(defReal), dimension(3), intent(in)             :: centroid
    real(defReal), dimension(:,:), intent(in)           :: faceNormalSigns
    integer(shortInt), intent(in)                       :: elementIdx

    ! if (elementIdx == 100) then
    !   print*, "currElementFaceIdxs", currElementFaceIdxs
    !   print*, "faceNormalSigns", faceNormalSigns
    !   call fatalError("done", "done")
    ! end if

    call testPolyhedronInclusion(faces, currElementFaceIdxs, centroid, faceNormalSigns, elementIdx, self % chi)

  end subroutine cellTestPolyhedronInclusion

  !!
  !!
  !!
  subroutine cellTestPolyhedronInclusion2(self, faces, currElementFaceIdxs, centroid, &
                                         elementIdx, removedFaceIdxsInArr, isOut, currLayer)
    class(cartesianCellFinest), intent(inout)                    :: self
    class(faceShelf), intent(in)                                 :: faces
    integer(shortInt), dimension(:), intent(in)                  :: currElementFaceIdxs
    real(defReal), dimension(3), intent(in)                      :: centroid
    integer(shortInt), intent(in)                                :: elementIdx, currLayer
    integer(shortInt), dimension(:), allocatable, intent(out)    :: removedFaceIdxsInArr
    logical(defBool), intent(out)                                :: isOut     


    ! if (elementIdx == 100) then
    !   print*, "currElementFaceIdxs", currElementFaceIdxs
    !   print*, "faceNormalSigns", faceNormalSigns
    !   call fatalError("done", "done")
    ! end if

    call testPolyhedronInclusion2New(faces, currElementFaceIdxs, centroid, elementIdx, self % chi, &
                                 removedFaceIdxsInArr, isOut, currLayer)

  end subroutine cellTestPolyhedronInclusion2

  !!
  !!
  !!
  subroutine cellTestFaceIntersection(self, vertices, edges, faces, currVertexIdxs, extraDistance, currFaceNormal, &
                                      centroid, cellSpacing, faceIdx, currFaceEdgeIdxs, targetDistance)
    class(cartesianCellFinest), intent(inout)           :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(in)                        :: edges
    class(faceShelf), intent(in)                        :: faces
    integer(shortInt), dimension(:), intent(in)         :: currVertexIdxs, currFaceEdgeIdxs
    real(defReal), dimension(3), intent(in)             :: currFaceNormal, centroid
    real(defReal), intent(in)                           :: extraDistance, cellSpacing, targetDistance
    integer(shortInt), intent(in)                       :: faceIdx

    ! if (faceIdx == 500) then
    !   print*, "faceIdx=", faceIdx
    !   print*, "currVertexIdxs", currVertexIdxs
    !   print*, "extraDistance", extraDistance
    !   print*, "currFaceNormal", currFaceNormal
    !   print*, "cellSpacing", cellSpacing
    !   print*, "currFaceEdgeIdxs", currFaceEdgeIdxs
    !   print*, "targetDistance", targetDistance

    !   call fatalError("done", "done")
    ! end if
    call testFaceIntersection(vertices, edges, faces, currVertexIdxs, extraDistance, currFaceNormal, &
                                  centroid, cellSpacing, faceIdx, currFaceEdgeIdxs, targetDistance, self % chi, &
                                  self % phi, self % phiCapital, self % faceIdxs)

  end subroutine cellTestFaceIntersection

  !!
  !!
  !!
  subroutine cellTestFaceIntersection2(self, vertices, edges, faces, currVertexIdxs, currFaceNormal, &
                                      centroid, cellSpacing, faceIdx, currFaceEdgeIdxs, targetDistance, currLayer)
    class(cartesianCellFinest), intent(inout)           :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(in)                        :: edges
    class(faceShelf), intent(in)                        :: faces
    integer(shortInt), dimension(:), intent(in)         :: currVertexIdxs, currFaceEdgeIdxs
    real(defReal), dimension(3), intent(in)             :: currFaceNormal, centroid
    real(defReal), intent(in)                           :: cellSpacing, targetDistance
    integer(shortInt), intent(in)                       :: faceIdx, currLayer

    ! if (faceIdx == 500) then
    !   print*, "faceIdx=", faceIdx
    !   print*, "currVertexIdxs", currVertexIdxs
    !   print*, "extraDistance", extraDistance
    !   print*, "currFaceNormal", currFaceNormal
    !   print*, "cellSpacing", cellSpacing
    !   print*, "currFaceEdgeIdxs", currFaceEdgeIdxs
    !   print*, "targetDistance", targetDistance

    !   call fatalError("done", "done")
    ! end if

    call testFaceIntersection2(vertices, edges, faces, currVertexIdxs, currFaceNormal, &
                                  centroid, cellSpacing, faceIdx, currFaceEdgeIdxs, targetDistance, self % chi, &
                                  self % phi, self % phiCapital, self % faceIdxs, currLayer)

  end subroutine cellTestFaceIntersection2

  !!
  !!
  !!
  subroutine cellConstructMapSingleFace(self, faces)
    class(cartesianCellFinest), intent(inout)           :: self
    class(faceShelf), intent(in)                        :: faces
    integer(shortInt), dimension(:), allocatable        :: currFaceEdgeIdxs!, currFaceVertexIdxs
    integer(shortInt)     :: temp

    ! call constructMapSingleFace(faces, self % faceIdxs, self % phiCapital)
    if (self % faceIdxs(1) /= 0 .AND. self % faceIdxs(2) == 0) then
      self % chi = -abs(self % faceIdxs(1))
    end if


  end subroutine cellConstructMapSingleFace

  !!
  !!
  !!
  subroutine setIsOutsideMesh(self, no) !"""
    class(cartesianCellFinest), intent(inout)           :: self
    integer(shortInt), intent(in)                       :: no

    ! if a cell intersects with neither any edge nor face, and it is not contained in a single polyhedron,
    ! then, this cell lies outside the computational domain for the unstructured mesh
    ! (needs to be changed) (construct nested if to allow early exit)
    if (self % chi == 0 .AND. self % phi == 0 .AND. self % phiCapital == 0) then
      ! if lies outside, set the chi value of the cell equal to -1. There are subroutines that test 
      ! if (chi != 0), but since this subroutine is called after all of those, they are unafftected.
      self % chi = no
    end if

  end subroutine setIsOutsideMesh

  !!
  !!
  !!
  subroutine cellFinitePrecision(self, faces, elements, centroid, elementIdx)
    class(cartesianCellFinest), intent(inout)           :: self
    class(faceShelf), intent(in)                        :: faces
    class(elementShelf), intent(in)                     :: elements
    real(defReal), dimension(3), intent(in)             :: centroid
    integer(shortInt), intent(in)                       :: elementIdx

    ! If the current Cartesian cell does not intersect with any faces nor included in a single mesh element,
    ! Test if centroid lies inside any mesh element. If yes, finite precision error messed it up. Hence, 
    ! update chi mapping info to that mesh element. If not, this cell lies within another element or is outside mesh domain.
    if (self%chi == 0) then                   ! The biggest proportion of cells will be included in a single mesh element
      !if (self%faceIdxs(1)==0) then           ! Then we test face intersection
      if (self%phi == 0) then
        if (self%phiCapital == 0) then
          call coverFinitePrecision(faces, elements, centroid, elementIdx, self%chi)  
        end if
      end if
    end if

  end subroutine cellFinitePrecision

  !!
  !!
  !!
  subroutine cellAllocateAttributes(self)
    class(cartesianCellFinest), intent(inout)           :: self

    !###
    allocate(self % faceIdxs, source=[0, 0])
    allocate(self % chi, source = 0)
    allocate(self % phi, source = 0)
    allocate(self %phiCapital, source = 0)
    !###

  end subroutine cellAllocateAttributes

  !!
  !!
  !!
  subroutine cellClearRedundancy(self)
    class(cartesianCellFinest), intent(inout)           :: self

    !###
    deallocate(self % faceIdxs)

    if (self % chi /= 0) then
      deallocate(self % phi)
      deallocate(self % phiCapital)
    else
      deallocate(self % chi)
      if (self % phi == 0) deallocate(self % phi)
      if (self % phiCapital == 0) deallocate(self % phiCapital)
    end if
    !###

  end subroutine cellClearRedundancy

  !!
  !!
  !!
  subroutine constructNRefineCell(self, vertices, edges, faces, elements, spacing, spacingInv, &
                   n_xyz, n_layers, currLayer, intersectedFaceIdxsOld, newGridBoundsMin, alpha, wStar, &
                   extraDistanceArrOld, candidateElementIdxs, normalSignsMatOld, &
                   circumscribedBallRadius, targetDistance, targetDistanceSqr)
    class(cartesianCellFinest), intent(inout)           :: self
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
    integer(shortInt)                                   :: i, j, currFaceIdx, currElementIdxInArr
    integer(shortInt), dimension(:), allocatable        :: currVertexIdxs, currFaceEdgeIdxs, &
                                                           intersectedFaceIdxs, removedFaceIdxs, &
                                                           currCandidateElementIdxs, testElementIdxs, &
                                                           currElementFaceIdxs,candidateEdgeIdxs
    real(defReal)                                       :: extraDistance, currFaceConst
    real(defReal), dimension(:), allocatable            :: extraDistanceArr
    real(defReal), dimension(:,:), allocatable          :: faceNormalSigns
    integer(shortInt), dimension(2)                     :: faceIdxsToBeTested

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

! if (currFaceIdx == 500) then
!   print*, "faceIdx=", currFaceIdx
!   print*, "currVertexIdxs", currVertexIdxs
!   print*, "extraDistance", extraDistance
!   print*, "currFaceNormal", currFaceNormal
!   print*, "cellSpacing", spacing(currLayer)
!   print*, "currFaceEdgeIdxs", currFaceEdgeIdxs
!   print*, "targetDistance", targetDistance
!   print*, "extraDistanceArr", extraDistanceArr

!   call fatalError("done", "done")
! end if

      call testFaceIntersectionNonCoarsest(vertices, edges, faces, currVertexIdxs, extraDistance, &
                    currFaceNormal, centroid, spacing(currLayer), currFaceIdx, currFaceEdgeIdxs, &
                    intersectedFaceIdxs, i, extraDistanceArr, removedFaceIdxs)

    end do

    ! print*, "intersectedFaceIdxs", intersectedFaceIdxs

    ! If there turns out to be no face intersecting with the cell, perform polyhedron inclusion test
    if (.NOT. allocated(intersectedFaceIdxs)) then

      ! print*, "removedFaceIdxs", removedFaceIdxs

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


          ! if (currCandidateElementIdxs(i) == 100) then
          !   print*, "currElementFaceIdxs", currElementFaceIdxs
          !   print*, "faceNormalSigns", faceNormalSigns
          !   call fatalError("done", "done")

          ! end if
          ! print*, currCandidateElementIdxs, "@@@@@@@@@@@@@@@@"

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

    ! If there still is faces intersecting the cell, then we refine further  
    else 

      ! If there is a single face remaining to intersect the cell, simply take one of the edges of 
      ! the face to phi mapping.
      if (size(intersectedFaceIdxs) == 1) then
        if (allocated(currFaceEdgeIdxs)) deallocate(currFaceEdgeIdxs)
        currFaceEdgeIdxs = faces % getFaceEdgeIdxs(intersectedFaceIdxs(1))
        self % phiCapital = currFaceEdgeIdxs(1)

      else

        ! If there is more than one intersected faces, start by testing edge intersection
        
        ! Retrieve all relevant edge indices (duplicates are allowed)
        if (allocated(candidateEdgeIdxs)) deallocate(candidateEdgeIdxs)
        do i = 1, size(intersectedFaceIdxs)
          call append(candidateEdgeIdxs, faces % getFaceEdgeIdxs(intersectedFaceIdxs(i)))
        end do

        ! Then, construct a sorted candidateEdgeIdxs based on their frequencies.
        ! More frequent edge index is more likely to be intersected when multiple faces are intersected.
        candidateEdgeIdxs = sortByHighestFrequency(candidateEdgeIdxs)

        ! Perform edge intersection tests for edges in candidateEdgeIdxs
        do i = 1, size(candidateEdgeIdxs)

          ! if (candidateEdgeIdxs(i) == 395) then
          !   print*, "circumscribedBallRadius", circumscribedBallRadius
          !   print*, "targetDistance", targetDistance
          !   print*, "EdgeVector", edges % getEdgeVector(candidateEdgeIdxs(i))
          !   print*, "EdgeVertices", edges % getEdgeVertexIdxs(candidateEdgeIdxs(i))
          !   print*, "the dot product", edges % getEdgeDotProductOfVector(candidateEdgeIdxs(i))

          !   call fatalError("done", "done")


          ! end if




          call testEdgeIntersection(vertices, edges, candidateEdgeIdxs(i), circumscribedBallRadius, &
                                    targetDistance, centroid, edges % getEdgeVector(candidateEdgeIdxs(i)), &
                                    edges % getEdgeVertexIdxs(candidateEdgeIdxs(i)), &
                                    edges % getEdgeDotProductOfVector(candidateEdgeIdxs(i)), &
                                    self % phi, self % phiCapital)

          ! (needs to be changed) Can instead check if (self % phi /= 0_shortInt)?
          ! if phiCapital == 0, phi == 0 but not vise versa. For generasity, phiCapital is used here.
          ! If edge intersection test against edge A gives phiCaptital = 0, can the same test against 
          ! other edges can give non-zero phiCapital? (Although it requires mathematical proof, I do
          ! not think that is the case). If my guess is right, we can instead use "phi" here!
          if (self % phiCapital /= 0_shortInt) return

        end do

        ! If phiCapital mapping is not constructed, we perform the last test using intersected faces.
        ! (needs to be changed) is it possible to get non-zero phiCapital from this test in case phiCapital
        ! was not set from edgeIntersection test? If not, we do not have to perform this part of the code.
        do i = 1, size(intersectedFaceIdxs)-1
          faceIdxsToBeTested(1) = intersectedFaceIdxs(i)
          do j = i+1, size(intersectedFaceIdxs)

            faceIdxsToBeTested(2) = intersectedFaceIdxs(j)
            call testTwoIntersectedFaces(vertices, edges, faces, centroid, targetDistanceSqr, self % phi, &
                                         self % phiCapital, faceIdxsToBeTested)

            ! (needs to be changed) (for the same reason a few lines above)
            if (self % phiCapital /= 0_shortInt) return

          end do
        end do

      end if

    end if

  end subroutine constructNRefineCell

  !!
  !!
  !!
  function getPhiCapital(self) result(phiCapital)
    class(cartesianCellFinest), intent(in)              :: self
    integer(shortInt)                                   :: phiCapital

    !###
    if(allocated(self%phiCapital)) then
      phiCapital = self % phiCapital
    else
      phiCapital = 0
    end if

    ! phiCapital = self % phiCapital
    !###

  end function getPhiCapital

  !!
  !!
  !!
  function getPhi(self) result(phi)
    class(cartesianCellFinest), intent(in)              :: self
    integer(shortInt)                                   :: phi

    !###
    if(allocated(self%phi)) then
      phi = self % phi
    else
      phi = 0
    end if

    ! phi = self % phi
    !###

  end function getPhi

  !!
  !!
  !!
  function getChi(self) result(chi)
    class(cartesianCellFinest), intent(in)              :: self
    integer(shortInt)                                   :: chi

    !###
    if(allocated(self%chi)) then
      chi = self % chi
    else
      chi = 0
    end if

    ! chi = self % chi
    !###

  end function getChi

  !!
  !!
  !!
  function cellGetNumberOfCells(self, n_layers) result(output)
    class(cartesianCellFinest), intent(in)              :: self
    integer(shortInt), intent(in)                       :: n_layers
    integer(shortInt), dimension(:), allocatable        :: output

    ! Allocate and initialise output array
    allocate(output(n_layers+1))
    output = 0

    output(n_layers) = 1

    if (.NOT. allocated(self % chi)) then
      if (.NOT. allocated(self % phiCapital)) then
        output(n_layers + 1) = 1
      end if
    end if
    
  end function cellGetNumberOfCells


end module cartesianCellFinest_class