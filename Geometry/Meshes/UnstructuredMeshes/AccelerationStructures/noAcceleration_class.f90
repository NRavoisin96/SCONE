module noAcceleration_class

  use accelerationStructure_inter,  only : accelerationStructure, initAccelerationStructurePayload
  use dictionary_class,             only : dictionary
  use face_class,                   only : face, faceBox
  use genericProcedures,            only : areEqual, fatalError, numToChar
  use element_class,                only : element, elementBox, inclusionTestResult
  use numPrecision
  use publicObjects,                only : coordData, intersectionTestPayload, intersectionTestResult, &
                                           newIntersectionTestPayload, newRationalIntersectionTestPayload, &
                                           rationalIntersectionTestPayload
  use ratint
  use topologicalObject_inter,      only : topologicalObjectBox
  use topologicalObjectShelf_class, only : topologicalObjectShelf
  use universalVariables,           only : INF, INSIDE_ELEMENT, NUDGE, ON_BOUNDARY_ELEMENT, OUTSIDE_ELEMENT, VALENCE

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(accelerationStructure) :: noAcceleration
    private
  contains
    procedure :: findEntranceBoundaryFace
    procedure :: findEntranceBoundaryFace_rational
    procedure :: findHostElement
    procedure :: init
    procedure :: kill
  end type noAcceleration

contains
  !!
  !!
  !!
  subroutine findEntranceBoundaryFace(self, faces, data, nIntersectedFaces, intersectedFaceIdxs)
    class(noAcceleration), intent(in)                  :: self
    type(topologicalObjectShelf), intent(in)           :: faces
    type(coordData), intent(inout)                     :: data
    integer(shortInt), intent(out)                     :: nIntersectedFaces
    integer(shortInt), dimension(VALENCE), intent(out) :: intersectedFaceIdxs
    integer(shortInt)                                  :: faceIdx, i, j
    real(defReal)                                      :: dMin
    real(defReal), dimension(VALENCE)                  :: distances
    type(faceBox)                                      :: testFace
    type(intersectionTestPayload)                      :: payload
    type(intersectionTestResult)                       :: intersectionResult
    character(*), parameter                            :: HERE = 'findEntranceBoundaryFace (noAcceleration_class.f90)'

    ! Initialise nIntersectedFaces = 0, intersectedFaceIdxs = 0, and dMin = INF.
    nIntersectedFaces = 0
    intersectedFaceIdxs = 0
    dMin = INF

    ! Initialise intersection test payload.
    payload = newIntersectionTestPayload(data % r, data % u, data % dMax)

    ! Assemble intersection test payload then loop through all faces in the geometry to collect a first set of
    ! potentially intersected boundary faces.
    do i = 1, faces % getObjectsNumber()
      testFace = faces % getFaceBox(i)
      ! Cycle to next face if current face is not active or not a boundary face.
      if(.not. (testFace % ptr % getIsActive() .and. testFace % ptr % getIsBoundary())) cycle

      ! Cycle if the current face is one of the faces the particle is currently on.
      faceIdx = testFace % ptr % getIdx()
      if(0 < data % front) then
        if(any(data % currentFaceIdxs(1:data % front) == faceIdx)) cycle

      end if

      ! Compute distance to boundary face and move onto the next face if it is not intersected.
      call testFace % ptr % intersects_Ray(payload, intersectionResult)

      ! Escalate to exact computation if needed.
      if(intersectionResult % needsRescue) then
        call self % findEntranceBoundaryFace_rational(faces, data, nIntersectedFaces, intersectedFaceIdxs)
        return

      elseif(.not. intersectionResult % intersects) then
        cycle

      end if

      if(nIntersectedFaces == 0) then
        ! First intersected face found.
        nIntersectedFaces = 1
        intersectedFaceIdxs(1) = faceIdx
        distances(1) = intersectionResult % d
        dMin = intersectionResult % d

      elseif(areEqual(intersectionResult % d, dMin)) then
        ! If distance to the current face is within tolerance of the best distance found so far,
        ! add the current face to the list of intersected faces.
        nIntersectedFaces = nIntersectedFaces + 1

        ! Call fatalError if the number of intersected faces exceeds VALENCE.
        if(VALENCE < nIntersectedFaces) &
          call fatalError(HERE, 'The number of intersected boundary faces exceeds '//numToChar(VALENCE)//'.')

        intersectedFaceIdxs(nIntersectedFaces) = faceIdx
        distances(nIntersectedFaces) = intersectionResult % d
        dMin = min(intersectionResult % d, dMin)

      elseif(intersectionResult % d < dMin) then
        ! Else, reset counter and list of intersected faces.
        nIntersectedFaces = 1
        intersectedFaceIdxs(1) = faceIdx
        distances(1) = intersectionResult % d
        dMin = intersectionResult % d

      end if

    end do

    ! Now prune the set of potentially intersected faces.
    j = 0
    do i = 1, nIntersectedFaces
      ! Reject this face as it is too far away from the minimum distance.
      if(.not. areEqual(distances(i), dMin)) cycle

      ! Update counter and set of intersected faces.
      j = j + 1
      intersectedFaceIdxs(j) = intersectedFaceIdxs(i)

    end do
    nIntersectedFaces = j

    ! If intersected faces have been found, update data % d.
    if(0 < nIntersectedFaces) data % d = dMin

  end subroutine findEntranceBoundaryFace

  !!
  !!
  !!
  subroutine findEntranceBoundaryFace_rational(self, faces, data, nIntersectedFaces, intersectedFaceIdxs)
    class(noAcceleration), intent(in)                  :: self
    type(topologicalObjectShelf), intent(in)           :: faces
    type(coordData), intent(inout)                     :: data
    integer(shortInt), intent(out)                     :: nIntersectedFaces
    integer(shortInt), dimension(VALENCE), intent(out) :: intersectedFaceIdxs
    integer(shortInt)                                  :: faceIdx, i, j
    type(faceBox)                                      :: testFace
    type(intersectionTestResult)                       :: intersectionResult
    type(ratint_t)                                     :: d
    type(ratint_t), dimension(VALENCE)                 :: distances
    type(rationalIntersectionTestPayload)              :: payload
    character(*), parameter                            :: HERE = 'findEntranceBoundaryFace_exact (noAcceleration_class.f90)'

    ! Initialise nIntersectedFaces = 0, intersectedFaceIdxs = 0, and dMin = INF.
    nIntersectedFaces = 0
    intersectedFaceIdxs = 0
    d = def_ratint_large()

    ! Initialise intersection test payload.
    payload = newRationalIntersectionTestPayload(convert_ieee(data % r), convert_ieee(data % u), convert_ieee(data % dMax))

    ! Assemble intersection test payload then loop through all faces in the geometry to collect a first set of
    ! potentially intersected boundary faces.
    do i = 1, faces % getObjectsNumber()
      testFace = faces % getFaceBox(i)
      ! Cycle to next face if current face is not active or not a boundary face.
      if(.not. (testFace % ptr % getIsActive() .and. testFace % ptr % getIsBoundary())) cycle

      ! Cycle if the current face is one of the faces the particle is currently on.
      faceIdx = testFace % ptr % getIdx()
      if(0 < data % front) then
        if(any(data % currentFaceIdxs(1:data % front) == faceIdx)) cycle

      end if

      ! Compute distance to boundary face and move onto the next face if it is not intersected.
      call testFace % ptr % intersects_Ray_rational(payload, intersectionResult)

      ! Escalate to exact computation if needed.
      if(.not. intersectionResult % intersects) cycle

      if(nIntersectedFaces == 0) then
        ! First intersected face found.
        nIntersectedFaces = 1
        intersectedFaceIdxs(1) = faceIdx
        distances(1) = intersectionResult % d_rational
        d = intersectionResult % d_rational

      elseif(intersectionResult % d_rational == d) then
        ! If distance to the current face is within tolerance of the best distance found so far,
        ! add the current face to the list of intersected faces.
        nIntersectedFaces = nIntersectedFaces + 1

        ! Call fatalError if the number of intersected faces exceeds VALENCE.
        if(VALENCE < nIntersectedFaces) &
          call fatalError(HERE, 'The number of intersected boundary faces exceeds '//numToChar(VALENCE)//'.')

        intersectedFaceIdxs(nIntersectedFaces) = faceIdx
        distances(nIntersectedFaces) = intersectionResult % d_rational

      elseif(d > intersectionResult % d_rational) then
        ! Else, reset counter and list of intersected faces.
        nIntersectedFaces = 1
        intersectedFaceIdxs(1) = faceIdx
        distances(1) = intersectionResult % d_rational
        d = intersectionResult % d_rational

      end if

    end do

    ! If intersected faces have been found, update data % d.
    if(0 < nIntersectedFaces) data % d = evaluate(d)

  end subroutine findEntranceBoundaryFace_rational

  !!
  !!
  !!
  subroutine findHostElement(self, elements, data, stopSearch)
    class(noAcceleration), intent(in)        :: self
    type(topologicalObjectShelf), intent(in) :: elements
    type(coordData), intent(inout)           :: data
    logical(defBool), intent(out)            :: stopSearch
    integer(shortInt)                        :: i
    type(elementBox)                         :: element
    type(inclusionTestResult)                :: testResult

    ! Perform brute-force search.
    stopSearch = .true.
    do i = 1, elements % getObjectsNumber()
      element = elements % getElementBox(i)
      ! Cycle to the next element if the current element is not active.
      if (.not. element % ptr % getIsActive()) cycle
      
      testResult = element % ptr % isPointInside(data % r)
      if (testResult % status == INSIDE_ELEMENT) then
        data % elementIdx = element % ptr % getIdx()
        data % localId = element % ptr % getLocalId()
        return 

      elseif (testResult % status == ON_BOUNDARY_ELEMENT) then
        ! If coordinates are on the element boundary (very rare), we need to push them off.
        do while (testResult % status == ON_BOUNDARY_ELEMENT)
          call element % ptr % pushFromBoundary(data % u, data % r)

          ! Perform containment test again.
          testResult = element % ptr % isPointInside(data % r)

        end do

        ! Now the coordinates are not on the boundary of the element anymore.
        if (testResult % status == INSIDE_ELEMENT) then
          ! If coordinates are now well inside the element, we have found our element.
          data % elementIdx = element % ptr % getIdx()
          data % localId = element % ptr % getLocalId()

        elseif (testResult % status == OUTSIDE_ELEMENT) then
          stopSearch = .false.

        end if
        return

      end if

    end do

  end subroutine findHostElement

  !!
  !!
  !!
  subroutine init(self, payload)
    class(noAcceleration), intent(inout)               :: self
    type(initAccelerationStructurePayload), intent(in) :: payload

    ! Do nothing.

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(noAcceleration), intent(inout) :: self

    ! Do nothing.

  end subroutine kill

end module noAcceleration_class