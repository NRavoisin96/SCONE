module noAcceleration_class

  use accelerationStructure_inter,  only : accelerationStructure, initAccelerationStructurePayload
  use dictionary_class,             only : dictionary
  use face_class,                   only : face, faceBox
  use genericProcedures,            only : fatalError, numToChar
  use element_class,                only : element, elementBox, inclusionTestResult
  use numPrecision
  use publicObjects,                only : coordData, intersectionTestPayload, intersectionTestResult, newIntersectionTestPayload
  use topologicalObject_inter,      only : topologicalObjectBox
  use topologicalObjectShelf_class, only : topologicalObjectShelf
  use universalVariables,           only : INF, INSIDE_ELEMENT, NUDGE, ON_BOUNDARY_ELEMENT, OUTSIDE_ELEMENT

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(accelerationStructure) :: noAcceleration
    private
  contains
    procedure :: findEntranceBoundaryFace
    procedure :: findHostElement
    procedure :: init
    procedure :: kill
  end type noAcceleration

contains
  !!
  !!
  !!
  subroutine findEntranceBoundaryFace(self, faces, data, boundaryFace)
    class(noAcceleration), intent(in)        :: self
    type(topologicalObjectShelf), intent(in) :: faces
    type(coordData), intent(inout)           :: data
    type(faceBox), intent(out)               :: boundaryFace
    type(faceBox)                            :: testFace
    integer(shortInt)                        :: i
    real(defReal)                            :: dMax
    type(intersectionTestResult)             :: intersectionResult

    ! Assemble intersection test payload then loop through all faces in the geometry.
    do i = 1, faces % getObjectsNumber()
      testFace = faces % getFaceBox(i)
      ! Cycle to next face if current face is not active or not a boundary face.
      if (.not. (testFace % ptr % getIsActive() .and. testFace % ptr % getIsBoundary())) cycle

      ! Compute distance to boundary face.
      dMax = min(data % dMax, data % d)
      call testFace % ptr % intersects(newIntersectionTestPayload(data % r, data % u, dMax), intersectionResult)
      if (intersectionResult % intersects .and. intersectionResult % d < data % d) then
        data % d = intersectionResult % d
        boundaryFace = testFace

      end if

    end do

  end subroutine findEntranceBoundaryFace

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