module noAcceleration_class

  use accelerationStructure_inter,  only : accelerationStructure
  use coord_class,                  only : coord
  use dictionary_class,             only : dictionary
  use face_class,                   only : face, faceBox
  use genericProcedures,            only : fatalError, numToChar
  use element_class,                only : element, elementBox, inclusionTestResult
  use numPrecision
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
  subroutine findEntranceBoundaryFace(self, faces, coords, d, boundaryFace)
    class(noAcceleration), intent(in)        :: self
    type(topologicalObjectShelf), intent(in) :: faces
    type(coord), intent(in)                  :: coords
    real(defReal), intent(inout)             :: d
    type(faceBox), intent(out)               :: boundaryFace
    type(faceBox)                            :: testFace
    integer(shortInt)                        :: i
    real(defReal)                            :: update

    do i = 1, faces % getObjectsNumber()
      testFace = faces % getFaceBox(i)
      ! Cycle to next face if current face is not active or not a boundary face.
      if (.not. (testFace % ptr % getIsActive() .and. testFace % ptr % getIsBoundary())) cycle

      ! Compute distance to boundary face.
      call testFace % ptr % computeIntersection(coords, update)
      if (update < d) then
        d = update
        boundaryFace = testFace

      end if

    end do

  end subroutine findEntranceBoundaryFace

  !!
  !!
  !!
  subroutine findHostElement(self, elements, coords, stopSearch)
    class(noAcceleration), intent(in)        :: self
    type(topologicalObjectShelf), intent(in) :: elements
    type(coord), intent(inout)               :: coords
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
      
      testResult = element % ptr % isPointInside(coords % getPositionToNudge())
      if (testResult % status == INSIDE_ELEMENT) then
        call coords % setElementIdx(element % ptr % getIdx())
        call coords % setParentElementIdx(element % ptr % getParentIdx())
        call coords % setLocalId(element % ptr % getLocalId())
        return 

      elseif (testResult % status == ON_BOUNDARY_ELEMENT) then
        ! If coordinates are on the element boundary (very rare), we need to push them off.
        do while (testResult % status == ON_BOUNDARY_ELEMENT)
          call element % ptr % pushFromBoundary(coords)

          ! Perform containment test again.
          testResult = element % ptr % isPointInside(coords % getPositionToNudge())

        end do

        ! Now the coordinates are not on the boundary of the element anymore.
        if (testResult % status == INSIDE_ELEMENT) then
          ! If coordinates are now well inside the element, we have found our element.
          call coords % setElementIdx(element % ptr % getIdx())
          call coords % setParentElementIdx(element % ptr % getParentIdx())
          call coords % setLocalId(element % ptr % getLocalId())

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
  subroutine init(self, dict, edges, elements, faces, vertices)
    class(noAcceleration), intent(inout)             :: self
    class(dictionary), intent(in)                    :: dict
    type(topologicalObjectShelf), target, intent(in) :: edges, elements, faces, vertices

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