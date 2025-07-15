module noAcceleration_class

  use accelerationStructure_inter,  only : accelerationStructure
  use coord_class,                  only : coord
  use dictionary_class,             only : dictionary
  use element_class,                only : elementBox, inclusionTestResult
  use numPrecision
  use topologicalObjectShelf_class, only : topologicalObjectShelf
  use universalVariables,           only : INSIDE_ELEMENT, ON_BOUNDARY_ELEMENT, OUTSIDE_ELEMENT

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(accelerationStructure) :: noAcceleration
    private
  contains
    procedure :: findHostElement
    procedure :: init
    procedure :: kill
  end type noAcceleration

contains
  !!
  !!
  !!
  subroutine findHostElement(self, elements, coords)
    class(noAcceleration), intent(in)        :: self
    type(topologicalObjectShelf), intent(in) :: elements
    type(coord), intent(inout)               :: coords
    integer(shortInt)                        :: i
    type(elementBox)                         :: element
    type(inclusionTestResult)                :: testResult

    ! Perform brute-force search.
    searchLoop: do
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
            return

          elseif (testResult % status == OUTSIDE_ELEMENT) then
            ! If the nudge has resulted in an overshoot, we cycle searchLoop and begin the entire process again.
            cycle searchLoop

          end if

        end if

      end do
      ! If reached here, the particle is not inside any element.
      return

    end do searchLoop

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