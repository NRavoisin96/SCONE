module octreeAcceleration_class

  use accelerationStructure_inter, only : accelerationStructure
  use coord_class,                 only : coord
  use dictionary_class,            only : dictionary
  use element_inter,               only : inclusionTestResult
  use elementShelf_class,          only : elementShelf
  use faceShelf_class,             only : faceShelf
  use numPrecision
  use octree_class,                only : octree
  use octreeNode_class,            only : octreeNode
  use vertexShelf_class,           only : vertexShelf
  use edgeShelf_class,             only : edgeShelf
  use universalVariables,          only : INSIDE_ELEMENT, ON_BOUNDARY_ELEMENT, OUTSIDE_ELEMENT

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(accelerationStructure) :: octreeAcceleration
    private
    type(octree) :: tree
  contains
    procedure :: findHostElement
    procedure :: init
    procedure :: kill
  end type octreeAcceleration

contains
  !!
  !!
  !!
  subroutine findHostElement(self, vertices, edges, faces, elements, coords)
    class(octreeAcceleration), intent(in)        :: self
    class(vertexShelf), intent(in)               :: vertices
    class(edgeShelf), intent(in)                 :: edges
    type(faceShelf), intent(in)                  :: faces
    type(elementShelf), intent(in)               :: elements
    type(coord), intent(inout)                   :: coords
    integer(shortInt), dimension(:), allocatable :: potentialElementIdxs
    integer(shortInt)                            :: i, nPotentialElements, potentialElementIdx
    type(octreeNode), pointer                    :: leaf
    type(inclusionTestResult)                    :: testResult

    searchLoop: do
      ! First search the acceleration structure for the indices of potential elements containing the coordinates.
      if (allocated(potentialElementIdxs)) deallocate(potentialElementIdxs)
      call self % tree % findLeaf(coords, leaf)

      ! If check if the leaf node pointer is associated.
      if (.not. associated(leaf)) return
      if (leaf % getIsOutside()) return
      
      ! Retrieve the element indices in the leaf.
      potentialElementIdxs = leaf % getElementIdxs()

      if (leaf % getIsInside()) then
        potentialElementIdx = potentialElementIdxs(1)
        call coords % setElementIdx(potentialElementIdx)
        call coords % setParentElementIdx(elements % getElementParentIdx(potentialElementIdx))
        return

      end if

      if (leaf % getIsIntersecting()) then
        nPotentialElements = size(potentialElementIdxs)
        do i = 1, nPotentialElements
          ! Perform inclusion test for the current element.
          potentialElementIdx = potentialElementIdxs(i)
          testResult = elements % isPointInside(potentialElementIdx, coords % getPositionToNudge(), faces)

          if (testResult % status == INSIDE_ELEMENT) then
            ! If coordinates are fully inside, we have found our element.
            call coords % setElementIdx(potentialElementIdx)
            call coords % setParentElementIdx(elements % getElementParentIdx(potentialElementIdx))
            return

          elseif (testResult % status == ON_BOUNDARY_ELEMENT) then
            ! If coordinates are on the element boundary (very rare), we need to push them off.
            do while (testResult % status == ON_BOUNDARY_ELEMENT)
              call elements % pushFromElementBoundary(potentialElementIdx, faces, coords)

              ! Perform containment test again.
              testResult = elements % isPointInside(potentialElementIdx, coords % getPositionToNudge(), faces)

            end do

            ! Now the coordinates are not on the boundary of the element anymore.
            if (testResult % status == INSIDE_ELEMENT) then
              ! If coordinates are now well inside the element, we have found our element.
              call coords % setElementIdx(potentialElementIdx)
              call coords % setParentElementIdx(elements % getElementParentIdx(potentialElementIdx))
              return

            elseif (testResult % status == OUTSIDE_ELEMENT) then
              ! If the nudge has resulted in an overshoot, we cycle searchLoop and begin the entire process again.
              cycle searchLoop

            end if

          end if

        end do
        ! If reached here, the coordinates are not inside any elements so they are outside the mesh. Simply return here.
        return

      end if

    end do searchLoop

  end subroutine findHostElement

  !!
  !!
  !!
  subroutine init(self, dict, vertices, edges, faces, elements)
    class(octreeAcceleration), intent(inout) :: self
    type(dictionary), intent(in)             :: dict
    type(vertexShelf), intent(in)            :: vertices
    type(faceShelf), intent(inout)           :: faces
    type(elementShelf), intent(in)           :: elements
    type(edgeShelf), intent(inout)           :: edges

    ! Simply initialise the octree.
    call self % tree % init(vertices, faces, elements)

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(octreeAcceleration), intent(inout) :: self

    ! Local.
    call self % tree % kill()

  end subroutine kill

end module octreeAcceleration_class