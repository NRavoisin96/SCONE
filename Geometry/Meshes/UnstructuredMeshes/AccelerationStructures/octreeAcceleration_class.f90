module octreeAcceleration_class

  use accelerationStructure_inter,  only : accelerationStructure
  use coord_class,                  only : coord
  use dictionary_class,             only : dictionary
  use element_class,                only : element, inclusionTestResult
  use genericProcedures,            only : fatalError, numToChar
  use kdTree_class,                 only : kdTree
  use kdTreeNode_class,             only : buildKDTreeNodePayload
  use node_inter,                   only : node
  use numPrecision
  use octree_class,                 only : octree
  use octreeNode_class,             only : buildOctreeNodePayload, octreeNode
  use topologicalObject_inter,      only : topologicalObjectBox
  use topologicalObjectShelf_class, only : topologicalObjectShelf
  use universalVariables,           only : INSIDE_ELEMENT, NUDGE, ON_BOUNDARY_ELEMENT, OUTSIDE_ELEMENT, TWO

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
  subroutine findHostElement(self, elements, coords)
    class(octreeAcceleration), intent(in)                 :: self
    type(topologicalObjectShelf), intent(in)              :: elements
    type(coord), intent(inout)                            :: coords
    type(topologicalObjectBox), dimension(:), allocatable :: objects
    integer(shortInt)                                     :: i, nPotentialElements, elementIdx
    class(node), pointer                                  :: genericLeaf
    type(octreeNode), pointer                             :: leaf
    type(inclusionTestResult)                             :: testResult
    character(*), parameter                               :: here = 'findHostElement (octreeAcceleration_class.f90)'

    searchLoop: do
      ! First search the acceleration structure for the indices of potential elements containing the coordinates.
      call self % tree % findLeaf(coords, genericLeaf, .true.)

      ! If check if the leaf node pointer is associated.
      if (.not. associated(genericLeaf)) return

      ! Downcast leaf to correct type.
      select type(ptr => genericLeaf)
        type is(octreeNode)
          leaf => ptr

        class default
          call fatalError(here, 'Invalid leaf node type.')

      end select

      ! Call fatalError if leaf is unchecked.
      if (leaf % getIsUnchecked()) call fatalError(here, 'Leaf node is neither fully inside, fully outside, nor intersecting.')
      if (leaf % getIsOutside()) return
      
      ! Retrieve the element in the leaf.
      objects = leaf % getContainingObjects()

      if (leaf % getIsInside()) then
        ! Downcast objects to correct type.
        select type(ptr => objects(1) % ptr)
          type is(element)
            call coords % setElementIdx(ptr % getIdx())
            call coords % setParentElementIdx(ptr % getParentIdx())
            call coords % setLocalId(ptr % getLocalId())
            return

          class default
            call fatalError(here, 'Object: '//numToChar(ptr % getIdx())//' is not an element.')

        end select

      end if

      if (leaf % getIsIntersecting()) then
        nPotentialElements = size(objects)
        do i = 1, nPotentialElements
          ! Downcast current object to correct type.
          select type(ptr => objects(i) % ptr)
            type is(element)
              ! Perform inclusion test for the current element.
              elementIdx = ptr % getIdx()
              testResult = ptr % isPointInside(coords % getPositionToNudge())

              if (testResult % status == INSIDE_ELEMENT) then
                ! If coordinates are fully inside, we have found our element.
                call coords % setElementIdx(elementIdx)
                call coords % setParentElementIdx(ptr % getParentIdx())
                call coords % setLocalId(ptr % getLocalId())
                return

              elseif (testResult % status == ON_BOUNDARY_ELEMENT) then
                ! If coordinates are on the element boundary (very rare), we need to push them off.
                do while (testResult % status == ON_BOUNDARY_ELEMENT)
                  call ptr % pushFromBoundary(coords)

                  ! Perform containment test again.
                  testResult = ptr % isPointInside(coords % getPositionToNudge())

                end do

                ! Now the coordinates are not on the boundary of the element anymore.
                if (testResult % status == INSIDE_ELEMENT) then
                  ! If coordinates are now well inside the element, we have found our element.
                  call coords % setElementIdx(elementIdx)
                  call coords % setParentElementIdx(ptr % getParentIdx())
                  call coords % setLocalId(ptr % getLocalId())
                  return

                elseif (testResult % status == OUTSIDE_ELEMENT) then
                  ! If the nudge has resulted in an overshoot, we cycle searchLoop and begin the entire process again.
                  cycle searchLoop

                end if

              end if

            class default
              call fatalError(here, 'Object: '//numToChar(ptr % getIdx())//' is not an element.')

          end select

        end do
        ! If reached here, the coordinates are not inside any elements so they are outside the mesh. Simply return here.
        return

      end if

    end do searchLoop

  end subroutine findHostElement

  !!
  !!
  !!
  subroutine init(self, dict, edges, elements, faces, vertices)
    class(octreeAcceleration), intent(inout)         :: self
    class(dictionary), intent(in)                    :: dict
    type(topologicalObjectShelf), target, intent(in) :: edges, elements, faces, vertices
    type(kdTree), target                             :: tree
    type(buildKDTreeNodePayload)                     :: kdTreePayload
    integer(shortInt)                                :: i
    type(buildOctreeNodePayload)                     :: octreePayload

    ! Initialise kd-tree using the faceShelf.
    kdTreePayload % shelf => faces
    kdTreePayload % idxs = kdTreePayload % shelf % getActiveObjectIdxs()
    kdTreePayload % lowerBound = 1
    kdTreePayload % upperBound = size(kdTreePayload % idxs)
    kdTreePayload % bucketSize = 4
    call tree % init(kdTreePayload)
    
    ! Build payload then initialise octree.
    octreePayload % shelf => faces
    octreePayload % computeLeafBoundingBox = .false.
    octreePayload % updateParentBoundingBox = .false.
    octreePayload % tree => tree
    octreePayload % bounds = octreePayload % tree % getRootBoundingBoxBounds() + &
                             reshape(TWO * [-NUDGE, -NUDGE, -NUDGE, NUDGE, NUDGE, NUDGE], [3, 2])
    call self % tree % init(octreePayload, dict)

    ! Assign non-intersecting cells to elements.
    call self % tree % assignElements(octreePayload)

    ! Kill the kd-tree as it is no longer required.
    call octreePayload % tree % kill()

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