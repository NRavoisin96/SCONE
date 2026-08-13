module octreeAcceleration_class

  use accelerationStructure_inter,  only : accelerationStructure, initAccelerationStructurePayload
  use dictionary_class,             only : dictionary
  use face_class,                   only : face, faceBox
  use element_class,                only : element, inclusionTestResult
  use genericProcedures,            only : fatalError, numToChar
  use kdTree_class,                 only : kdTree
  use kdTreeNode_class,             only : buildKDTreeNodePayload
  use node_inter,                   only : node
  use numPrecision
  use publicObjects,                only : coordData, intersectionTestResult, newIntersectionTestPayload
  use octree_class,                 only : octree
  use octreeNode_class,             only : buildOctreeNodePayload, octreeNode
  use topologicalObject_inter,      only : topologicalObjectBox
  use topologicalObjectShelf_class, only : topologicalObjectShelf
  use universalVariables,           only : INF, INSIDE_ELEMENT, NUDGE, ON_BOUNDARY_ELEMENT, OUTSIDE_ELEMENT, &
                                           TWO, VALENCE

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(accelerationStructure) :: octreeAcceleration
    private
    type(octree) :: tree
  contains
    procedure :: findEntranceBoundaryFace
    procedure :: findHostElement
    procedure :: init
    procedure :: kill
  end type octreeAcceleration

contains
  !!
  !!
  !!
  subroutine findEntranceBoundaryFace(self, faces, data, nIntersectedFaces, intersectedFaceIdxs)
    class(octreeAcceleration), intent(in)              :: self
    type(topologicalObjectShelf), intent(in)           :: faces
    type(coordData), intent(inout)                     :: data
    integer(shortInt), intent(out)                     :: nIntersectedFaces
    integer(shortInt), dimension(VALENCE), intent(out) :: intersectedFaceIdxs
    character(*), parameter                            :: HERE = 'distanceToBoundaryFace (octreeAcceleration_class.f90)'
    
    ! Call fatalError for now.
    nIntersectedFaces = 0
    intersectedFaceIdxs = 0
    call fatalError(HERE, 'Unsupported procedure.')

  end subroutine findEntranceBoundaryFace

  !!
  !!
  !!
  subroutine findHostElement(self, elements, data, stopSearch)
    class(octreeAcceleration), intent(in)                 :: self
    type(topologicalObjectShelf), intent(in)              :: elements
    type(coordData), intent(inout)                        :: data
    logical(defBool), intent(out)                         :: stopSearch
    type(topologicalObjectBox), dimension(:), allocatable :: objects
    integer(shortInt)                                     :: i, nPotentialElements
    class(node), pointer                                  :: genericLeaf
    type(octreeNode), pointer                             :: leaf
    type(inclusionTestResult)                             :: testResult
    character(*), parameter                               :: here = 'findHostElement (octreeAcceleration_class.f90)'

    ! First search the acceleration structure for the indices of potential elements containing the coordinates.
    stopSearch = .true.
    call self % tree % findLeaf(data % r, data % u, genericLeaf)

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
          data % elementIdx = ptr % getIdx()
          data % localId = ptr % getLocalId()
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
            testResult = ptr % isPointInside(data % r)

            if (testResult % status == INSIDE_ELEMENT) then
              ! If coordinates are fully inside, we have found our element.
              data % elementIdx = ptr % getIdx()
              data % localId = ptr % getLocalId()
              return

            elseif (testResult % status == ON_BOUNDARY_ELEMENT) then
              ! If coordinates are on the element boundary (very rare), we need to push them off.
              do while (testResult % status == ON_BOUNDARY_ELEMENT)
                call ptr % pushFromBoundary(data % u, data % r)

                ! Perform containment test again.
                testResult = ptr % isPointInside(data % r)

              end do

              ! Now the coordinates are not on the boundary of the element anymore.
              if (testResult % status == INSIDE_ELEMENT) then
                ! If coordinates are now well inside the element, we have found our element.
                data % elementIdx = ptr % getIdx()
                data % localId = ptr % getLocalId()

              elseif (testResult % status == OUTSIDE_ELEMENT) then
                ! If the nudge has resulted in an overshoot, we cycle searchLoop and begin the entire process again.
                stopSearch = .false.

              end if
              return

            end if

          class default
            call fatalError(here, 'Object: '//numToChar(ptr % getIdx())//' is not an element.')

        end select

      end do

    end if

  end subroutine findHostElement

  !!
  !!
  !!
  subroutine init(self, payload)
    class(octreeAcceleration), intent(inout)           :: self
    type(initAccelerationStructurePayload), intent(in) :: payload
    type(kdTree), target                               :: tree
    type(buildKDTreeNodePayload)                       :: kdTreePayload
    type(buildOctreeNodePayload)                       :: octreePayload

    ! Initialise kd-tree using the faceShelf.
    kdTreePayload % shelf => payload % faces
    kdTreePayload % idxs = kdTreePayload % shelf % getActiveObjectIdxs()
    kdTreePayload % lowerBound = 1
    kdTreePayload % upperBound = size(kdTreePayload % idxs)
    kdTreePayload % bucketSize = 4
    call tree % init(kdTreePayload)
    
    ! Build payload then initialise octree.
    octreePayload % shelf => payload % faces
    octreePayload % computeLeafBoundingBox = .false.
    octreePayload % updateParentBoundingBox = .false.
    octreePayload % tree => tree
    octreePayload % bounds = octreePayload % tree % getRootBoundingBoxBounds() + &
                             reshape(TWO * [-NUDGE, -NUDGE, -NUDGE, NUDGE, NUDGE, NUDGE], [3, 2])
    call self % tree % init(octreePayload, payload % dict)

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