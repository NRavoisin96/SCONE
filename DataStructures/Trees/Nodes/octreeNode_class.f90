module octreeNode_class

  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use element_class,                only : element, elementBox, inclusionTestResult
  use face_class,                   only : face
  use genericProcedures,            only : append, areEqual, fatalError, numToChar
  use kdTree_class,                 only : kdTree
  use node_inter,                   only : buildNodePayload, kill_super => kill, node, nodeBox
  use numPrecision
  use publicObjects,                only : coordData, intersectionTestPayload, intersectionTestResult, newIntersectionTestPayload
  use topologicalObject_inter,      only : getUniqueSharingElements, topologicalObjectBox
  use topologicalObjectShelf_class, only : topologicalObjectShelf
  use universalVariables,           only : HALF, INF, INSIDE_ELEMENT, ON_BOUNDARY_ELEMENT, OUTSIDE_ELEMENT

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(buildNodePayload)  :: buildOctreeNodePayload
    real(defReal), dimension(3, 2)         :: bounds
    type(kdTree), pointer                  :: tree => null()
  end type buildOctreeNodePayload

  !!
  !!
  !!
  type, public, extends(node) :: octreeNode
    private
    logical(defBool)          :: isUnchecked = .true., isInside = .false., isOutside = .false., isIntersecting = .false.
  contains
    ! Build procedures.
    procedure :: allocateChild
    procedure :: assignElements
    procedure :: build
    procedure :: getChildrenNumber
    procedure :: kill
    procedure :: preparePayloadForChild
    ! Runtime procedures.
    procedure :: countInside
    procedure :: countOutside
    procedure :: findFirstIntersectedObject
    procedure :: findNearestObject
    procedure :: getDescentChildIdx
    procedure :: getIsInside
    procedure :: getIsIntersecting
    procedure :: getIsOutside
    procedure :: getIsUnchecked
  end type octreeNode

contains
  !!
  !!
  !!
  subroutine allocateChild(self, ptr)
    class(octreeNode), intent(in)     :: self
    class(node), pointer, intent(out) :: ptr

    allocate(octreeNode :: ptr)

  end subroutine allocateChild

  !!
  !!
  !!
  recursive subroutine assignElements(self, payload)
    class(octreeNode), intent(inout)                      :: self
    class(buildNodePayload), intent(in)                   :: payload
    type(nodeBox), dimension(8)                           :: children
    integer(shortInt)                                     :: i, nElements
    type(buildOctreeNodePayload), pointer                 :: payloadPtr
    type(topologicalObjectBox), dimension(:), allocatable :: elements
    real(defReal), dimension(3)                           :: boundingBoxCentre
    type(topologicalObjectBox)                            :: nearestFace
    type(inclusionTestResult)                             :: insideResult
    character(*), parameter                               :: here = 'assignElements (octreeNode_class.f90)'

    ! If the cell is not a leaf, descend deeper into the tree.
    if (.not. self % getIsLeaf()) then
      children = self % getChildren()
      do i = 1, 8
        if (.not. associated(children(i) % ptr)) call fatalError(here, 'Unassociated pointer for child: '//numToChar(i)//'.')
        select type(ptr => children(i) % ptr)
          type is(octreeNode)
            call ptr % assignElements(payload)

          class default
            call fatalError(here, 'Invalid node type for child: '//numToChar(i)//'.')

        end select

      end do
      return

    end if

    ! Downcast payload to correct type.
    select type(ptr => payload)
      type is(buildOctreeNodePayload)
        payloadPtr => ptr

      class default
        call fatalError(here, 'Invalid payload type.')

    end select

    ! If the cell has already been checked simply return.
    if (.not. self % isUnchecked) return

    ! If the cell is unchecked, find the nearest mesh face from the tree.
    self % isUnchecked = .false.
    boundingBoxCentre = self % getBoundingBoxCentre()
    nearestFace = payloadPtr % tree % findNearestObject(boundingBoxCentre)

    ! Retrieve the elements associated with the nearest face.
    elements = nearestFace % ptr % getSharingElements()
    nElements = size(elements)
    if (nElements == 0) &
    call fatalError(here, 'Unable to retrieve elements associated with face: '//numToChar(nearestFace % ptr % getIdx())//'.')
    
    do i = 1, nElements
      ! Downcast current element to correct type.
      select type(ptr => elements(i) % ptr)
        type is(element)
          ! Check for inclusion in the current element.
          insideResult = ptr % isPointInside(boundingBoxCentre)

          ! If the current element contains the bounding box, assign it to the cell
          ! and return.
          if (insideResult % status == INSIDE_ELEMENT) then
            call self % addContainingObject(elements(i))
            self % isInside = .true.
            return

          end if

        class default
          call fatalError(here, 'Element :'//numToChar(ptr % getIdx())//' associated with face: '&
                          //numToChar(nearestFace % ptr % getIdx())//' is not an element.')


      end select

    end do

    ! If reached here, the cell is not in any element, so set it as being outside.
    self % isOutside = .true.

  end subroutine assignElements

  !!
  !!
  !!
  subroutine build(self, payload, stop)
    class(octreeNode), intent(inout)                      :: self
    class(buildNodePayload), intent(inout)                :: payload
    logical(defBool), intent(out)                         :: stop
    type(buildOctreeNodePayload), pointer                 :: payloadPtr
    integer(shortInt)                                     :: i, nFaces, nIntersectedFaces
    type(topologicalObjectBox), dimension(:), allocatable :: elements, intersectedFaces
    character(*), parameter                               :: here = 'build (octreeNode_class.f90)'

    ! Initialise stop = .false.
    stop = .false.
    
    ! Downgrade payload type.
    select type(ptr => payload)
      type is (buildOctreeNodePayload)
        payloadPtr => ptr

      class default
        call fatalError(here, 'Invalid payload type.')

    end select

    ! Check that accelerator tree is correctly associated.
    if (.not. associated(payloadPtr % tree)) call fatalError(here, 'Unassociated tree pointer.')
    
    ! Set the node's bounding box.
    call self % initBoundingBox(payloadPtr % bounds)

    ! Use simplified logic for the root node.
    if (self % getDepth() == 1) then
      self % isUnchecked = .false.
      self % isIntersecting = .true.
      ! Retrieve the number of faces in the tree.
      nFaces = payloadPtr % tree % getObjectsNumber()

      ! If nFaces <= maxFacesNumber (very simple unstructured mesh geometries), there is no need
      ! to refine the cell and we can simply return.
      if (nFaces <= self % getBucketSize()) then
        intersectedFaces = payloadPtr % shelf % getShelf()
        elements = getUniqueSharingElements(intersectedFaces)
        do i = 1, size(elements)
          ! Downcast element to correct type.
          select type(ptr => elements(i) % ptr)
            type is(element)
              ! Do nothing.

            class default
              call fatalError(here, 'Element: '//numToChar(ptr % getIdx())//' is not an element.')

          end select

        end do
        call self % addContainedObject(intersectedFaces)
        call self % addContainingObject(elements)
        stop = .true.

      end if
      return

    end if

    ! Check the number of intersections between the current cell and the faces in the mesh by traversing the
    ! k-d tree starting from the root node.
    intersectedFaces = payloadPtr % tree % findIntersectedObjects(self % getBoundingBoxPtr())
    nIntersectedFaces = size(intersectedFaces)
    if (nIntersectedFaces == 0) then
      stop = .true.

    else
      self % isUnchecked = .false.
      self % isIntersecting = .true.
      if (self % getDepth() == payloadPtr % maxDepth .or. nIntersectedFaces <= self % getBucketSize()) then
        elements = getUniqueSharingElements(intersectedFaces)
        do i = 1, size(elements)
          select type(ptr => elements(i) % ptr)
            type is(element)
              ! Ok, do nothing.

            class default
              call fatalError(here, 'Element: '//numToChar(ptr % getIdx())//' is not an element.')

          end select

        end do
        call self % addContainedObject(intersectedFaces)
        call self % addContainingObject(elements)
        stop = .true.

      end if

    end if

  end subroutine build

  !!
  !!
  !!
  recursive subroutine countInside(self, nInside)
    class(octreeNode), intent(in)    :: self
    integer(shortInt), intent(inout) :: nInside
    type(nodeBox), dimension(8)      :: children
    integer(shortInt)                :: i
    character(*), parameter          :: here = 'countInside (octreeNode_class.f90)'

    if (.not. self % getIsLeaf()) then
      children = self % getChildren()
      do i = 1, 8
        if (.not. associated(children(i) % ptr)) call fatalError(here, 'Unassociated pointer for child: '//numToChar(i)//'.')
        select type(ptr => children(i) % ptr)
          type is(octreeNode)
            call ptr % countInside(nInside)

          class default
            call fatalError(here, 'Invalid node type for child: '//numToChar(i)//'.')

        end select

      end do
      return

    end if

    if (self % isInside) nInside = nInside + 1

  end subroutine countInside

  !!
  !!
  !!
  recursive subroutine countOutside(self, nOutside)
    class(octreeNode), intent(in)    :: self
    integer(shortInt), intent(inout) :: nOutside
    type(nodeBox), dimension(8)      :: children
    integer(shortInt)                :: i
    character(*), parameter          :: here = 'countOutside (octreeNode_class.f90)'

    if (.not. self % getIsLeaf()) then
      children = self % getChildren()
      do i = 1, 8
        if (.not. associated(children(i) % ptr)) call fatalError(here, 'Unassociated pointer for child: '//numToChar(i)//'.')
        select type(ptr => children(i) % ptr)
          type is(octreeNode)
            call ptr % countOutside(nOutside)

          class default
            call fatalError(here, 'Invalid node type for child: '//numToChar(i)//'.')

        end select

      end do
      return

    end if

    if (self % isOutside) nOutside = nOutside + 1

  end subroutine countOutside

  !!
  !!
  !!
  subroutine findFirstIntersectedObject(self, data, firstIntersectedObject)
    class(octreeNode), intent(in)                         :: self
    type(coordData), intent(inout)                        :: data
    type(topologicalObjectBox), intent(out)               :: firstIntersectedObject
    type(axisAlignedBoundingBox), pointer                 :: boundingBoxPtr
    type(intersectionTestResult)                          :: boundingBoxIntersectionResult, intersectionResult
    real(defReal), dimension(3)                           :: rCurrent
    class(node), pointer                                  :: genericLeaf
    type(octreeNode), pointer                             :: leafPtr
    type(topologicalObjectBox), dimension(:), allocatable :: containedObjects
    integer(shortInt)                                     :: i
    type(face), pointer                                   :: facePtr
    real(defReal)                                         :: dTravelled, dMax, dMin
    character(*), parameter                               :: here = 'findFirstIntersectedObject (octreeNode_class.f90)'

    ! First check if a leaf node already contains the particle.
    dMax = data % dMax
    call self % findLeaf(data % r, data % u, genericLeaf, checkContainment = .true.)
    if (associated(genericLeaf)) then
      ! Downcast genericLeaf to correct type then set the initial traversal position as the particle's current location.
      select type(ptr => genericLeaf)
        type is(octreeNode)
          rCurrent = data % r
          dTravelled = ZERO

        class default
          ! Call fatalError if the leaf is not an octree node.
          call fatalError(here, 'Invalid leaf node type.')

      end select

    else
      ! Compute the intersection between the ray and the root node's bounding box.
      boundingBoxPtr => self % getBoundingBoxPtr()
      call boundingBoxPtr % intersects(newIntersectionTestPayload(data % r, data % u, data % dMax), boundingBoxIntersectionResult)
      if (.not. boundingBoxIntersectionResult % intersects) return

      ! At this point, the ray intersects the bounding box. Nudge intersection coordinates slightly
      ! then find leaf node containing them.
      dTravelled = boundingBoxIntersectionResult % d
      dMax = data % dMax - dTravelled

      ! If distance to bounding box intersection is already greater than maximum allowed distance return early,
      ! else find the leaf node containing the intersection point.
      if (dMax <= ZERO) return
      rCurrent = data % r + dTravelled * data % u
      call self % findLeaf(rCurrent, data % u, genericLeaf, checkContainment = .true.)

    end if

    do
      ! If there is no leaf node containing the intersection point, the intersection has been nudged outside the root bounding box
      ! and we simply return.
      if (.not. associated(genericLeaf)) return

      ! Downcast genericLeaf to correct type.
      select type(ptr => genericLeaf)
        type is(octreeNode)
          leafPtr => ptr

        class default
          call fatalError(here, 'Invalid node type.')

      end select

      ! Call fatalError if leaf is unchecked here.
      if (leafPtr % isUnchecked) call fatalError(here, 'Unchecked leaf.')

      ! Retrieve the leaf's bounding box and compute the distance to exit.
      boundingBoxPtr => leafPtr % getBoundingBoxPtr()
      call boundingBoxPtr % intersects(newIntersectionTestPayload(rCurrent, data % u, dMax), boundingBoxIntersectionResult)

      ! Retrieve all contained objects within the leaf and test them all for an intersection.
      dTravelled = dTravelled + boundingBoxIntersectionResult % d
      dMax = data % dMax - dTravelled
      dMin = dTravelled

      if (leafPtr % isIntersecting) then
        containedObjects = leafPtr % getContainedObjects()
        do i = 1, size(containedObjects)
          select type(ptr => containedObjects(i) % ptr)
            type is(face)
              facePtr => ptr

            class default
              call fatalError(here, 'Object with index: '//numToChar(ptr % getIdx())//' is not a face.')

          end select

          ! Skip test if current face is not a boundary face.
          if (.not. facePtr % getIsBoundary()) cycle
          call facePtr % intersects_Ray(newIntersectionTestPayload(data % r, data % u, data % dMax), intersectionResult)
          if (intersectionResult % intersects .and. intersectionResult % d < dMin) then
            dMin = intersectionResult % d
            data % d = dMin
            firstIntersectedObject % ptr => facePtr

          end if

        end do

      end if

      if (associated(firstIntersectedObject % ptr)) return
      if (dMax <= ZERO) return
      rCurrent = data % r + dTravelled * data % u
      call leafPtr % findLeaf(rCurrent, data % u, genericLeaf, .true., .true.)

    end do

  end subroutine findFirstIntersectedObject

  !!
  !!
  !!
  recursive subroutine findNearestObject(self, r, radiusSquared, nearestObject)
    class(octreeNode), intent(in)             :: self
    real(defReal), dimension(3), intent(in)   :: r
    real(defReal), intent(inout)              :: radiusSquared
    type(topologicalObjectBox), intent(inout) :: nearestObject
    character(*), parameter                   :: here = 'findNearestObject (octreeNode_class.f90)'

    call fatalError(here, 'Octrees do not support nearest neighbour searches.')

  end subroutine findNearestObject

  !!
  !!
  !!
  elemental function getChildrenNumber(self) result(nChildren)
    class(octreeNode), intent(in) :: self
    integer(shortInt)             :: nChildren

    nChildren = 8

  end function getChildrenNumber

  !!
  !!
  !!
  function getDescentChildIdx(self, r, u) result(childIdx)
    class(octreeNode), intent(in)           :: self
    real(defReal), dimension(3), intent(in) :: r, u
    integer(shortInt)                       :: childIdx
    real(defReal), dimension(3)             :: boundingBoxCentre
    integer(shortInt)                       :: i
    logical(defBool)                        :: incrementIdx

    boundingBoxCentre = self % getBoundingBoxCentre()
    childIdx = 1
    do i = 1, 3
      incrementIdx = (areEqual(r(i), boundingBoxCentre(i)) .and. ZERO < u(i)) .or. &
                     ((.not. areEqual(r(i), boundingBoxCentre(i))) .and. boundingBoxCentre(i) < r(i))
      childIdx = merge(childIdx + 2 ** (i - 1), childIdx, incrementIdx)

    end do

  end function getDescentChildIdx

  !!
  !!
  !!
  elemental function getIsInside(self) result(isInside)
    class(octreeNode), intent(in) :: self
    logical(defBool)              :: isInside

    isInside = self % isInside

  end function getIsInside

  !!
  !!
  !!
  elemental function getIsIntersecting(self) result(isIntersecting)
    class(octreeNode), intent(in) :: self
    logical(defBool)              :: isIntersecting

    isIntersecting = self % isIntersecting

  end function getIsIntersecting

  !!
  !!
  !!
  elemental function getIsOutside(self) result(isOutside)
    class(octreeNode), intent(in) :: self
    logical(defBool)              :: isOutside

    isOutside = self % isOutside

  end function getIsOutside

  !!
  !!
  !!
  elemental function getIsUnchecked(self) result(isUnchecked)
    class(octreeNode), intent(in) :: self
    logical(defBool)              :: isUnchecked

    isUnchecked = self % isUnchecked

  end function getIsUnchecked

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  pure recursive subroutine kill(self)
    class(octreeNode), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % isUnchecked = .true.
    self % isInside = .false.
    self % isOutside = .false.
    self % isIntersecting = .false.

  end subroutine kill

  !!
  !!
  !!
  subroutine preparePayloadForChild(self, childNumber, payload)
    class(octreeNode), target, intent(in)  :: self
    integer(shortInt), intent(in)          :: childNumber
    class(buildNodePayload), intent(inout) :: payload
    type(buildOctreeNodePayload), pointer  :: payloadPtr
    real(defReal), dimension(3, 2)         :: boundingBoxBounds
    real(defReal), dimension(3)            :: boundingBoxCentre
    integer(shortInt)                      :: i
    character(*), parameter                :: here = 'preparePayloadForChild (octreeNode_class.f90)'

    select type(ptr => payload)
      type is(buildOctreeNodePayload)
        payloadPtr => ptr

      class default
        call fatalError(here, 'Invalid parent payload type.')

    end select

    ! Calculate the child's bounding box bounds.
    boundingBoxBounds = self % getBoundingBoxBounds()
    boundingBoxCentre = self % getBoundingBoxCentre()
    do i = 1, 3
      if (btest(childNumber - 1, i - 1)) then
        payloadPtr % bounds(i, 1) = boundingBoxCentre(i)
        payloadPtr % bounds(i, 2) = boundingBoxBounds(i, 2)

      else
        payloadPtr % bounds(i, 1) = boundingBoxBounds(i, 1)
        payloadPtr % bounds(i, 2) = boundingBoxCentre(i)

      end if

    end do

  end subroutine preparePayloadForChild

end module octreeNode_class