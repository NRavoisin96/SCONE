module octreeNode_class

  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use coord_class,                  only : coord
  use element_class,                only : element, elementBox, inclusionTestResult
  use face_class,                   only : face
  use genericProcedures,            only : append, areEqual, fatalError, numToChar
  use kdTree_class,                 only : kdTree
  use node_inter,                   only : buildNodePayload, kill_super => kill, node, nodeBox
  use numPrecision
  use topologicalObject_inter,      only : topologicalObjectBox
  use topologicalObjectShelf_class, only : topologicalObjectShelf
  use universalVariables,           only : HALF, INF, INSIDE_ELEMENT, NUDGE, ON_BOUNDARY_ELEMENT, OUTSIDE_ELEMENT

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
    integer(shortInt)                                     :: i, nearestFaceIdx, nElements
    type(buildOctreeNodePayload), pointer                 :: payloadPtr
    type(topologicalObjectBox), dimension(:), allocatable :: elements
    real(defReal), dimension(3)                           :: boundingBoxCentre
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
    nearestFaceIdx = payloadPtr % tree % findNearestObject(boundingBoxCentre)

    ! Retrieve the elements associated with the nearest face.
    elements = payloadPtr % shelf % getObjectElements(nearestFaceIdx)
    nElements = size(elements)
    if (nElements == 0) call fatalError(here, 'Unable to retrieve elements associated with face: '//numToChar(nearestFaceIdx)//'.')
    
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
                          //numToChar(nearestFaceIdx)//' is not an element.')


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
    integer(shortInt), dimension(:), allocatable          :: faceIdxs
    type(topologicalObjectBox), dimension(:), allocatable :: elements
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
        allocate(faceIdxs(nFaces))
        do i = 1, nFaces
          faceIdxs(i) = i

        end do
        elements = payloadPtr % shelf % getObjectElements(faceIdxs)
        do i = 1, size(elements)
          ! Downcast element to correct type.
          select type(ptr => elements(i) % ptr)
            type is(element)
              ! Do nothing.

            class default
              call fatalError(here, 'Element: '//numToChar(ptr % getIdx())//' is not an element.')

          end select

        end do
        call self % addContainedObject(payloadPtr % shelf % getShelf())
        call self % addContainingObject(elements)
        stop = .true.

      end if
      return

    end if

    ! Check the number of intersections between the current cell and the faces in the mesh by traversing the
    ! k-d tree starting from the root node.
    faceIdxs = payloadPtr % tree % findIntersectedObjects(self % getBoundingBoxPtr())
    nIntersectedFaces = size(faceIdxs)
    if (nIntersectedFaces == 0) then
      stop = .true.

    else
      self % isUnchecked = .false.
      self % isIntersecting = .true.
      if (self % getDepth() == payloadPtr % maxDepth .or. nIntersectedFaces <= self % getBucketSize()) then
        elements = payloadPtr % shelf % getObjectElements(faceIdxs)
        do i = 1, size(elements)
          select type(ptr => elements(i) % ptr)
            type is(element)
              ! Ok, do nothing.

            class default
              call fatalError(here, 'Element: '//numToChar(ptr % getIdx())//' is not an element.')

          end select

        end do
        call self % addContainedObject(payloadPtr % shelf % getObjectBox(faceIdxs))
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
  recursive subroutine findNearestObject(self, r, radiusSquared, idx)
    class(octreeNode), intent(in)           :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), intent(inout)            :: radiusSquared
    integer(shortInt), intent(inout)        :: idx
    character(*), parameter                 :: here = 'findNearestObject (octreeNode_class.f90)'

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
  function getDescentChildIdx(self, r) result(childIdx)
    class(octreeNode), intent(in)           :: self
    real(defReal), dimension(3), intent(in) :: r
    integer(shortInt)                       :: childIdx
    real(defReal), dimension(3)             :: boundingBoxCentre
    integer(shortInt)                       :: i

    boundingBoxCentre = self % getBoundingBoxCentre()
    childIdx = 1
    do i = 1, 3
      if (boundingBoxCentre(i) < r(i)) childIdx = childIdx + 2 ** (i - 1)

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