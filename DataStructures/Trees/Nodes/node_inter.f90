module node_inter
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use element_class,                only : elementBox
  use genericProcedures,            only : append, areEqual, fatalError, numToChar, quickSort, removeDuplicates, swap
  use numPrecision
  use publicObjects,                only : coordData
  use topologicalObject_inter,      only : topologicalObjectBox
  use topologicalObjectShelf_class, only : topologicalObjectShelf
  use universalVariables,           only : ZERO, HALF, SURF_TOL, INF
  
  implicit none
  private

  ! Extendable procedures.
  public :: kill

  !!
  !!
  !!
  type, public :: buildNodePayload
    type(topologicalObjectShelf), pointer :: shelf => null()
    integer(shortInt)                     :: bucketSize = 0, depth = 0, maxDepth = 0, nLeaves = 0, nNodes = 0
    logical(defBool)                      :: computeLeafBoundingBox = .true., updateParentBoundingBox = .true.
    class(node), pointer                  :: parent => null()
  end type buildNodePayload

  !!
  !!
  !!
  type, public :: nodeBox
    class(node), pointer :: ptr => null()
  end type nodeBox
  
  !!
  !! Node of a kd-tree. Used to recursively subdivide space into smaller and smaller numbers of
  !! vertices. Each node subdivides space along an axis-aligned dimension known as the cut
  !! dimension, which is taken as the dimension of greatest extent among the vertices in the node.
  !! The node is then split into (generally) two children nodes along the cut dimension. The value
  !! used to split the node along the cut dimension is taken as the average of the vertices'
  !! coordinates along the cut dimension. The child node containing vertices whose coordinates along
  !! the cut dimension less than the cut value is termed the left node and vice versa.
  !!
  !! Private members:
  !!   cutDimension   -> Index of the dimension (1 = x, 2 = y, 3 = z) used to split the node into
  !!                     children nodes.
  !!   lowerBound     -> Lowest-index vertex in the node.
  !!   upperBound     -> Greatest-index vertex in the node.
  !!   bucketSize     -> Maximum number of vertices in a terminal node.
  !!   hasLeft        -> Does the node have a left child node?
  !!   hasRight       -> Does the node have a right child node?
  !!   cutValue       -> Value along the cutDimension used to split the node into children nodes.
  !!   cutValue_left  -> CutValue of the node's left child node.
  !!   cutValue_right -> CutValue of the node's right child node.
  !!   left           -> Left child node.
  !!   right          -> Right child node.
  !!   boundingBox    -> Axis-aligned bounding box (AABB) of the node.
  !!
  type, public, abstract                                  :: node
    private
    integer(shortInt)                                     :: bucketSize = 0, depth = 0, idx = 0
    logical(defBool)                                      :: isLeaf = .false.
    type(nodeBox), dimension(:), allocatable              :: children
    class(node), pointer                                  :: parent => null()
    type(topologicalObjectBox), dimension(:), allocatable :: containedObjects, containingObjects
    type(axisAlignedBoundingBox)                          :: boundingBox
  contains
    ! Build procedures.
    generic                                     :: addContainedObject => addContainedObject_box, addContainedObject_boxes
    procedure, private                          :: addContainedObject_box
    procedure, private                          :: addContainedObject_boxes
    generic                                     :: addContainingObject => addContainingObject_box, addContainingObject_boxes
    procedure, private                          :: addContainingObject_box
    procedure, private                          :: addContainingObject_boxes
    procedure(allocateChild), deferred          :: allocateChild
    procedure                                   :: assignElements
    procedure(build), deferred                  :: build
    procedure(getChildrenNumber), deferred      :: getChildrenNumber
    procedure, non_overridable                  :: init
    procedure, non_overridable                  :: initBoundingBox
    procedure                                   :: kill
    procedure(preparePayloadForChild), deferred :: preparePayloadForChild
    ! Runtime procedures.
    procedure, non_overridable                      :: distanceSquared
    procedure(findFirstIntersectedObject), deferred :: findFirstIntersectedObject
    generic                                     :: findIntersectedObjects => findIntersectedObjects_BoundingBox
    procedure, private                          :: findIntersectedObjects_BoundingBox 
    procedure, non_overridable                  :: findLeaf
    procedure(findNearestObject), deferred      :: findNearestObject
    procedure, non_overridable                  :: getBoundingBoxBounds
    procedure, non_overridable                  :: getBoundingBoxCentre
    procedure, non_overridable                  :: getBoundingBoxPtr
    procedure, non_overridable                  :: getBucketSize
    procedure, non_overridable                  :: getChildren
    procedure, non_overridable                  :: getDepth
    procedure(getDescentChildIdx), deferred     :: getDescentChildIdx
    procedure, non_overridable                  :: getIdx
    procedure, non_overridable                  :: getIsLeaf
    procedure, non_overridable                  :: getContainedObjects
    procedure, non_overridable                  :: getContainingObjects
    procedure, non_overridable                  :: getParentPtr
  end type node

  abstract interface
    !!
    !!
    !!
    subroutine allocateChild(self, ptr)
      import                            :: node
      class(node), intent(in)           :: self
      class(node), pointer, intent(out) :: ptr
    end subroutine allocateChild

    !!
    !!
    !!
    subroutine build(self, payload, stop)
      import                                 :: buildNodePayload, defBool, node
      class(node), intent(inout)             :: self
      class(buildNodePayload), intent(inout) :: payload
      logical(defBool), intent(out)          :: stop
    end subroutine build

    !!
    !!
    !!
    subroutine findFirstIntersectedObject(self, data, firstIntersectedObject)
      import                                  :: coordData, node, topologicalObjectBox
      class(node), intent(in)                 :: self
      type(coordData), intent(inout)          :: data
      type(topologicalObjectBox), intent(out) :: firstIntersectedObject
    end subroutine findFirstIntersectedObject

    !!
    !!
    !!
    recursive subroutine findNearestObject(self, r, radiusSquared, nearestObject)
      import                                    :: defReal, node, topologicalObjectBox
      class(node), intent(in)                   :: self
      real(defReal), dimension(3), intent(in)   :: r
      real(defReal), intent(inout)              :: radiusSquared
      type(topologicalObjectBox), intent(inout) :: nearestObject
    end subroutine findNearestObject

    !!
    !!
    !!
    elemental function getChildrenNumber(self) result(nChildren)
      import                  :: node, shortInt
      class(node), intent(in) :: self
      integer(shortInt)       :: nChildren
    end function getChildrenNumber

    !!
    !!
    !!
    function getDescentChildIdx(self, r, u) result(childIdx)
      import                                  :: defReal, node, shortInt
      class(node), intent(in)                 :: self
      real(defReal), dimension(3), intent(in) :: r, u
      integer(shortInt)                       :: childIdx
    end function getDescentChildIdx

    !!
    !!
    !!
    subroutine preparePayloadForChild(self, childNumber, payload)
      import                                 :: buildNodePayload, node, shortInt
      class(node), target, intent(in)        :: self
      integer(shortInt), intent(in)          :: childNumber
      class(buildNodePayload), intent(inout) :: payload
    end subroutine preparePayloadForChild

  end interface

contains
  !!
  !!
  !!
  subroutine addContainedObject_box(self, box)
    class(node), intent(inout)                            :: self
    type(topologicalObjectBox), intent(in)                :: box
    integer(shortInt)                                     :: nObjects
    type(topologicalObjectBox), dimension(:), allocatable :: tempObjects

    if (allocated(self % containedObjects)) then
      nObjects = size(self % containedObjects)
      allocate(tempObjects(nObjects + 1))
      tempObjects(1:nObjects) = self % containedObjects
      tempObjects(nObjects + 1) = box
      call move_alloc(tempObjects, self % containedObjects)

    else
      allocate(self % containedObjects(1))
      self % containedObjects(1) = box

    end if

  end subroutine addContainedObject_box

  !!
  !!
  !!
  subroutine addContainedObject_boxes(self, boxes)
    class(node), intent(inout)                            :: self
    type(topologicalObjectBox), dimension(:), intent(in)  :: boxes
    integer(shortInt)                                     :: nNewObjects, nObjects
    type(topologicalObjectBox), dimension(:), allocatable :: tempObjects

    nNewObjects = size(boxes)
    if (allocated(self % containedObjects)) then
      nObjects = size(self % containedObjects)
      allocate(tempObjects(nObjects + nNewObjects))
      tempObjects(1:nObjects) = self % containedObjects
      tempObjects(nObjects + 1:nObjects + nNewObjects) = boxes
      call move_alloc(tempObjects, self % containedObjects)

    else
      allocate(self % containedObjects(nNewObjects))
      self % containedObjects = boxes

    end if

  end subroutine addContainedObject_boxes

  !!
  !!
  !!
  subroutine addContainingObject_box(self, box)
    class(node), intent(inout)                            :: self
    type(topologicalObjectBox), intent(in)                :: box
    integer(shortInt)                                     :: nObjects
    type(topologicalObjectBox), dimension(:), allocatable :: tempObjects

    if (allocated(self % containingObjects)) then
      nObjects = size(self % containingObjects)
      allocate(tempObjects(nObjects + 1))
      tempObjects(1:nObjects) = self % containingObjects
      tempObjects(nObjects + 1) = box
      call move_alloc(tempObjects, self % containingObjects)

    else
      allocate(self % containingObjects(1))
      self % containingObjects(1) = box

    end if

  end subroutine addContainingObject_box

  !!
  !!
  !!
  subroutine addContainingObject_boxes(self, boxes)
    class(node), intent(inout)                            :: self
    type(topologicalObjectBox), dimension(:), intent(in)  :: boxes
    integer(shortInt)                                     :: nNewObjects, nObjects
    type(topologicalObjectBox), dimension(:), allocatable :: tempObjects

    nNewObjects = size(boxes)
    if (allocated(self % containingObjects)) then
      nObjects = size(self % containingObjects)
      allocate(tempObjects(nObjects + nNewObjects))
      tempObjects(1:nObjects) = self % containingObjects
      tempObjects(nObjects + 1:nObjects + nNewObjects) = boxes
      call move_alloc(tempObjects, self % containingObjects)

    else
      allocate(self % containingObjects(nNewObjects))
      self % containingObjects = boxes

    end if

  end subroutine addContainingObject_boxes

  !!
  !!
  !!
  subroutine assignElements(self, payload)
    class(node), intent(inout)          :: self
    class(buildNodePayload), intent(in) :: payload
    character(*), parameter             :: here = 'assignElements (node_inter.f90)'

    ! Call fatalError.
    call fatalError(here, 'Unsupported procedure.')

  end subroutine assignElements

  !!
  !!
  !!
  pure function distanceSquared(self, r) result(dSquared)
    class(node), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: dSquared

    dSquared = self % boundingBox % distanceSquared(r)

  end function distanceSquared

  !!
  !!
  !!
  recursive subroutine findIntersectedObjects_BoundingBox(self, boundingBox, intersectedObjects)
    class(node), intent(in)                                              :: self
    type(axisAlignedBoundingBox), intent(in)                             :: boundingBox
    type(topologicalObjectBox), dimension(:), allocatable, intent(inout) :: intersectedObjects
    logical(defBool)                                                     :: doesIt
    integer(shortInt)                                                    :: i, nIntersectedObjects
    type(topologicalObjectBox), dimension(:), allocatable                :: tempIntersectedObjects
    character(*), parameter :: here = 'findIntersectedObjects_BoundingBox (node_inter.f90)'

    ! Check if node is a leaf and if so simply test all its testObjects for an intersection.
    if (self % isLeaf) then
      if (.not. allocated(self % containedObjects)) return
      do i = 1, size(self % containedObjects)
        call self % containedObjects(i) % ptr % intersects_BoundingBox(boundingBox, doesIt)
        if (doesIt) then
          if (allocated(intersectedObjects)) then
            nIntersectedObjects = size(intersectedObjects)
            allocate(tempIntersectedObjects(nIntersectedObjects + 1))
            tempIntersectedObjects(1:nIntersectedObjects) = intersectedObjects
            tempIntersectedObjects(nIntersectedObjects + 1) = self % containedObjects(i)
            call move_alloc(tempIntersectedObjects, intersectedObjects)

          else
            allocate(intersectedObjects(1))
            intersectedObjects(1) = self % containedObjects(i)

          end if

        end if

      end do
      return

    end if

    ! Check for intersection in each of the node's children.
    if (.not. allocated(self % children)) call fatalError(here, 'Internal node has no children.')
    do i = 1, size(self % children)
      call self % children(i) % ptr % findIntersectedObjects_BoundingBox(boundingBox, intersectedObjects)

    end do

  end subroutine findIntersectedObjects_BoundingBox

  !!
  !!
  !!
  recursive subroutine findLeaf(self, r, u, leaf, checkContainment, skipFirstNode)
    class(node), intent(in), target         :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(node), intent(out), pointer       :: leaf
    logical(defBool), intent(in), optional  :: checkContainment, skipFirstNode
    logical(defBool)                        :: check, skip
    integer(shortInt)                       :: childIdx

    ! If skipFirstNode immediately call the parent.
    skip = .false.
    if (present(skipFirstNode)) skip = skipFirstNode

    if (skip) then
      if (associated(self % parent)) then
        call self % parent % findLeaf(r, u, leaf, checkContainment = .true.)

      else
        leaf => null()

      end if
      return

    end if

    ! Only perform containment check if it has been required.
    check = .false.
    if (present(checkContainment)) check = checkContainment

    if (check) then
      if (.not. self % boundingBox % contains(r, u)) then
        if (associated(self % parent)) then
          call self % parent % findLeaf(r, u, leaf, checkContainment = .true.)

        else
          leaf => null()

        end if
        return

      end if

    end if

    ! If cell is a leaf, simply associate the leaf pointer and return.
    if (self % getIsLeaf()) then
      leaf => self
      return

    end if

    ! Descend into correct child node.
    childIdx = self % getDescentChildIdx(r, u)
    call self % children(childIdx) % ptr % findLeaf(r, u, leaf)

  end subroutine findLeaf

  !!
  !!
  !!
  recursive subroutine init(self, payload)
    class(node), target, intent(inout)          :: self
    class(buildNodePayload), intent(inout)      :: payload
    logical(defBool)                            :: stop
    integer(shortInt)                           :: i, nChildren, nObjects
    real(defReal), dimension(:, :), allocatable :: boundingBoxBounds
    character(*), parameter                     :: here = 'init (node_inter.f90)'
    
    ! Check if payload has a valid topologicalObjectShelf.
    if (.not. associated(payload % shelf)) &
    call fatalError(here, 'Payload contains an unassociated topologicalObjectShelf.')

    ! Increment depth and set this node's parent if it was supplied.
    payload % depth = payload % depth + 1
    self % depth = payload % depth
    if (associated(payload % parent)) self % parent => payload % parent

    ! Build current node.
    payload % nNodes = payload % nNodes + 1
    self % idx = payload % nNodes
    self % bucketSize = payload % bucketSize
    call self % build(payload, stop)

    ! Check if termination criterion has been met.
    if (payload % depth == payload % maxDepth .or. stop) then
      ! Set the node as a leaf.
      self % isLeaf = .true.
      payload % nLeaves = payload % nLeaves + 1

      ! Compute leaf node bounding box if requested.
      if (payload % computeLeafBoundingBox) then
        if (.not. allocated(self % containedObjects)) &
        call fatalError(here, 'Unallocated topologicalObjectBox array for leaf with index: '//numToChar(self % idx)//'.')

        nObjects = size(self % containedObjects)
        allocate(boundingBoxBounds(3, 2 * nObjects))
        do i = 1, nObjects
          boundingBoxBounds(:, 2 * (i - 1) + 1:2 * i) = self % containedObjects(i) % ptr % getBoundingBoxBounds()

        end do
        call self % boundingBox % computeBounds(boundingBoxBounds)

      end if
      return

    end if

    ! Allocate the number of children depending on the specific type of the node.
    allocate(self % children(self % getChildrenNumber()))

    ! Build each child.
    nChildren = size(self % children)
    do i = 1, nChildren
      ! Allocate child to correct type.
      call self % allocateChild(self % children(i) % ptr)

      ! Reset payload % depth and payload % parent then prepare payload for next child and initialise.
      payload % depth = self % depth
      payload % parent => self
      call self % preparePayloadForChild(i, payload)
      call self % children(i) % ptr % init(payload)

    end do

    ! Update bounding box from children nodes if requested.
    if (payload % updateParentBoundingBox) then
      allocate(boundingBoxBounds(3, 2 * nChildren))
      do i = 1, nChildren
        boundingBoxBounds(:, 2 * (i - 1) + 1:2 * i) = self % children(i) % ptr % getBoundingBoxBounds()

      end do
      call self % boundingBox % computeBounds(boundingBoxBounds)

    end if

  end subroutine init

  !!
  !!
  !!
  pure subroutine initBoundingBox(self, bounds)
    class(node), intent(inout)              :: self
    real(defReal), dimension(6), intent(in) :: bounds

    call self % boundingBox % init(bounds)

  end subroutine initBoundingBox

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  pure recursive subroutine kill(self)
    class(node), intent(inout) :: self
    integer(shortInt)          :: i
    
    ! Local.
    self % bucketSize = 0
    self % idx = 0
    self % isLeaf = .false.

    if (allocated(self % children)) then
      do i = 1, size(self % children)
        call self % children(i) % ptr % kill()
        deallocate(self % children(i) % ptr)
        nullify(self % children(i) % ptr)

      end do
      deallocate(self % children)

    end if

    if (associated(self % parent)) nullify(self % parent)

    if (allocated(self % containedObjects)) then
      do i = 1, size(self % containedObjects)
        nullify(self % containedObjects(i) % ptr)

      end do
      deallocate(self % containedObjects)

    end if

    if (allocated(self % containingObjects)) then
      do i = 1, size(self % containingObjects)
        nullify(self % containingObjects(i) % ptr)

      end do
      deallocate(self % containingObjects)

    end if

    call self % boundingBox % kill()

  end subroutine kill

  !!
  !!
  !!
  pure function getBoundingBoxBounds(self) result(boundingBoxBounds)
    class(node), intent(in)        :: self
    real(defReal), dimension(3, 2) :: boundingBoxBounds

    boundingBoxBounds = self % boundingBox % getBounds()

  end function getBoundingBoxBounds

  !!
  !!
  !!
  pure function getBoundingBoxCentre(self) result(boundingBoxCentre)
    class(node), intent(in)     :: self
    real(defReal), dimension(3) :: boundingBoxCentre

    boundingBoxCentre = self % boundingBox % getCentre()

  end function getBoundingBoxCentre

  !! Function 'getBoundingBox'
  !!
  !! Basic description:
  !!   Returns the bounding box of the node.
  !!
  !! Result:
  !!   boundingBox -> Bounding box of the node.
  !!
  function getBoundingBoxPtr(self) result(ptr)
    class(node), target, intent(in)       :: self
    type(axisAlignedBoundingBox), pointer :: ptr

    ptr => self % boundingBox

  end function getBoundingBoxPtr

  !!
  !!
  !!
  elemental function getBucketSize(self) result(bucketSize)
    class(node), intent(in) :: self
    integer(shortInt)       :: bucketSize

    bucketSize = self % bucketSize

  end function getBucketSize

  !!
  !!
  !!
  function getChildren(self) result(children)
    class(node), intent(in)                  :: self
    type(nodeBox), dimension(:), allocatable :: children

    if (allocated(self % children)) then
      children = self % children

    else
      allocate(children(0))

    end if

  end function getChildren

  !!
  !!
  !!
  elemental function getDepth(self) result(depth)
    class(node), intent(in) :: self
    integer(shortInt)       :: depth

    depth = self % depth

  end function getDepth

  !!
  !!
  !!
  elemental function getIdx(self) result(idx)
    class(node), intent(in) :: self
    integer(shortInt)       :: idx

    idx = self % idx

  end function getIdx
  
  !!
  !!
  !!
  elemental function getIsLeaf(self) result(isLeaf)
    class(node), intent(in) :: self
    logical(defBool)        :: isLeaf

    isLeaf = self % isLeaf

  end function getIsLeaf

  !!
  !!
  !!
  function getContainedObjects(self) result(testObjects)
    class(node), intent(in)                               :: self
    type(topologicalObjectBox), dimension(:), allocatable :: testObjects

    if (allocated(self % containedObjects)) then
      testObjects = self % containedObjects

    else
      allocate(testObjects(0))

    end if

  end function getContainedObjects

  !!
  !!
  !!
  function getContainingObjects(self) result(containingObjects)
    class(node), intent(in)                               :: self
    type(topologicalObjectBox), dimension(:), allocatable :: containingObjects

    if (allocated(self % containingObjects)) then
      containingObjects = self % containingObjects

    else
      allocate(containingObjects(0))

    end if

  end function getContainingObjects

  !!
  !!
  !!
  function getParentPtr(self) result(parentPtr)
    class(node), intent(in) :: self
    class(node), pointer    :: parentPtr

    parentPtr => self % parent

  end function getParentPtr
  
end module node_inter