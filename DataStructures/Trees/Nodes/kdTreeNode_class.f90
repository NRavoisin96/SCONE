module kdTreeNode_class

  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use genericProcedures,            only : fatalError, numToChar, quickSort
  use node_inter,                   only : buildNodePayload, kill_super => kill, node, nodeBox
  use numPrecision,
  use topologicalObject_inter,      only : topologicalObjectBox
  use universalVariables,           only : HALF, INF, ZERO

  implicit none
  private

  ! Extendable procedures.
  public :: kill

  !!
  !!
  !!
  type, public, extends(buildNodePayload)        :: buildKDTreeNodePayload
    integer(shortInt)                            :: lowerBound = 0, upperBound = 0
    integer(shortInt), dimension(:), allocatable :: idxs
  end type buildKDTreeNodePayload

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
  type, public, extends(node) :: kdTreeNode
    private
    integer(shortInt)         :: cutDimension = 0, cutIdx = 0, lowerBound = 0, upperBound = 0
    real(defReal)             :: cutValue = ZERO
  contains
    ! Build procedures.
    procedure :: allocateChild
    procedure :: build
    procedure :: getChildrenNumber
    procedure :: kill
    procedure :: preparePayloadForChild
    ! Runtime procedures.
    procedure :: findNearestObject
    procedure :: getCutDimension
    procedure :: getCutIdx
    procedure :: getCutValue
    procedure :: getDataIdxs
    procedure :: getDescentChildIdx
    procedure :: getLowerBound
    procedure :: getUpperBound
    procedure :: process
  end type kdTreeNode

contains
  !!
  !!
  !!
  subroutine allocateChild(self, ptr)
    class(kdTreeNode), intent(in)     :: self
    class(node), pointer, intent(out) :: ptr

    allocate(kdTreeNode :: ptr)

  end subroutine allocateChild

  !!
  !!
  !!
  subroutine build(self, payload, stop)
    class(kdTreeNode), intent(inout)             :: self
    class(buildNodePayload), intent(inout)       :: payload
    logical(defBool), intent(out)                :: stop
    type(buildKDTreeNodePayload), pointer        :: payloadPtr
    integer(shortInt)                            :: cutDimension, i, j, middleIdx, nObjects
    integer(shortInt), dimension(:), allocatable :: objectIdxs, sortedObjectIdxs
    real(defReal)                                :: cutValue, diff, maxVariance, mean, variance
    real(defReal), dimension(:, :), allocatable  :: centroids
    real(defReal), dimension(:), allocatable     :: sortedCentroids
    character(*), parameter                      :: here = 'build (kdTreeNode_class.f90)'

    ! Initialise stop = .false.
    stop = .false.
    
    ! Downgrade payload type.
    select type(ptr => payload)
      type is (buildKDTreeNodePayload)
        payloadPtr => ptr

      class default
        call fatalError(here, 'Invalid payload type.')

    end select

    ! Set the node's lower and upper bounds, and compute the number of data points.
    self % lowerBound = payloadPtr % lowerBound
    self % upperBound = payloadPtr % upperBound
    nObjects = self % upperBound - self % lowerBound + 1
    objectIdxs = payloadPtr % idxs(self % lowerBound:self % upperBound)
    
    ! If nData <= bucketSize, the node is a leaf and there is no need to further subdivide. 
    ! Simply compute the node's bounding box along each dimension and return.
    if (self % getDepth() == payloadPtr % maxDepth .or. nObjects <= self % getBucketSize()) then
      call self % addContainedObject(payloadPtr % shelf % getObjectBox(objectIdxs))
      stop = .true.
      return
    
    end if

    allocate(centroids(3, nObjects))
    allocate(sortedCentroids(nObjects))
    centroids = payloadPtr % shelf % getObjectCentroid(objectIdxs)

    ! Initialise maxVariance = -INF, determine if the current node is a child node, then loop over
    ! all dimensions.
    maxVariance = -INF
    do i = 1, 3
      ! If the parent node is allocated and the current dimension is not equal to its cut 
      ! dimension, then the bounding box for the child node is set to its parent's bounding box. 

      ! Compute the mean value along the current dimension and initialise variance = ZERO.
      mean = sum(centroids(i, :)) / nObjects
      variance = ZERO
      ! Loop over all vertices and update the variance along the current dimension.
      do j = 1, nObjects
        diff = centroids(i, j) - mean
        variance = variance + diff * diff

      end do

      ! If variance > maxVariance, update cutDimension and maxVariance.
      if (variance > maxVariance) then
        cutDimension = i
        maxVariance = variance

      end if

    end do
    ! Set the node's cut dimension.
    self % cutDimension = cutDimension
    
    ! Compute the cut value by averaging the coordinates of all the vertices in the node along the
    ! cut dimension, then set the node's cut value.

    sortedCentroids = centroids(self % cutDimension, :)
    sortedObjectIdxs = objectIdxs
    call quickSort(sortedCentroids, sortedObjectIdxs)
    if (mod(nObjects, 2) == 0) then
      middleIdx = nObjects / 2
      cutValue = HALF * sum(sortedCentroids([middleIdx, middleIdx + 1]))

    else
      middleIdx = nObjects / 2 + 1
      cutValue = sortedCentroids(middleIdx)

    end if
    self % cutValue = cutValue
    payloadPtr % idxs(self % lowerBound:self % upperBound) = sortedObjectIdxs
    self % cutIdx = self % lowerBound + middleIdx - 1

  end subroutine build

  !!
  !!
  !!
  recursive subroutine findNearestObject(self, r, radiusSquared, idx)
    class(kdTreeNode), intent(in)            :: self
    real(defReal), dimension(3), intent(in)  :: r
    real(defReal), intent(inout)             :: radiusSquared
    integer(shortInt), intent(inout)         :: idx
    type(nodeBox), dimension(2)              :: children
    integer(shortInt)                        :: nChildren
    class(node), pointer                     :: genericFar, genericNear
    class(kdTreeNode), pointer               :: farNode, nearNode
    character(*), parameter                  :: here = 'search (kdTreeNode_class.f90)'

    ! If the current node is a leaf simply process it.
    if (self % getIsLeaf()) then
      call self % process(r, radiusSquared, idx)
      return

    end if

    ! Determine which node is near and which is far based on the current node's cut value.
    children = self % getChildren()
    nChildren = size(children)
    if (nChildren /= 2) call fatalError(here, 'Invalid number of children for k-d tree node: '//numToChar(nChildren)//'.')
    if (r(self % cutDimension) < self % cutValue) then
      genericNear => children(1) % ptr
      genericFar => children(2) % ptr

    else
      genericNear => children(2) % ptr
      genericFar => children(1) % ptr

    end if

    ! Downcast pointers.
    if (associated(genericNear)) then
      select type(ptr => genericNear)
        type is(kdTreeNode)
          nearNode => ptr

        class default
          call fatalError(here, 'Invalid k-d tree node type.')

      end select

    else
      call fatalError(here, 'Unable to retrieve child k-d tree node.')

    end if

    if (associated(genericFar)) then
      select type(ptr => genericFar)
        type is(kdTreeNode)
          farNode => ptr

        class default
          call fatalError(here, 'Invalid k-d tree node type.')

      end select

    else
      call fatalError(here, 'Unable to retrieve child k-d tree node.')

    end if

    ! Always search the nearer node first.
    if (associated(nearNode)) call nearNode % findNearestObject(r, radiusSquared, idx)

    ! Search the further node only if the distance to its bounding box is less than the current
    ! best distance.
    if (associated(farNode)) then
      if (farNode % distanceSquared(r) < radiusSquared) call farNode % findNearestObject(r, radiusSquared, idx)

    end if

  end subroutine findNearestObject

  !!
  !!
  !!
  elemental function getChildrenNumber(self) result(nChildren)
    class(kdTreeNode), intent(in) :: self
    integer(shortInt)             :: nChildren

    nChildren = 2

  end function getChildrenNumber

  !!
  !!
  !!
  elemental function getCutDimension(self) result(cutDimension)
    class(kdTreeNode), intent(in) :: self
    integer(shortInt)             :: cutDimension

    cutDimension = self % cutDimension

  end function getCutDimension
  
  !!
  !!
  !!
  elemental function getCutIdx(self) result(cutIdx)
    class(kdTreeNode), intent(in) :: self
    integer(shortInt)             :: cutIdx

    cutIdx = self % cutIdx

  end function getCutIdx

  !!
  !!
  !!
  elemental function getCutValue(self) result(cutValue)
    class(kdTreeNode), intent(in) :: self
    real(defReal)                 :: cutValue

    cutValue = self % cutValue

  end function getCutValue
  
  !! Function 'getDataIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the data points in the node.
  !!
  !! Result:
  !!   idxs -> Indices of the data points in the node. Note: these indices are internally sorted so
  !!           they will not correspond to the indices of the original data points!
  !!
  pure function getDataIdxs(self) result(idxs)
    class(kdTreeNode), intent(in)                                           :: self
    integer(shortInt), dimension(self % upperBound - self % lowerBound + 1) :: idxs
    integer(shortInt)                                                       :: i
    
    do i = 0, self % upperBound - self % lowerBound
      idxs(i + 1) = self % lowerBound + i

    end do

  end function getDataIdxs

  !!
  !!
  !!
  pure function getDescentChildIdx(self, r) result(childIdx)
    class(kdTreeNode), intent(in)           :: self
    real(defReal), dimension(3), intent(in) :: r
    integer(shortInt)                       :: childIdx

    childIdx = 1
    if (self % cutValue <= r(self % cutDimension)) childIdx = 2

  end function getDescentChildIdx

  !!
  !!
  !!
  elemental function getLowerBound(self) result(lowerBound)
    class(kdTreeNode), intent(in) :: self
    integer(shortInt)             :: lowerBound

    lowerBound = self % lowerBound

  end function getLowerBound

  !!
  !!
  !!
  elemental function getUpperBound(self) result(upperBound)
    class(kdTreeNode), intent(in) :: self
    integer(shortInt)             :: upperBound

    upperBound = self % upperBound

  end function getUpperBound

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  pure recursive subroutine kill(self)
    class(kdTreeNode), intent(inout) :: self

    ! Superclass.
    call kill_super(self)
    
    ! Local.
    self % cutDimension = 0
    self % cutIdx = 0
    self % lowerBound = 0
    self % upperBound = 0
    self % cutValue = ZERO

  end subroutine kill

  !!
  !!
  !!
  subroutine preparePayloadForChild(self, childNumber, payload)
    class(kdTreeNode), target, intent(in)  :: self
    integer(shortInt), intent(in)          :: childNumber
    class(buildNodePayload), intent(inout) :: payload
    type(buildKDTreeNodePayload), pointer  :: payloadPtr
    character(*), parameter                :: here = 'preparePayloadForChild (kdTreeNode_inter.f90)'

    ! Downcast payload to correct type.
    select type(ptr => payload)
      type is(buildKDTreeNodePayload)
        payloadPtr => ptr

      class default
        call fatalError(here, 'Invalid parent payload type.')

    end select

    if (childNumber == 1) then
      ! Left node.
      payloadPtr % lowerBound = self % lowerBound
      payloadPtr % upperBound = self % cutIdx

    else
      ! Right node.
      payloadPtr % lowerBound = self % cutIdx + 1
      payloadPtr % upperBound = self % upperBound

    end if

  end subroutine preparePayloadForChild

  !!
  !!
  !!
  subroutine process(self, r, radiusSquared, idx)
    class(kdTreeNode), intent(in)                         :: self
    real(defReal), dimension(3), intent(in)               :: r
    real(defReal), intent(inout)                          :: radiusSquared
    integer(shortInt), intent(inout)                      :: idx
    type(topologicalObjectBox), dimension(:), allocatable :: testObjects
    integer(shortInt)                                     :: i
    real(defReal)                                         :: dSquared

    ! Loop over all objects in the leaf.
    testObjects = self % getTestObjects()
    do i = 1, size(testObjects)
      ! Compute the distance squared to the current face and update minimum distance.
      dSquared = testObjects(i) % ptr % distanceSquared(r)

      if (dSquared < radiusSquared) then
        ! Set idx to the index of the vertex corresponding to the current lowest distance and ballSize
        ! to said lowest distance.
        idx = testObjects(i) % ptr % getIdx()
        radiusSquared = dSquared

      end if

    end do

  end subroutine process

end module kdTreeNode_class