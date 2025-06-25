module kdTreeNode_inter

  use node_inter,         only : node, kill_super => kill
  use genericProcedures,  only : quickSort
  use numPrecision,
  use universalVariables, only : HALF, INF, ZERO

  implicit none
  private

  ! Extendable procedures.
  public :: kill

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
  type, public, abstract, extends(node) :: kdTreeNode
    private
    integer(shortInt)                   :: cutDimension = 0, cutIdx = 0, lowerBound = 0, upperBound = 0
    logical(defBool)                    :: hasLeft = .false., hasRight = .false.
    real(defReal)                       :: cutValue = ZERO
  contains
    ! Build procedures.
    procedure                    :: kill
    procedure                    :: setupBase
    ! Runtime procedures.
    procedure                    :: getCutDimension
    procedure                    :: getCutIdx
    procedure                    :: getCutValue
    procedure                    :: getDataIdxs
    procedure                    :: getHasLeft
    procedure                    :: getHasRight
    procedure                    :: getLowerBound
    procedure                    :: getUpperBound
    procedure                    :: setHasLeft
    procedure                    :: setHasRight
  end type kdTreeNode

contains
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
  elemental function getHasLeft(self) result(hasLeft)
    class(kdTreeNode), intent(in) :: self
    logical(defBool)              :: hasLeft

    hasLeft = self % hasLeft

  end function getHasLeft

  !!
  !!
  !!
  elemental function getHasRight(self) result(hasRight)
    class(kdTreeNode), intent(in) :: self
    logical(defBool)              :: hasRight

    hasRight = self % hasRight

  end function getHasRight

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
    self % hasLeft = .false.
    self % hasRight = .false.

  end subroutine kill

  !!
  !!
  !!
  elemental subroutine setHasLeft(self)
    class(kdTreeNode), intent(inout) :: self

    self % hasLeft = .true.

  end subroutine setHasLeft

  !!
  !!
  !!
  elemental subroutine setHasRight(self)
    class(kdTreeNode), intent(inout) :: self

    self % hasRight = .true.

  end subroutine setHasRight
  
  !! Subroutine 'init'
  !!
  !! Basic description:
  !!   Recursively initialises the node and its children nodes.
  !!
  !! Detailed description:
  !!   If the node is a leaf, the subroutine simply computes its bounding box and sets its lower 
  !!   and upper bounds. Else, it computes an initial, approximate bounding box for the node from 
  !!   that of its parent node (if present). Then, the cut dimension for the node is computed from 
  !!   the coordinate with the greatest extent amongst the node's vertices; the cut value is then 
  !!   the average coordinate from all of the node's vertices along the cut dimension. From this cut 
  !!   value, the vertex to use to split the current node into children nodes is computed; finally, 
  !!   the present node's bounding box is recomputed exactly from that of its children nodes.
  !!
  !! Arguments:
  !!   coords [in]          -> Coordinates of all the vertices in the tree.
  !!   verticesIdxs [inout] -> Internal re-ordering of the vertices indices in the tree.
  !!   lowerBound [in]      -> Vertex of smallest index in the node.
  !!   upperBound [in]      -> Vertex of greatest index in the node.
  !!   parent [in]          -> Parent node.
  !!
  subroutine setupBase(self, data, idxs, lowerBound, upperBound, nNodes, nLeaves, parentIdx)
    class(kdTreeNode), intent(inout)                                 :: self
    real(defReal), dimension(:, :), intent(in)                       :: data
    integer(shortInt), dimension(:), intent(inout)                   :: idxs
    integer(shortInt), intent(in)                                    :: lowerBound, upperBound
    integer(shortInt), intent(inout)                                 :: nNodes, nLeaves
    integer(shortInt), intent(in), optional                          :: parentIdx
    integer(shortInt)                                                :: i, j, cutDimension, cutDatumIdx, nData, middleIdx
    integer(shortInt), dimension(:), allocatable                     :: indicesArray
    integer(shortInt), dimension(upperBound - lowerBound + 1)        :: nodeDataIdxs, sortedNodeDataIdxs
    real(defReal)                                                    :: cutValue, mean, variance, maxVariance, difference
    logical(defBool)                                                 :: isChild
    real(defReal), dimension(upperBound - lowerBound + 1)            :: sortedCoords
    real(defReal), dimension(6)                                      :: bounds

    ! Initialise superclass.
    nNodes = nNodes + 1
    call self % setIdx(nNodes)
    if (present(parentIdx)) call self % setParentIdx(parentIdx)
    
    ! Set the node's lower and upper bounds, and compute the number of data points.
    self % lowerBound = lowerBound
    self % upperBound = upperBound
    nData = upperBound - lowerBound + 1
    nodeDataIdxs = idxs(lowerBound:upperBound)
    
    ! If nData <= bucketSize, the node is a leaf and there is no need to further subdivide. 
    ! Simply compute the node's bounding box along each dimension and return.
    if (nData <= self % getBucketSize()) then
      nLeaves = nLeaves + 1
      call self % setIsLeaf()
      return
    
    end if

    ! Initialise maxVariance = -INF, determine if the current node is a child node, then loop over
    ! all dimensions.
    maxVariance = -INF
    do i = 1, 3
      ! If the parent node is allocated and the current dimension is not equal to its cut 
      ! dimension, then the bounding box for the child node is set to its parent's bounding box. 

      ! Compute the mean value along the current dimension and initialise variance = ZERO.
      mean = sum(data(i, nodeDataIdxs)) / nData
      variance = ZERO
      ! Loop over all vertices and update the variance along the current dimension.
      do j = 1, nData
        difference = data(i, nodeDataIdxs(j)) - mean
        variance = variance + difference * difference

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
    sortedCoords = data(cutDimension, nodeDataIdxs)
    sortedNodeDataIdxs = nodeDataIdxs
    call quickSort(sortedCoords, sortedNodeDataIdxs)
    if (mod(nData, 2) == 0) then
      middleIdx = nData / 2
      cutValue = HALF * (sortedCoords(middleIdx) + sortedCoords(middleIdx + 1))

    else
      middleIdx = nData / 2 + 1
      cutValue = sortedCoords(middleIdx)

    end if
    self % cutValue = cutValue
    idxs(lowerBound:upperBound) = sortedNodeDataIdxs
    self % cutIdx = lowerBound + middleIdx - 1

  end subroutine setupBase

end module kdTreeNode_inter