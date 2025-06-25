module node_inter
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use coord_class,                  only : coord
  use faceShelf_class,              only : faceShelf
  use genericProcedures,            only : append, areEqual, quickSort, removeDuplicates, swap
  use numPrecision
  use universalVariables,           only : ZERO, HALF, SURF_TOL, INF
  use vertexShelf_class,            only : vertexShelf
  
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
  type, public, abstract         :: node
    private
    integer(shortInt)            :: bucketSize = 4, idx = 0, parentIdx = 0
    logical(defBool)             :: isLeaf = .false.
    type(axisAlignedBoundingBox) :: boundingBox
  contains
    ! Build procedures.
    procedure :: kill
    procedure :: setBoundingBox
    procedure :: setIdx
    procedure :: setIsLeaf
    procedure :: setParentIdx
    ! Runtime procedures.
    procedure :: getBoundingBox
    procedure :: getBucketSize
    procedure :: getIdx
    procedure :: getIsLeaf
    procedure :: getParentIdx
  end type node

contains
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  pure recursive subroutine kill(self)
    class(node), intent(inout) :: self
    
    ! Local.
    self % idx = 0
    self % parentIdx = 0
    self % isLeaf = .false.
    call self % boundingBox % kill()

  end subroutine kill

  !! Function 'getBoundingBox'
  !!
  !! Basic description:
  !!   Returns the bounding box of the node.
  !!
  !! Result:
  !!   boundingBox -> Bounding box of the node.
  !!
  pure function getBoundingBox(self) result(boundingBox)
    class(node), intent(in)      :: self
    type(axisAlignedBoundingBox) :: boundingBox

    boundingBox = self % boundingBox

  end function getBoundingBox

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
  elemental function getParentIdx(self) result(parentIdx)
    class(node), intent(in) :: self
    integer(shortInt)             :: parentIdx

    parentIdx = self % parentIdx

  end function getParentIdx

  !!
  !!
  !!
  elemental subroutine setBoundingBox(self, boundingBox)
    class(node), intent(inout)               :: self
    type(axisAlignedBoundingBox), intent(in) :: boundingBox

    self % boundingBox = boundingBox

  end subroutine setBoundingBox

  !!
  !!
  !!
  elemental subroutine setIdx(self, idx)
    class(node), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx

    self % idx = idx

  end subroutine setIdx

  !!
  !!
  !!
  elemental subroutine setIsLeaf(self)
    class(node), intent(inout) :: self

    self % isLeaf = .true.

  end subroutine setIsLeaf

  !!
  !!
  !!
  elemental subroutine setParentIdx(self, parentIdx)
    class(node), intent(inout)    :: self
    integer(shortInt), intent(in) :: parentIdx

    self % parentIdx = parentIdx

  end subroutine setParentIdx
  
end module node_inter