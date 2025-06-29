module objectNode_class

  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use faceShelf_class,              only : faceShelf
  use genericProcedures,            only : append
  use kdTreeNode_inter,             only : kdTreeNode, kill_super => kill
  use numPrecision
  use vertexShelf_class,            only : vertexShelf

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(kdTreeNode) :: objectNode
    type(objectNode), pointer       :: left => null(), right => null()
  contains
    ! Build procedures.
    procedure          :: kill
    procedure          :: init
    ! Runtime procedures.
    generic            :: findPotentiallyIntersectedObjects => findPotentiallyIntersectedObjects_BoundingBox
    procedure, private :: findPotentiallyIntersectedObjects_BoundingBox
    procedure          :: process
    procedure          :: search
  end type objectNode

contains
  !!
  !!
  !!
  recursive subroutine findPotentiallyIntersectedObjects_BoundingBox(self, boundingBox, objectIdxs)
    class(objectNode), intent(in)                               :: self
    type(axisAlignedBoundingBox), intent(in)                    :: boundingBox
    integer(shortInt), dimension(:), allocatable, intent(inout) :: objectIdxs
    type(axisAlignedBoundingBox)                                :: nodeBoundingBox

    ! If current node is a leaf, check for intersection against each face in the node.
    if (self % getIsLeaf()) then
      call append(objectIdxs, self % getDataIdxs())
      return

    end if
    
    ! Check if the bounding box intersects the current node's bounding box.
    nodeBoundingBox = self % getBoundingBox()
    if (.not. nodeBoundingBox % intersects(boundingBox)) return

    ! If the bounding box intersects the bounding box of the current node, descend deeper into the tree.
    call self % left % findPotentiallyIntersectedObjects_BoundingBox(boundingBox, objectIdxs)
    call self % right % findPotentiallyIntersectedObjects_BoundingBox(boundingBox, objectIdxs)

  end subroutine findPotentiallyIntersectedObjects_BoundingBox

  !!
  !!
  !!
  pure recursive subroutine kill(self)
    class(objectNode), intent(inout) :: self

    ! Local.
    if (associated(self % left)) then
      call self % left % kill()
      deallocate(self % left)
      nullify(self % left)

    end if

    if (associated(self % right)) then
      call self % right % kill()
      deallocate(self % right)
      nullify(self % right)

    end if

    ! Superclass.
    call kill_super(self)

  end subroutine kill

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
  recursive subroutine init(self, data, idxs, lowerBound, upperBound, nNodes, nLeaves, boundingBoxes, parentIdx)
    class(objectNode), intent(inout)                       :: self
    real(defReal), dimension(:, :), intent(in)             :: data
    integer(shortInt), dimension(:), intent(inout)         :: idxs
    integer(shortInt), intent(in)                          :: lowerBound, upperBound
    integer(shortInt), intent(inout)                       :: nNodes, nLeaves
    type(axisAlignedBoundingBox), dimension(:), intent(in) :: boundingBoxes
    integer(shortInt), intent(in), optional                :: parentIdx
    type(axisAlignedBoundingBox)                           :: boundingBox
    integer(shortInt)                                      :: cutIdx
    
    ! Initialise superclass.
    call self % setupBase(data, idxs, lowerBound, upperBound, nNodes, nLeaves, parentIdx)

    ! If node has been identified as a leaf, compute its bounding box here.
    if (self % getIsLeaf()) then
      boundingBox = self % getBoundingBox()
      call boundingBox % computeBounds(boundingBoxes(idxs(self % getLowerBound():self % getUpperBound())))
      call self % setBoundingBox(boundingBox)
      return

    end if
    
    ! Build new children nodes. Catch degenerate cases for which there are no vertices on the left
    ! or right, in which case only a single child node is built. Then, recompute the current 
    ! node's bounding box exactly from its children nodes. If one such children is missing, then 
    ! the current node's bounding box is simply that of its only child node, and the cut value is 
    ! set accordingly. If the two children nodes are present, then the current node's bounding box 
    ! is taken as the average of its two children's bounding boxes.
    cutIdx = self % getCutIdx()
    call self % setHasLeft()
    allocate(self % left)
    call self % left % init(data, idxs, lowerBound, cutIdx, nNodes, nLeaves, boundingBoxes, self % getIdx())
    
    call self % setHasRight()
    allocate(self % right)
    call self % right % init(data, idxs, cutIdx + 1, upperBound, nNodes, nLeaves, boundingBoxes, self % getIdx())

    ! Update bounding box from children bounding boxes.
    boundingBox = self % getBoundingBox()
    call boundingBox % computeBounds([self % left % getBoundingBox(), self % right % getBoundingBox()])
    call self % setBoundingBox(boundingBox)

  end subroutine init

  !! Subroutine 'process'
  !!
  !! Basic description:
  !!   Returns the index of the vertex in the terminal node which is closest to the supplied 3-D
  !!   coordinates.
  !!
  !! Arguments:
  !!   treeData [in]    -> 3-D coordinates of all the vertices in the tree.
  !!   coordinates [in] -> Supplied 3-D coordinates.
  !!   ballsize [inout] -> Smallest search radius reached up to this point.
  !!   idx [out]        -> Index of the vertex closest to the supplied 3-D coordinates.
  !!
  subroutine process(self, r, radiusSquared, idxs, faces, idx)
    class(objectNode), intent(in)               :: self
    real(defReal), dimension(3), intent(in)     :: r
    real(defReal), intent(inout)                :: radiusSquared
    integer(shortInt), dimension(:), intent(in) :: idxs
    type(faceShelf), intent(in)                 :: faces
    integer(shortInt), intent(inout)            :: idx
    integer(shortInt)                           :: i
    real(defReal)                               :: dSquared
    
    ! Loop over all vertices in the terminal node.
    do i = self % getLowerBound(), self % getUpperBound()
      ! Compute the distance squared to the current face and update minimum distance.
      dSquared = faces % distanceSquaredFromFace(idxs(i), r)

      if (dSquared < radiusSquared) then
        ! Set idx to the index of the vertex corresponding to the current lowest distance and ballSize
        ! to said lowest distance.
        idx = i
        radiusSquared = dSquared

      end if

    end do

  end subroutine process
  
  !! Subroutine 'search'
  !!
  !! Basic description:
  !!   Searches a node and its children node for the vertex which is closest to the supplied 3-D
  !!   coordinates.
  !!
  !! Detailed description:
  !!   First checks whether the current node is a terminal node. If it is, the node is simply
  !!   processed. If not, the subroutine checks which child node of the current node is closest to
  !!   the supplied coordinates. If it exists, the closer node is then searched. In some cases,
  !!   however, the further node may actually contain a vertex which is closer than any vertices in
  !!   the closer node. In this case, the further node is searched only provided that the square of
  !!   the distance between the coordinates and the further node's bounding box is lower than the
  !!   minimum ball size already returned when searching the closer node.
  !!
  !! Arguments:
  !!   treeData [in]    -> 3-D coordinates of all the vertices in the tree.
  !!   coordinates [in] -> Supplied 3-D coordinates.
  !!   ballsize [inout] -> Smallest search radius reached up to this point.
  !!   idx [out]        -> Index of the vertex closest to the supplied 3-D coordinates. 
  !!
  recursive subroutine search(self, r, radiusSquared, idxs, vertices, faces, idx)
    class(objectNode), intent(in)               :: self
    real(defReal), dimension(3), intent(in)     :: r
    real(defReal), intent(inout)                :: radiusSquared
    integer(shortInt), dimension(:), intent(in) :: idxs
    type(vertexShelf), intent(in)               :: vertices
    type(faceShelf), intent(in)                 :: faces
    integer(shortInt), intent(out)              :: idx
    class(objectNode), pointer                  :: nearNode, farNode
    type(axisAlignedBoundingBox)                :: boundingBox
    
    ! If the current node is a leaf simply process it.
    if (self % getIsLeaf()) then
      call self % process(r, radiusSquared, idxs, faces, idx)
      return

    end if

    ! Determine which node is near and which is far based on the current node's cut value.
    if (r(self % getCutDimension()) < self % getCutValue()) then
      nearNode => self % left
      farNode => self % right

    else
      nearNode => self % right
      farNode => self % left

    end if

    ! Always search the nearer node first.
    if (associated(nearNode)) call nearNode % search(r, radiusSquared, idxs, vertices, faces, idx)

    ! Search the further node only if the distance to its bounding box is less than the current
    ! best distance.
    if (associated(farNode)) then
      boundingBox = farNode % getBoundingBox()
      if (boundingBox % distanceSquared(r) < radiusSquared) call farNode % search(r, radiusSquared, idxs, vertices, faces, idx)

    end if

  end subroutine search

end module objectNode_class