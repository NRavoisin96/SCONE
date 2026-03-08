module octree_class

  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use dictionary_class,             only : dictionary
  use edgeShelf_class,              only : edgeShelf
  use elementShelf_class,           only : elementShelf
  use errors_mod,                   only : fatalError
  use faceShelf_class,              only : faceShelf
  use genericProcedures,            only : append, areEqual, numToChar
  use objectKDTree_class,           only : objectKDTree
  use mesh_inter,                   only : mesh
  use numPrecision
  use octreeNode_class,             only : octreeNode
  use universalVariables,           only : INSIDE_ELEMENT, NUDGE
  use vertexShelf_class,            only : vertexShelf

  implicit none
  private

  type, public :: octree
    private
    integer(shortInt)                           :: depth = 0, nChildren = 0, nLeaves = 0, newChildIdx = 1, nMaxFaces = 0
    type(octreeNode)                            :: root
    type(octreeNode), dimension(:), allocatable :: children
  contains
    procedure :: findLeaf
    procedure :: getStorageSize
    procedure :: init
    procedure :: kill
    procedure :: refineNode
  end type octree

contains
  !!
  !!
  !!
  subroutine findLeaf(self, r, leaf, requiresContainmentCheck)
    class(octree), target, intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r
    type(octreeNode), pointer, intent(out)  :: leaf
    logical(defBool), intent(in), optional  :: requiresContainmentCheck
    integer(shortInt)                       :: firstChildIdx
    logical(defBool)                        :: checkContainment
    real(defReal), dimension(3)             :: boundingBoxCentre
    type(axisAlignedBoundingBox), pointer   :: currentNodeBoundingBoxPtr
    type(octreeNode), pointer               :: currentNodePtr

    ! Initialise pointers.
    leaf => null()
    currentNodePtr => self % root
    currentNodeBoundingBoxPtr => currentNodePtr % getBoundingBoxPtr()

    ! Only perform containment check if it has been required.
    checkContainment = .false.
    if(present(requiresContainmentCheck)) checkContainment = requiresContainmentCheck

    ! If containment check has been requested, check that the point is inside the root node.
    if(checkContainment .and. .not. currentNodeBoundingBoxPtr % contains(r)) return

    ! Start descending the tree from the root node.
    do
      firstChildIdx = currentNodePtr % getFirstChildIdx()
      if(firstChildIdx == 0) then
        ! We have found our leaf node so return.
        leaf => currentNodePtr
        return

      else
        ! Compute the appropriate child index then descend.
        boundingBoxCentre = currentNodeBoundingBoxPtr % getCentre()
        currentNodePtr => self % children(firstChildIdx + merge(4, 0, boundingBoxCentre(1) <= r(1)) + &
                                          merge(2, 0, boundingBoxCentre(2) <= r(2)) + merge(1, 0, boundingBoxCentre(3) <= r(3)))
        currentNodeBoundingBoxPtr => currentNodePtr % getBoundingBoxPtr()

      end if
      
    end do

  end subroutine findLeaf

  !!
  !!
  !!
  elemental function getStorageSize(self) result(storageSize)
    class(octree), intent(in) :: self
    integer(longInt)          :: storageSize
    integer(shortInt)         :: i

    storageSize = storage_size(self) / 8
    storageSize = storageSize + self % root % getStorageSize()
    if(allocated(self % children)) then
      do i = 1, size(self % children)
        storageSize = storageSize + self % children(i) % getStorageSize()

      end do

    end if

  end function getStorageSize

  !!
  !!
  !!
  subroutine init(self, dict, edges, elements, vertices, faces)
    class(octree), intent(inout)                :: self
    type(dictionary), intent(in)                :: dict
    type(edgeShelf), intent(in)                 :: edges
    type(elementShelf), intent(in)              :: elements
    type(vertexShelf), intent(in)               :: vertices
    type(faceShelf), intent(inout)              :: faces
    integer(shortInt)                           :: i, nChildren
    real(defReal), dimension(6)                 :: boundingBoxBounds
    type(objectKDTree)                          :: tree
    type(octreeNode), dimension(:), allocatable :: temp
    character(*), parameter                     :: HERE = 'init (octree_class.f90)'

    call dict % getOrDefault(self % depth, 'depth', 10)
    if(self % depth < 1) &
    call fatalError(HERE, 'Depth must be at least 1. Is: '//numToChar(self % depth)//'.')
    
    call dict % getOrDefault(self % nMaxFaces, 'nMaxFaces', 4)
    if(self % nMaxFaces < 1) &
    call fatalError(HERE, 'nMaxFaces must be at least 1. Is: '//numToChar(self % nMaxFaces)//'.')

    ! Initialise SAT cache data for faces.
    do i = 1, faces % getSize()
      call faces % computeFaceSATData(i, edges, vertices)

    end do

    ! Initialise k-d trees from the unstructured mesh faces and elements, then build the Cartesian grid's 
    ! root cell and all its children cells.
    call tree % init(faces % getAllFaceCentroids(), faces % getAllFaceBoundingBoxes())

    ! Initialise root node then refine it.
    boundingBoxBounds = tree % getRootBoundingBoxBounds() + [-NUDGE, -NUDGE, -NUDGE, NUDGE, NUDGE, NUDGE]
    call self % root % init(1, 0, boundingBoxBounds)
    call self % refineNode(0, elements, faces, tree, vertices, self % root)

    ! Now resize the shelf of children.
    nChildren = size(self % children)
    if(self % nChildren < nChildren) then
      allocate(temp(self % nChildren))
      temp = self % children(1:self % nChildren)
      call move_alloc(temp, self % children)

    end if

    ! Kill the k-d tree as it is no longer needed.
    call tree % kill()

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(octree), intent(inout) :: self
    integer(shortInt)            :: i

    ! Local.
    self % depth = 0
    self % nLeaves = 0
    self % nMaxFaces = 0
    call self % root % kill()
    if(allocated(self % children)) then
      do i = 1, size(self % children)
        call self % children(i) % kill()

      end do
      deallocate(self % children)

    end if

  end subroutine kill

  !!
  !!
  !!
  recursive subroutine refineNode(self, parentIdx, elements, faces, tree, vertices, node)
    class(octree), intent(inout)                 :: self
    integer(shortInt), intent(in)                :: parentIdx
    type(elementShelf), intent(in)               :: elements
    type(faceShelf), intent(in)                  :: faces
    type(objectKDTree), intent(in)               :: tree
    type(vertexShelf), intent(in)                :: vertices
    type(octreeNode), intent(inout)              :: node
    integer(shortInt)                            :: childIdx, depth, i, j, k, nChildren, nIntersectedFaces
    integer(shortInt), dimension(:), allocatable :: elementIdxs, intersectedFaceIdxs, potentialFaceIdxs
    real(defReal), dimension(3)                  :: boundingBoxCentre
    real(defReal), dimension(6)                  :: boundingBoxBounds, childBoundingBoxBounds
    type(axisAlignedBoundingBox), pointer        :: boundingBoxPtr
    type(octreeNode), dimension(:), allocatable  :: temp

    depth = node % getDepth()
    boundingBoxPtr => node % getBoundingBoxPtr()
    boundingBoxBounds = boundingBoxPtr % getBounds()

    ! Use simplified logic for the root Cartesian grid cell.
    if(depth == 1) then
      nIntersectedFaces = faces % getSize()
      intersectedFaceIdxs = [(i, i = 1, nIntersectedFaces)]

    else
      ! Check the number of intersections between the current cell and the faces in the mesh by traversing the
      ! k-d tree starting from the root node.
      call tree % findPotentiallyIntersectedObjects(boundingBoxPtr, potentialFaceIdxs)
      nIntersectedFaces = 0
      do i = 1, size(potentialFaceIdxs)
        ! Use the separating axis theorem to determine if the bounding box intersects the face.
        if (.not. faces % intersectsFace(potentialFaceIdxs(i), boundingBoxPtr)) cycle
        
        ! For now, just increment nIntersectedPrimitives and append the index of the face to the list.
        nIntersectedFaces = nIntersectedFaces + 1
        call append(intersectedFaceIdxs, potentialFaceIdxs(i))

      end do

    end if

    boundingBoxCentre = boundingBoxPtr % getCentre()
    if(nIntersectedFaces <= self % nMaxFaces .or. depth == self % depth) then
      self % nLeaves = self % nLeaves + 1
      if(0 < nIntersectedFaces) then
        call node % setIntersectedFaceIdxs(intersectedFaceIdxs)
        call node % setElementIdxs(faces % getFaceElementIdxs(intersectedFaceIdxs))

      else
        ! Check for element containment.
        elementIdxs = faces % getFaceElementIdxs(tree % findNearestObject(boundingBoxCentre, vertices, faces))
        do i = 1, size(elementIdxs)
          ! Check for inclusion in the current element.
          if(elements % isPointInsideElementNoBoundaryCheck(elementIdxs(i), boundingBoxCentre, faces)) then
            call node % setElementIdxs([elementIdxs(i)])
            return

          end if

        end do

      end if
      return

    end if

    ! If reached here, split the current node into eight children and recurse.
    childIdx = self % nChildren
    call node % setFirstChildIdx(childIdx + 1)
    self % nChildren = self % nChildren + 8
    
    if(.not. allocated(self % children)) then
      allocate(self % children(100))

    else
      nChildren = size(self % children)
      if(nChildren < self % nChildren) then
        allocate(temp(2 * nChildren))
        temp(1:nChildren) = self % children
        call move_alloc(temp, self % children)

      end if

    end if
    
    do i = 1, 2
      do j = 1, 2
        do k = 1, 2
          childIdx = childIdx + 1
          childBoundingBoxBounds(1) = merge(boundingBoxBounds(1), boundingBoxCentre(1), i == 1)
          childBoundingBoxBounds(2) = merge(boundingBoxBounds(2), boundingBoxCentre(2), j == 1)
          childBoundingBoxBounds(3) = merge(boundingBoxBounds(3), boundingBoxCentre(3), k == 1)
          childBoundingBoxBounds(4) = merge(boundingBoxCentre(1), boundingBoxBounds(4), i == 1)
          childBoundingBoxBounds(5) = merge(boundingBoxCentre(2), boundingBoxBounds(5), j == 1)
          childBoundingBoxBounds(6) = merge(boundingBoxCentre(3), boundingBoxBounds(6), k == 1)
          call self % children(childIdx) % init(depth + 1, parentIdx, childBoundingBoxBounds)
          call self % refineNode(childIdx, elements, faces, tree, vertices, self % children(childIdx))

        end do

      end do

    end do

  end subroutine refineNode

end module octree_class