module octreeNode_class

  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use coord_class,                  only : coord
  use element_inter,                only : inclusionTestResult
  use elementShelf_class,           only : elementShelf
  use face_inter,                   only : face
  use faceShelf_class,              only : faceShelf
  use genericProcedures,            only : append, areEqual, fatalError
  use objectKDTree_class,           only : objectKDTree
  use numPrecision
  use node_inter,                   only : node, kill_super => kill
  use universalVariables,           only : HALF, INF, INSIDE_ELEMENT, NUDGE, ON_BOUNDARY_ELEMENT, OUTSIDE_ELEMENT
  use vertexShelf_class,            only : vertexShelf

  implicit none
  private

  type, public, extends(node) :: octreeNode
    private
    integer(shortInt)                            :: level = 0, nIntersectingFaces = 0
    integer(shortInt), dimension(:), allocatable :: elementIdxs
    logical(defBool)                             :: isUnchecked = .true., isInside = .false., isOutside = .false., &
                                                    isIntersecting = .false.
    type(octreeNode), dimension(:), allocatable  :: children
    type(octreeNode), pointer                    :: parent => null()
  contains
    ! Build procedures.
    procedure :: assignElement
    procedure :: init
    procedure :: kill
    procedure :: refine
    procedure :: split
    ! Runtime procedures.
    procedure :: countInside
    procedure :: countOutside
    procedure :: findLeaf
    procedure :: getElementIdxs
    procedure :: getIsInside
    procedure :: getIsIntersecting
    procedure :: getIsOutside
  end type octreeNode

contains
  !!
  !!
  !!
  recursive subroutine assignElement(self, tree, vertices, faces, elements)
    class(octreeNode), intent(inout)             :: self
    type(objectKDTree), intent(in)               :: tree
    type(vertexShelf), intent(in)                :: vertices
    type(faceShelf), intent(in)                  :: faces
    type(elementShelf), intent(in)               :: elements
    integer(shortInt)                            :: elementIdx, i, nearestFaceIdx
    integer(shortInt), dimension(:), allocatable :: elementIdxs
    real(defReal), dimension(3)                  :: boundingBoxCentre
    type(inclusionTestResult)                    :: insideResult

    ! If the cell is not a leaf, descend deeper into the tree.
    if (.not. self % getIsLeaf()) then
      do i = 1, 8
        call self % children(i) % assignElement(tree, vertices, faces, elements)

      end do
      return

    end if

    ! If the cell has already been checked simply return.
    if (.not. self % isUnchecked) return

    ! If the cell is unchecked, find the nearest mesh face from the tree.
    self % isUnchecked = .false.
    boundingBoxCentre = self % getBoundingBoxCentre()
    nearestFaceIdx = tree % findNearestObject(boundingBoxCentre, vertices, faces)
    
    ! Retrieve the elements associated with the nearest face.
    elementIdxs = faces % getFaceElementIdxs(nearestFaceIdx)
    do i = 1, size(elementIdxs)
      ! Check for inclusion in the current element.
      elementIdx = elementIdxs(i)
      insideResult = elements % isPointInside(elementIdx, boundingBoxCentre, faces)

      ! If the current element contains the bounding box, assign it to the cell
      ! and return.
      if (insideResult % status == INSIDE_ELEMENT) then
        allocate(self % elementIdxs(1))
        self % elementIdxs(1) = elementIdx
        self % isInside = .true.
        return

      end if

    end do

    ! If reached here, the cell is not in any element, so set it as being outside.
    self % isOutside = .true.

  end subroutine assignElement

  !!
  !!
  !!
  pure recursive subroutine countInside(self, nInside)
    class(octreeNode), intent(in)    :: self
    integer(shortInt), intent(inout) :: nInside
    integer(shortInt)                :: i

    if (.not. self % getIsLeaf()) then
      do i = 1, 8
        call self % children(i) % countInside(nInside)

      end do
      return

    end if

    if (self % isInside) nInside = nInside + 1

  end subroutine countInside

  !!
  !!
  !!
  pure recursive subroutine countOutside(self, nOutside)
    class(octreeNode), intent(in)    :: self
    integer(shortInt), intent(inout) :: nOutside
    integer(shortInt)                :: i

    if (.not. self % getIsLeaf()) then
      do i = 1, 8
        call self % children(i) % countOutside(nOutside)

      end do
      return

    end if

    if (self % isOutside) nOutside = nOutside + 1

  end subroutine countOutside

  !!
  !!
  !!
  recursive subroutine findLeaf(self, coords, leaf, requiresContainmentCheck)
    class(octreeNode), intent(in), target  :: self
    type(coord), intent(inout)             :: coords
    type(octreeNode), intent(out), pointer :: leaf
    logical(defBool), intent(in), optional :: requiresContainmentCheck
    logical(defBool)                       :: checkContainment, inside
    real(defReal), dimension(3)            :: boundingBoxCentre, r
    integer(shortInt)                      :: i, idx

    ! Only perform containment check if it has been required.
    checkContainment = .false.
    if (present(requiresContainmentCheck)) checkContainment = requiresContainmentCheck

    ! Perform containment check if needed.
    if (checkContainment) then
      if (.not. self % boundingBoxContains(coords % getPositionToNudge())) then
        if (associated(self % parent)) then
          call self % parent % findLeaf(coords, leaf, .true.)

        else
          leaf => null()

        end if
        return

      end if

    end if

    ! Push coordinates from boundary of bounding box if applicable.
    call self % pushFromBoundingBoxBoundary(coords, inside)

    ! Check for overshoot.
    if (.not. inside) then
      if (associated(self % parent)) then
        call self % parent % findLeaf(coords, leaf, .true.)

        else
            ! We are at the root and overshot. Particle is outside the domain.
            leaf => null()

        end if
        return

    end if

    ! If cell is a leaf, simply associate the leaf pointer and return.
    if (self % getIsLeaf()) then
      leaf => self
      return

    end if

    ! Retrieve coordinates position and descend into correct child node.
    r = coords % getPositionToNudge()
    boundingBoxCentre = self % getBoundingBoxCentre()
    idx = 1
    if (r(1) >= boundingBoxCentre(1)) idx = idx + 4
    if (r(2) >= boundingBoxCentre(2)) idx = idx + 2
    if (r(3) >= boundingBoxCentre(3)) idx = idx + 1

    ! Check the correct child cell.
    call self % children(idx) % findLeaf(coords, leaf)

  end subroutine findLeaf

  !!
  !!
  !!
  pure function getElementIdxs(self) result(elementIdxs)
    class(octreeNode), intent(in)                :: self
    integer(shortInt), dimension(:), allocatable :: elementIdxs

    ! Check if the elementIdxs component is allocated and return empty array if not.
    if (.not. allocated(self % elementIdxs)) then
      allocate(elementIdxs(0))

    else
      elementIdxs = self % elementIdxs

    end if

  end function getElementIdxs

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

  !! Subroutine 'init'
  !!
  !! Basic description:
  !!   
  subroutine init(self, boundingBoxBounds, tree, vertices, faces, level, maxFacesNumber, maxRefinementLevel, &
                  nLeaves, parent)
    class(octreeNode), intent(inout)               :: self
    real(defReal), dimension(6), intent(in)        :: boundingBoxBounds
    type(objectKDTree), intent(in)                 :: tree
    type(vertexShelf), intent(in)                  :: vertices
    type(faceShelf), intent(in)                    :: faces
    integer(shortInt), intent(in)                  :: level, maxFacesNumber, maxRefinementLevel
    integer(shortInt), intent(inout)               :: nLeaves
    type(octreeNode), intent(in), target, optional :: parent

    ! Set the cell's bounding box and level, then begin the recursive refinement procedure.
    call self % initBoundingBox(boundingBoxBounds)
    self % level = level
    if (present(parent)) self % parent => parent
    call self % refine(tree, vertices, faces, maxFacesNumber, maxRefinementLevel, nLeaves)

  end subroutine init

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  pure recursive subroutine kill(self)
    class(octreeNode), intent(inout) :: self
    integer(shortInt)                :: i

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % level = 0
    self % nIntersectingFaces = 0
    if (allocated(self % elementIdxs)) deallocate(self % elementIdxs)
    self % isUnchecked = .true.
    self % isInside = .false.
    self % isOutside = .false.
    self % isIntersecting = .false.
    if (associated(self % parent)) nullify(self % parent)
    if (allocated(self % children)) then
      do i = 1, size(self % children)
        call self % children(i) % kill()

      end do
      deallocate(self % children)

    end if

  end subroutine kill

  !! Subroutine 'refine'
  !!
  !! Basic description:
  !!   Recursively refines a Cartesian grid cell based on the number of mesh faces it intersects.
  !!
  recursive subroutine refine(self, tree, vertices, faces, maxFacesNumber, maxRefinementLevel, nLeaves)
    class(octreeNode), intent(inout)      :: self
    type(objectKDTree), intent(in)               :: tree
    type(vertexShelf), intent(in)                :: vertices
    type(faceShelf), intent(in)                  :: faces
    integer(shortInt), intent(in)                :: maxFacesNumber, maxRefinementLevel
    integer(shortInt), intent(inout)             :: nLeaves
    integer(shortInt)                            :: potentialFaceIdx, i, nFaces, nIntersectedFaces
    integer(shortInt), dimension(:), allocatable :: potentialFaceIdxs, intersectedFaceIdxs
    type(axisAlignedBoundingBox)                 :: boundingBox

    ! Use simplified logic for the root Cartesian grid cell.
    if (self % level == 1) then
      self % isUnchecked = .false.
      self % isIntersecting = .true.
      ! Retrieve the number of faces in the tree.
      nFaces = tree % getDataNumber()

      ! If nFaces <= maxFacesNumber (very simple unstructured mesh geometries), there is no need
      ! to refine the cell and we can simply return.
      if (nFaces <= maxFacesNumber) then
        call self % setIsLeaf()
        nLeaves = nLeaves + 1
        self % nIntersectingFaces = faces % getSize()
        allocate(intersectedFaceIdxs(self % nIntersectingFaces))
        do i = 1, self % nIntersectingFaces
          intersectedFaceIdxs(i) = i

        end do
        self % elementIdxs = faces % getFaceElementIdxs(intersectedFaceIdxs)
        return

      end if

      ! Else, split the cell and return.
      call self % split(tree, vertices, faces, maxFacesNumber, maxRefinementLevel, nLeaves)
      return

    end if

    ! Check the number of intersections between the current cell and the faces in the mesh by traversing the
    ! k-d tree starting from the root node.
    boundingBox = self % getBoundingBox()
    call tree % findPotentiallyIntersectedObjects(boundingBox, potentialFaceIdxs)
    nIntersectedFaces = 0
    do i = 1, size(potentialFaceIdxs)
      ! First check if bounding box intersects the current face's bounding box.
      potentialFaceIdx = potentialFaceIdxs(i)
      if (.not. faces % intersectsFaceBoundingBox(potentialFaceIdx, boundingbox)) cycle

      ! If the bounding boxes intersect, perform a test based on the separating axis theorem to determine if the bounding box actually
      ! intersects the face.
      if (.not. faces % intersectsFace(potentialFaceIdx, boundingbox)) cycle
      
      ! For now, just increment nIntersectedPrimitives and append the index of the face to the list.
      nIntersectedFaces = nIntersectedFaces + 1
      call append(intersectedFaceIdxs, potentialFaceIdx)

    end do

    self % nIntersectingFaces = nIntersectedFaces
    if (nIntersectedFaces == 0) then
      call self % setIsLeaf()
      nLeaves = nLeaves + 1

    else
      self % isUnchecked = .false.
      self % isIntersecting = .true.
      if (self % level == maxRefinementLevel .or. nIntersectedFaces <= maxFacesNumber) then
        call self % setIsLeaf()
        nLeaves = nLeaves + 1
        self % elementIdxs = faces % getFaceElementIdxs(intersectedFaceIdxs)

      else
        call self % split(tree, vertices, faces, maxFacesNumber, maxRefinementLevel, nLeaves)

      end if

    end if

  end subroutine refine

  !! Subroutine 'split'
  !!
  !! Basic description:
  !!   Splits the current Cartesian grid cell into eight children. Computes the bounding box of each child
  !!   from the current cell's bounding box, then initialises each child cell.
  !!
  subroutine split(self, tree, vertices, faces, maxFacesNumber, maxRefinementLevel, nLeaves)
    class(octreeNode), intent(inout) :: self
    type(objectKDTree), intent(in)          :: tree
    type(vertexShelf), intent(in)           :: vertices
    type(faceShelf), intent(in)             :: faces
    integer(shortInt), intent(in)           :: maxFacesNumber, maxRefinementLevel
    integer(shortInt), intent(inout)        :: nLeaves
    integer(shortInt)                       :: i, j, k, childIdx
    real(defReal), dimension(3)             :: boundingBoxCentre
    real(defReal), dimension(6)             :: boundingBoxBounds, childBoundingBoxBounds

    ! Allocate 8 children cells for the current cell.
    allocate(self % children(8))

    ! Compute the centre coordinates of the current cell's bounding box.
    boundingBoxBounds = self % getBoundingBoxBounds()
    boundingBoxCentre = self % getBoundingBoxCentre()

    ! Loop through all children cells and initialise them.
    childIdx = 0
    do i = 1, 2
      do j = 1, 2
        do k = 1, 2
          ! Increment childIdx and compute the bounding box of the new cell.
          childIdx = childIdx + 1
          if (i == 1) then
            childBoundingBoxBounds(1) = boundingBoxBounds(1)
            childBoundingBoxBounds(4) = boundingBoxCentre(1)

          else
            childBoundingBoxBounds(1) = boundingBoxCentre(1)
            childBoundingBoxBounds(4) = boundingBoxBounds(4)

          end if

          if (j == 1) then
            childBoundingBoxBounds(2) = boundingBoxBounds(2)
            childBoundingBoxBounds(5) = boundingBoxCentre(2)

          else
            childBoundingBoxBounds(2) = boundingBoxCentre(2)
            childBoundingBoxBounds(5) = boundingBoxBounds(5)

          end if

          if (k == 1) then
            childBoundingBoxBounds(3) = boundingBoxBounds(3)
            childBoundingBoxBounds(6) = boundingBoxCentre(3)

          else
            childBoundingBoxBounds(3) = boundingBoxCentre(3)
            childBoundingBoxBounds(6) = boundingBoxBounds(6)

          end if

          ! Initialise the new cell with the computed bounding box.
          call self % children(childIdx) % init(childBoundingBoxBounds, tree, vertices, faces, &
                                                self % level + 1, maxFacesNumber, maxRefinementLevel, nLeaves, self)

        end do

      end do

    end do

  end subroutine split

end module octreeNode_class