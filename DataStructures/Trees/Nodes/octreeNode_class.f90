module octreeNode_class

  use element_inter,      only : inclusionTestResult
  use elementShelf_class, only : elementShelf
  use faceShelf_class,    only : faceShelf
  use node_inter,         only : node, kill_super => kill
  use numPrecision
  use universalVariables, only : INSIDE_ELEMENT, ON_BOUNDARY_ELEMENT, OUTSIDE_ELEMENT

  implicit none
  private

  type, public, extends(node) :: octreeNode
    private
    integer(shortInt)                            :: firstChildIdx = 0, level = 0
    integer(shortInt), dimension(:), allocatable :: elementIdxs, intersectedFaceIdxs
  contains
    procedure :: findHostElementIdx
    procedure :: getDepth
    procedure :: getElementIdxs
    procedure :: getFirstChildIdx
    procedure :: getIntersectedFaceIdxs
    procedure :: getStorageSize
    procedure :: init
    procedure :: isInside
    procedure :: isOutside
    procedure :: kill
    procedure :: setElementIdxs
    procedure :: setFirstChildIdx
    procedure :: setIntersectedFaceIdxs
  end type octreeNode

contains
  !!
  !!
  !!
  pure subroutine findHostElementIdx(self, u, elements, faces, elementIdx, r, exitLoop)
    class(octreeNode), intent(in)              :: self
    real(defReal), dimension(3), intent(in)    :: u
    type(elementShelf), intent(in)             :: elements
    type(faceShelf), intent(in)                :: faces
    integer(shortInt), intent(inout)           :: elementIdx
    real(defReal), dimension(3), intent(inout) :: r
    logical(defBool), intent(out)              :: exitLoop
    integer(shortInt)                          :: i, nElements
    type(inclusionTestResult)                  :: testResult

    exitLoop = .true.
    if(self % isOutside()) return

    nElements = size(self % elementIdxs)
    if(nElements == 1) then
      elementIdx = self % elementIdxs(1)
      return

    elseif(size(self % intersectedFaceIdxs) == 1) then
      call faces % testFaceHalfSpace(self % intersectedFaceIdxs(1), r, elementIdx)
      return

    end if

    do i = 1, nElements
      ! Perform inclusion test for the current element.
      testResult = elements % isPointInside(self % elementIdxs(i), r, faces)

      if(testResult % status == INSIDE_ELEMENT) then
        ! If coordinates are fully inside, we have found our element.
        elementIdx = self % elementIdxs(i)
        return

      elseif(testResult % status == ON_BOUNDARY_ELEMENT) then
        ! If coordinates are on the element boundary (very rare), we need to push them off.
        do while (testResult % status == ON_BOUNDARY_ELEMENT)
          call elements % pushFromElementBoundary(self % elementIdxs(i), u, faces, r)

          ! Perform containment test again.
          testResult = elements % isPointInside(self % elementIdxs(i), r, faces)

        end do

        ! Now the coordinates are not on the boundary of the element anymore.
        if(testResult % status == INSIDE_ELEMENT) then
          ! If coordinates are now well inside the element, we have found our element.
          elementIdx = self % elementIdxs(i)
          return

        elseif(testResult % status == OUTSIDE_ELEMENT) then
          ! If the nudge has resulted in an overshoot, we cycle searchLoop and begin the entire process again.
          exitLoop = .false.
          return

        end if

      end if

    end do

  end subroutine findHostElementIdx

  !!
  !!
  !!
  elemental function getDepth(self) result(depth)
    class(octreeNode), intent(in) :: self
    integer(shortInt)             :: depth

    depth = self % level

  end function getDepth

  !!
  !!
  !!
  pure function getElementIdxs(self) result(elementIdxs)
    class(octreeNode), intent(in)                :: self
    integer(shortInt), dimension(:), allocatable :: elementIdxs

    if(allocated(self % elementIdxs)) then
      elementIdxs = self % elementIdxs

    else
      allocate(elementIdxs(0))

    end if

  end function getElementIdxs

  !!
  !!
  !!
  elemental function getFirstChildIdx(self) result(firstChildIdx)
    class(octreeNode), intent(in) :: self
    integer(shortInt)             :: firstChildIdx

    firstChildIdx = self % firstChildIdx

  end function getFirstChildIdx

  !!
  !!
  !!
  pure function getIntersectedFaceIdxs(self) result(intersectedFaceIdxs)
    class(octreeNode), intent(in)                :: self
    integer(shortInt), dimension(:), allocatable :: intersectedFaceIdxs

    if(allocated(self % intersectedFaceIdxs)) then
      intersectedFaceIdxs = self % intersectedFaceIdxs

    else
      allocate(intersectedFaceIdxs(0))

    end if

  end function getIntersectedFaceIdxs

  !!
  !!
  !!
  elemental function getStorageSize(self) result(storageSize)
    class(octreeNode), intent(in) :: self
    integer(longInt)              :: storageSize

    storageSize = storage_size(self) / 8
    if(allocated(self % elementIdxs)) storageSize = storageSize + 4 * size(self % elementIdxs)
    if(allocated(self % intersectedFaceIdxs)) storageSize = storageSize + 4 * size(self % intersectedFaceIdxs)

  end function getStorageSize

  !! Subroutine 'init'
  !!
  !! Basic description:
  !!   
  subroutine init(self, depth, parentIdx, boundingBoxBounds)
    class(octreeNode), intent(inout)               :: self
    integer(shortInt), intent(in)                  :: depth, parentIdx
    real(defReal), dimension(6), intent(in)        :: boundingBoxBounds

    ! Set the node's depth, parent index, and bounding box.
    self % level = depth
    call self % setParentIdx(parentIdx)
    call self % initBoundingBox(boundingBoxBounds)

  end subroutine init

  !!
  !!
  !!
  elemental function isInside(self) result(isIt)
    class(octreeNode), intent(in) :: self
    logical(defBool)              :: isIt

    if(allocated(self % elementIdxs)) then
      isIt = size(self % elementIdxs) == 1

    else
      isIt = .false.

    end if

  end function isInside

  !!
  !!
  !!
  elemental function isOutside(self) result(isIt)
    class(octreeNode), intent(in) :: self
    logical(defBool)              :: isIt

    isIt = .not. allocated(self % elementIdxs) .and. .not. allocated(self % intersectedFaceIdxs)

  end function isOutside

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
    self % firstChildIdx = 0
    self % level = 0
    if(allocated(self % elementIdxs)) deallocate(self % elementIdxs)
    if(allocated(self % intersectedFaceIdxs)) deallocate(self % intersectedFaceIdxs)

  end subroutine kill

  !!
  !!
  !!
  pure subroutine setElementIdxs(self, elementIdxs)
    class(octreeNode), intent(inout)            :: self
    integer(shortInt), dimension(:), intent(in) :: elementIdxs

    self % elementIdxs = elementIdxs

  end subroutine setElementIdxs

  !!
  !!
  !!
  elemental subroutine setFirstChildIdx(self, firstChildIdx)
    class(octreeNode), intent(inout) :: self
    integer(shortInt), intent(in)    :: firstChildIdx

    self % firstChildIdx = firstChildIdx

  end subroutine setFirstChildIdx

  !!
  !!
  !!
  pure subroutine setIntersectedFaceIdxs(self, intersectedFaceIdxs)
    class(octreeNode), intent(inout)            :: self
    integer(shortInt), dimension(:), intent(in) :: intersectedFaceIdxs

    self % intersectedFaceIdxs = intersectedFaceIdxs

  end subroutine setIntersectedFaceIdxs

end module octreeNode_class