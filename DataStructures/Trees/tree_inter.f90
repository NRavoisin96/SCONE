module tree_inter

  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use coord_class,                  only : coord
  use dictionary_class,             only : dictionary
  use genericProcedures,            only : fatalError
  use node_inter,                   only : buildNodePayload, node
  use numPrecision
  use topologicalObject_inter,      only : topologicalObjectBox
  use universalVariables,           only : INF

  implicit none
  private

  ! Extendable procedures.
  public :: kill

  !!
  !!
  !!
  type, public, abstract :: tree
    private
    integer(shortInt)    :: maxDepth = 0, nLeaves = 0, nNodes = 0, nObjects = 0
    class(node), pointer :: root
  contains
    ! Build procedures.
    procedure(allocateRoot), deferred :: allocateRoot
    procedure                         :: assignElements
    procedure                         :: init
    procedure                         :: kill
    ! Runtime procedures.
    procedure                         :: findLeaf
    generic                           :: findIntersectedObjects => findIntersectedObjects_BoundingBox
    procedure, private                :: findIntersectedObjects_BoundingBox
    procedure                         :: findNearestObject
    procedure                         :: getLeavesNumber
    procedure                         :: getNodesNumber
    procedure                         :: getObjectsNumber
    procedure                         :: getRootBoundingBoxBounds
  end type tree

  abstract interface
    !!
    !!
    !!
    subroutine allocateRoot(self, ptr)
      import                            :: node, tree
      class(tree), intent(in)           :: self
      class(node), pointer, intent(out) :: ptr
    end subroutine allocateRoot

  end interface

contains
  !!
  !!
  !!
  subroutine assignElements(self, payload)
    class(tree), intent(in)             :: self
    class(buildNodePayload), intent(in) :: payload

    call self % root % assignElements(payload)

  end subroutine assignElements

  !!
  !!
  !!
  subroutine findLeaf(self, coords, leaf, requiresContainmentCheck)
    class(tree), intent(in)                :: self
    type(coord), intent(inout)             :: coords
    class(node), pointer, intent(out)      :: leaf
    logical(defBool), intent(in), optional :: requiresContainmentCheck

    call self % root % findLeaf(coords, leaf, requiresContainmentCheck)

  end subroutine findLeaf

  !!
  !!
  !!
  function findIntersectedObjects_BoundingBox(self, boundingBox) result(intersectedObjects)
    class(tree), intent(in)                               :: self
    type(axisAlignedBoundingBox), intent(in)              :: boundingBox
    type(topologicalObjectBox), dimension(:), allocatable :: intersectedObjects

    ! Initialise idxs to a zero-sized array then query the root node.
    allocate(intersectedObjects(0))
    call self % root % findIntersectedObjects(boundingBox, intersectedObjects)

  end function findIntersectedObjects_BoundingBox

  !!
  !!
  !!
  function findNearestObject(self, r) result(nearestObject)
    class(tree), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    type(topologicalObjectBox)              :: nearestObject
    real(defReal)                           :: radiusSquared

    ! Initialise radiusSquared = INF then search root node.
    radiusSquared = INF
    call self % root % findNearestObject(r, radiusSquared, nearestObject)

  end function findNearestObject

  !!
  !!
  !!
  elemental function getLeavesNumber(self) result(nLeaves)
    class(tree), intent(in) :: self
    integer(shortInt)       :: nLeaves

    nLeaves = self % nLeaves

  end function getLeavesNumber

  !!
  !!
  !!
  elemental function getNodesNumber(self) result(nNodes)
    class(tree), intent(in) :: self
    integer(shortInt)       :: nNodes

    nNodes = self % nNodes

  end function getNodesNumber

  !!
  !!
  !!
  elemental function getObjectsNumber(self) result(nObjects)
    class(tree), intent(in) :: self
    integer(shortInt)       :: nObjects

    nObjects = self % nObjects

  end function getObjectsNumber

  !!
  !!
  !!
  pure function getRootBoundingBoxBounds(self) result(bounds)
    class(tree), intent(in)        :: self
    real(defReal), dimension(3, 2) :: bounds

    bounds = self % root % getBoundingBoxBounds()

  end function getRootBoundingBoxBounds

  !!
  !!
  !!
  subroutine init(self, payload, dict)
    class(tree), intent(inout)              :: self
    class(buildNodePayload), intent(inout)  :: payload
    class(dictionary), intent(in), optional :: dict
    character(*), parameter                 :: here = 'init (tree_inter.f90)'

    ! Check if dict was supplied.
    if (present(dict)) then
      ! Get maximum depth from dictionary. Use 0 if no maximum depth was specified.
      call dict % getOrDefault(payload % maxDepth, 'maxDepth', 0)

      ! Get bucket size from dictionary. Use 4 if no bucket size was specified.
      call dict % getOrDefault(payload % bucketSize, 'bucketSize', 4)

    end if

    ! Check payload inputs.
    if (payload % bucketSize == 0) call fatalError(here, 'Bucket size must be positive.')

    ! First allocate the root node.
    call self % allocateRoot(self % root)

    ! Now build the root node.
    call self % root % init(payload)

    ! Unpack finalised information from payload.
    self % maxDepth = payload % maxDepth
    self % nLeaves = payload % nLeaves
    self % nNodes = payload % nNodes
    self % nObjects = payload % shelf % getObjectsNumber()

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(tree), intent(inout) :: self

    ! Local.
    call self % root % kill()
    deallocate(self % root)
    nullify(self % root)
    self % maxDepth = 0
    self % nLeaves = 0
    self % nNodes = 0
    self % nObjects = 0

  end subroutine

end module tree_inter