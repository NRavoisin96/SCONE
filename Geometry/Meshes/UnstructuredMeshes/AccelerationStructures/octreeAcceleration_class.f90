module octreeAcceleration_class

  use accelerationStructure_inter, only : accelerationStructure, getStorageSize_super => getStorageSize
  use dictionary_class,            only : dictionary
  use edgeShelf_class,             only : edgeShelf
  use element_inter,               only : inclusionTestResult
  use elementShelf_class,          only : elementShelf
  use errors_mod,                  only : fatalError
  use faceShelf_class,             only : faceShelf
  use numPrecision
  use octree_class,                only : octree
  use octreeNode_class,            only : octreeNode
  use vertexShelf_class,           only : vertexShelf
  use universalVariables,          only : INSIDE_ELEMENT, ON_BOUNDARY_ELEMENT, OUTSIDE_ELEMENT

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(accelerationStructure) :: octreeAcceleration
    private
    type(octree) :: tree
  contains
    procedure :: findHostElementIdx
    procedure :: getStorageSize
    procedure :: init
    procedure :: kill
  end type octreeAcceleration

contains
  !!
  !!
  !!
  subroutine findHostElementIdx(self, u, edges, elements, faces, vertices, elementIdx, r)
    class(octreeAcceleration), intent(in)        :: self
    real(defReal), dimension(3), intent(in)      :: u
    type(edgeShelf), intent(in)                  :: edges
    type(elementShelf), intent(in)               :: elements
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    integer(shortInt), intent(inout)             :: elementIdx
    real(defReal), dimension(3), intent(inout)   :: r
    logical(defBool)                             :: exitLoop
    type(octreeNode), pointer                    :: leaf

    searchLoop: do
      ! First search the acceleration structure for the indices of potential elements containing the coordinates.
      call self % tree % findLeaf(r, leaf)

      ! If check if the leaf node pointer is associated.
      if(.not. associated(leaf)) return

      ! If leaf is entirely contained within an element or is outside the mesh, we are done.
      call leaf % findHostElementIdx(u, elements, faces, elementIdx, r, exitLoop)
      if(exitLoop) return

    end do searchLoop

  end subroutine findHostElementIdx

  !!
  !!
  !!
  elemental function getStorageSize(self) result(storageSize)
    class(octreeAcceleration), intent(in) :: self
    integer(longInt)                      :: storageSize

    storageSize = getStorageSize_super(self) + self % tree % getStorageSize()

  end function getStorageSize

  !!
  !!
  !!
  subroutine init(self, dict, vertices, edges, faces, elements)
    class(octreeAcceleration), intent(inout) :: self
    type(dictionary), intent(in)             :: dict
    type(vertexShelf), intent(in)            :: vertices
    type(faceShelf), intent(inout)           :: faces
    type(elementShelf), intent(in)           :: elements
    type(edgeShelf), intent(inout)           :: edges

    ! Simply initialise the octree.
    call self % tree % init(dict, edges, elements, vertices, faces)

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(octreeAcceleration), intent(inout) :: self

    ! Local.
    call self % tree % kill()

  end subroutine kill

end module octreeAcceleration_class