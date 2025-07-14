module octree_class

  use coord_class,        only : coord
  use dictionary_class,   only : dictionary
  use elementShelf_class, only : elementShelf
  use faceShelf_class,    only : faceShelf
  use genericProcedures,  only : areEqual, fatalError
  use kdTree_class,       only : kdTree
  use mesh_inter,         only : mesh
  use node_inter,         only : node
  use numPrecision
  use octreeNode_class,   only : buildOctreeNodePayload, octreeNode
  use tree_inter,         only : tree
  use universalVariables, only : NUDGE
  use vertexShelf_class,  only : vertexShelf

  implicit none
  private

  type, public, extends(tree) :: octree
    private
  contains
    ! Build procedures.
    procedure :: allocateRoot
    ! Runtime procedures.
  end type octree

contains
  !!
  !!
  !!
  subroutine allocateRoot(self, ptr)
    class(octree), intent(in)         :: self
    class(node), pointer, intent(out) :: ptr

    allocate(octreeNode :: ptr)

  end subroutine allocateRoot

end module octree_class