module octree_class

  use node_inter,         only : node
  use octreeNode_class,   only : octreeNode
  use tree_inter,         only : tree

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