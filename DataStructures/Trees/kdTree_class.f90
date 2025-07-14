module kdTree_class

  use kdTreeNode_class, only : kdTreeNode
  use node_inter,       only : node
  use tree_inter,       only : tree

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(tree) :: kdTree
    private
  contains
    procedure :: allocateRoot
  end type kdTree

contains
  !!
  !!
  !!
  subroutine allocateRoot(self, ptr)
    class(kdTree), intent(in)         :: self
    class(node), pointer, intent(out) :: ptr

    allocate(kdTreeNode :: ptr)

  end subroutine allocateRoot

end module kdTree_class