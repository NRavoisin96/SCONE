module tree_inter

  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use node_inter,                   only : node
  use objectNode_class,             only : objectNode
  use numPrecision

  implicit none
  private

  ! Extendable procedures.
  public :: kill

  !!
  !!
  !!
  type, abstract             :: tree
    integer(shortInt)        :: nLeaves = 0, nNodes = 0
    class(node), allocatable :: root
  contains
    ! Build procedures.
    procedure :: kill
    ! Runtime procedures.
    procedure :: getLeavesNumber
    procedure :: getNodesNumber
  end type tree

contains
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
  elemental subroutine kill(self)
    class(tree), intent(inout) :: self

    ! Local.
    call self % root % kill()
    self % nLeaves = 0
    self % nNodes = 0

  end subroutine

end module tree_inter