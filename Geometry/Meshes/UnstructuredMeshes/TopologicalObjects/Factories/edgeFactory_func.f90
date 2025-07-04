module edgeFactory_func

  use edge_class,   only : buildEdgeInfo, edge, edgeBox
  use numPrecision

  implicit none
  private

  ! Public interface.
  public :: newEdgeBox

contains
  !!
  !!
  !!
  subroutine newEdgeBox(info, box)
    type(buildEdgeInfo), intent(in) :: info
    type(edgeBox), intent(out)      :: box

    allocate(edge :: box % ptr)
    call box % ptr % init(info)

  end subroutine newEdgeBox

end module edgeFactory_func