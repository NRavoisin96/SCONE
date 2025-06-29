module edgeFactory_func

  use edge_class,   only : edge, edgeBox
  use numPrecision
  use vertex_class, only : vertexBox

  implicit none
  private

  ! Public interface.
  public :: newEdgeBox

contains
  !!
  !!
  !!
  subroutine newEdgeBox(idx, vertices, box)
    integer(shortInt), intent(in)             :: idx
    type(vertexBox), dimension(2), intent(in) :: vertices
    type(edgeBox), intent(out)                :: box

    allocate(edge :: box % ptr)
    call box % ptr % init(idx, vertices)

  end subroutine newEdgeBox

end module edgeFactory_func