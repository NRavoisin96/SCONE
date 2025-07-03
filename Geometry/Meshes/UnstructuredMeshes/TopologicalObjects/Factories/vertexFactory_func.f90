module vertexFactory_func

  use numPrecision
  use vertex_class, only : vertex, vertexBox

  implicit none
  private

  ! Public interface.
  public :: newVertexBox

contains
  !!
  !!
  !!
  subroutine newVertexBox(idx, coords, box)
    integer(shortInt), intent(in)           :: idx
    real(defReal), dimension(3), intent(in) :: coords
    type(vertexBox), intent(out)            :: box

    ! Initialise vertex.
    allocate(vertex :: box % ptr)
    call box % ptr % init(idx, coords)

  end subroutine newVertexBox

end module vertexFactory_func