module vertexFactory_func

  use numPrecision
  use vertex_class, only : vertex

  implicit none
  private

  ! Public interface.
  public :: newVertexPtr

contains
  !!
  !!
  !!
  subroutine newVertexPtr(idx, coords, vertexPtr)
    integer(shortInt), intent(in)           :: idx
    real(defReal), dimension(3), intent(in) :: coords
    type(vertex), pointer, intent(out)      :: vertexPtr

    ! Initialise vertex.
    allocate(vertex :: vertexPtr)
    call vertexPtr % init(idx, coords)

  end subroutine newVertexPtr

end module vertexFactory_func