module faceFactory_func

  use numPrecision
  use edge_class,    only : edgeBox
  use face_class,    only : buildFaceInfo, face, faceBox
  use vertex_class,  only : vertexBox

  implicit none
  private

  ! Public interface.
  public :: newFaceBox

contains
  !!
  !!
  !!
  subroutine newFaceBox(info, box)
    type(buildFaceInfo), intent(in) :: info
    type(faceBox), intent(out)      :: box

    ! Initialise face.
    allocate(face :: box % ptr)
    call box % ptr % init(info)

  end subroutine newFaceBox

end module faceFactory_func