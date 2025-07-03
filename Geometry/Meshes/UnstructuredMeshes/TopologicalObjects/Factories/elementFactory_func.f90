module elementFactory_func

  use element_class, only : buildElementInfo, element, elementBox
  use numPrecision

  implicit none
  private

  ! Public interface.
  public :: newElementBox

contains
  !!
  !!
  !!
  subroutine newElementBox(info, box)
    type(buildElementInfo), intent(in) :: info
    type(elementBox), intent(out)      :: box

    allocate(element :: box % ptr)
    call box % ptr % init(info)

  end subroutine newElementBox

end module elementFactory_func