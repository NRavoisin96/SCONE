module cartesianGridCell_class

  implicit none
  private

  !!
  !!
  !!
  type, public :: cartesianGridCell
    private
  contains
    procedure :: kill
  end type cartesianGridCell

contains
  !!
  !!
  !!
  elemental subroutine kill(self)
    class(cartesianGridCell), intent(inout) :: self

    ! Local.

  end subroutine kill

end module cartesianGridCell_class