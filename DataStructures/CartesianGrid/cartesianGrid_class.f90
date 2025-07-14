module cartesianGrid_class

  use cartesianGridCell_class, only : cartesianGridCell
  use numPrecision

  implicit none
  private

  !!
  !!
  !!
  type, public :: cartesianGrid
    private
    type(cartesianGridCell), dimension(:), allocatable :: cells
  contains
    procedure :: init
    procedure :: kill
  end type cartesianGrid

contains
  !!
  !!
  !!
  subroutine init(self)
    class(cartesianGrid), intent(inout) :: self

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(cartesianGrid), intent(inout) :: self
    integer(shortInt)                   :: i

    ! Local.
    if (allocated(self % cells)) then
      do i = 1, size(self % cells)
        call self % cells(i) % kill()

      end do
      deallocate(self % cells)

    end if

  end subroutine kill

end module cartesianGrid_class