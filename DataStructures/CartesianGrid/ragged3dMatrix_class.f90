module ragged3dMatrix_class

  use numPrecision,             only : defReal

  implicit none
  private

  !!
  !!
  !!
  type :: ragged2d
     real(defReal), allocatable :: data(:,:)
  end type ragged2d

  !!
  !!
  !!
  type, public :: ragged3d
    private
    type(ragged2d), allocatable :: slice(:)
  contains
    procedure :: init       ! construct empty or with capacity
    procedure :: kill       ! free all storage
    procedure :: append     ! append one 2D matrix
    procedure :: scale      ! scale all stored elements
    procedure :: get        ! retrieve a copy of slice k
    procedure :: take       ! move-allocate slice k out
    procedure :: nslices    ! returns the number of slices
  end type ragged3d

contains

  !!
  !!
  !!
  subroutine init(self, n)
    class(ragged3d), intent(inout) :: self
    integer, intent(in), optional  :: n

    if (allocated(self%slice)) call kill(self)
    if (present(n) .and. n > 0) then
      allocate(self%slice(n))
    end if

  end subroutine init

  !!
  !!
  !!
  subroutine kill(self)
    class(ragged3d), intent(inout) :: self
    integer :: i

    if (allocated(self%slice)) then
      do i = 1, size(self%slice)
        if (allocated(self%slice(i)%data)) deallocate(self%slice(i)%data)
      end do
      deallocate(self%slice)
    end if

  end subroutine kill

  !!
  !!
  !!
  subroutine append(self, src)
    class(ragged3d), intent(inout) :: self
    real(defReal),  intent(in)     :: src(:,:)
    type(ragged2d), allocatable    :: tmp(:)
    integer :: n

    n = nslices(self)
    allocate(tmp(n+1))
    if (n > 0) tmp(1:n) = self%slice

    allocate(tmp(n+1)%data, mold=src)
    tmp(n+1)%data = src

    call move_alloc(tmp, self%slice)

  end subroutine append

  !!
  !!
  !!
  subroutine scale(self, alpha)
    class(ragged3d), intent(inout) :: self
    real(defReal),  intent(in)     :: alpha
    integer :: i

    do i = 1, nslices(self)
      if (allocated(self%slice(i)%data)) self%slice(i)%data = alpha * self%slice(i)%data
    end do

  end subroutine scale

  !!
  !!
  !!
  function get(self, k) result(out)
    class(ragged3d), intent(in) :: self
    integer,         intent(in) :: k
    real(defReal), allocatable  :: out(:,:)

    if (.not. allocated(self%slice)) return
    if (k < 1 .or. k > size(self%slice)) return
    if (.not. allocated(self%slice(k)%data)) return

    allocate(out, mold=self%slice(k)%data)
    out = self%slice(k)%data

  end function get

  !!
  !!
  !!
  subroutine take(self, k, out)
    class(ragged3d),       intent(inout) :: self
    integer,               intent(in)    :: k
    real(defReal), allocatable, intent(out) :: out(:,:)

    if (.not. allocated(self%slice)) return   
    if (k < 1 .or. k > size(self%slice)) return  
    if (.not. allocated(self%slice(k)%data)) return

    call move_alloc(self%slice(k)%data, out)

  end subroutine take

  !!
  !!
  !!
  pure function nslices(self) result(n)
    class(ragged3d), intent(in) :: self
    integer :: n
    if (allocated(self%slice)) then
      n = size(self%slice)
    else
      n = 0
    end if
    
  end function nslices

end module ragged3dMatrix_class
