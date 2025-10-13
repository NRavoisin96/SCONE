module dynamic2dMatSet_class

  use numPrecision,        only : defReal, shortInt
  use genericProcedures,   only : fatalError

  implicit none
  private

  !!
  !! Tiny linear multiset for small N: keys(:), cnts(:), n.
  !!
  type :: smallIntCounts
     integer(shortInt), allocatable :: keys(:)
     integer(shortInt), allocatable :: cnts(:)
     integer(shortInt) :: n = 0
  contains
     procedure :: init         => sic_init
     procedure :: kill         => sic_kill
     procedure :: add          => sic_add
     procedure :: sub          => sic_sub
     procedure :: unique_count => sic_unique_count
     procedure :: to_list      => sic_to_list
     procedure, private :: grow => sic_grow
  end type smallIntCounts

  !!
  !! Single 2D matrix (3, ncols) with per-column tags.
  !!
  type :: dynamic2dMat
     real(defReal),       allocatable :: data(:,:)
     integer(shortInt),   allocatable :: tags(:)
  end type dynamic2dMat

  !!
  !! Dynamic set of 2D matrices of shape (3, n_i).
  !!
  type, public :: dynamic2dMatSet
    private
    type(dynamic2dMat), allocatable :: slice(:)
    integer(shortInt) :: nused = 0
    type(smallIntCounts) :: gids
  contains
    procedure :: init
    procedure :: kill
    procedure :: append
    procedure :: nslices
    procedure :: scale
    procedure :: update_replace
    procedure :: filter_columns
    procedure :: delete
    procedure :: delete_columns
    procedure :: get_copy
    procedure :: unique_count
    procedure :: is_singleton
    procedure :: unique_list
  end type dynamic2dMatSet

contains

!========================
! smallIntCounts methods
!========================

  subroutine sic_init(self, cap_hint)
    class(smallIntCounts), intent(inout) :: self
    integer(shortInt), intent(in), optional :: cap_hint
    call self%kill()
    if (present(cap_hint) .and. cap_hint > 0) then
      allocate(self%keys(cap_hint), self%cnts(cap_hint))
    else
      allocate(self%keys(8), self%cnts(8))
    end if
    self%keys = 0; self%cnts = 0; self%n = 0
  end subroutine sic_init

  subroutine sic_kill(self)
    class(smallIntCounts), intent(inout) :: self
    if (allocated(self%keys))  deallocate(self%keys)
    if (allocated(self%cnts))  deallocate(self%cnts)
    self%n = 0
  end subroutine sic_kill

  subroutine sic_grow(self)
    class(smallIntCounts), intent(inout) :: self
    integer(shortInt), allocatable :: k2(:), c2(:)
    integer :: newcap      ! default int for size arithmetic
    if (.not. allocated(self%keys)) then
      allocate(self%keys(8), self%cnts(8))
      self%keys = 0; self%cnts = 0; self%n = 0
      return
    end if
    newcap = max(8, 2*size(self%keys))
    allocate(k2(newcap), c2(newcap))
    k2 = 0; c2 = 0
    if (self%n > 0) then
      k2(1:self%n) = self%keys(1:self%n)
      c2(1:self%n) = self%cnts(1:self%n)
    end if
    call move_alloc(k2, self%keys)
    call move_alloc(c2, self%cnts)
  end subroutine sic_grow

  subroutine sic_add(self, key, delta)
    class(smallIntCounts), intent(inout) :: self
    integer(shortInt), intent(in) :: key
    integer(shortInt), intent(in), optional :: delta
    integer(shortInt) :: d
    integer :: i
    d = merge(delta, 1_shortInt, present(delta))
    if (.not. allocated(self%keys)) call self%init()
    do i = 1, self%n
      if (self%keys(i) == key) then
        self%cnts(i) = self%cnts(i) + d
        return
      end if
    end do
    if (self%n == size(self%keys)) call self%grow()
    self%n = self%n + 1_shortInt
    self%keys(self%n) = key
    self%cnts(self%n) = d
  end subroutine sic_add

  subroutine sic_sub(self, key, delta)
    class(smallIntCounts), intent(inout) :: self
    integer(shortInt), intent(in) :: key
    integer(shortInt), intent(in), optional :: delta
    integer(shortInt) :: d
    integer :: i
    d = merge(delta, 1_shortInt, present(delta))
    if (.not. allocated(self%keys)) return
    do i = 1, self%n
      if (self%keys(i) == key) then
        self%cnts(i) = self%cnts(i) - d
        if (self%cnts(i) <= 0_shortInt) then
          self%keys(i) = self%keys(self%n)
          self%cnts(i) = self%cnts(self%n)
          self%n = self%n - 1_shortInt
        end if
        return
      end if
    end do
  end subroutine sic_sub

  pure integer(shortInt) function sic_unique_count(self) result(nu)
    class(smallIntCounts), intent(in) :: self
    nu = self%n
  end function sic_unique_count

  function sic_to_list(self) result(out)
    class(smallIntCounts), intent(in) :: self
    integer(shortInt), allocatable :: out(:)
    if (self%n <= 0_shortInt) then
      allocate(out(0))
    else
      allocate(out(self%n))
      out = self%keys(1:self%n)
    end if
  end function sic_to_list

!==============================
! dynamic2dMatSet: user methods
!==============================

  subroutine init(self, n_slices, hint_unique)
    class(dynamic2dMatSet), intent(inout) :: self
    integer(shortInt), intent(in) :: n_slices
    integer(shortInt), intent(in), optional :: hint_unique

    call self%kill()

    if (n_slices < 0_shortInt) then
       call fatalError("dynamic2dMatSet:init", "n_slices < 0")
    end if

    allocate(self%slice(max(0, n_slices)))
    self%nused = 0_shortInt
    call self%gids%init(merge(hint_unique, 8_shortInt, present(hint_unique)))
  end subroutine init

  subroutine kill(self)
    class(dynamic2dMatSet), intent(inout) :: self
    integer :: i
    if (allocated(self%slice)) then
      do i = 1, size(self%slice)
        if (allocated(self%slice(i)%data)) deallocate(self%slice(i)%data)
        if (allocated(self%slice(i)%tags)) deallocate(self%slice(i)%tags)
      end do
      deallocate(self%slice)
    end if
    self%nused = 0_shortInt
    call self%gids%kill()
  end subroutine kill

  subroutine append(self, src, tags)
    class(dynamic2dMatSet), intent(inout) :: self
    real(defReal),     intent(in) :: src(:,:)
    integer(shortInt), intent(in) :: tags(:)
    integer :: ncols
    integer(shortInt) :: k
    integer :: j

    if (.not. allocated(self%slice)) then
       call fatalError("dynamic2dMatSet:append", "not initialised")
    end if
    if (size(src,1) /= 3) then
       call fatalError("dynamic2dMatSet:append", "first dimension must be 3")
    end if

    ncols = size(src,2)
    if (size(tags) /= ncols) then
       call fatalError("dynamic2dMatSet:append", "tags size mismatch with columns")
    end if
    if (self%nused >= size(self%slice)) then
       call fatalError("dynamic2dMatSet:append", "capacity exceeded")
    end if

    k = self%nused + 1_shortInt
    allocate(self%slice(k)%data, mold=src); self%slice(k)%data = src
    allocate(self%slice(k)%tags(ncols));    self%slice(k)%tags = tags
    self%nused = k

    do j = 1, ncols
      call self%gids%add(tags(j), 1_shortInt)
    end do
  end subroutine append

  pure function nslices(self) result(n)
    class(dynamic2dMatSet), intent(in) :: self
    integer(shortInt) :: n
    n = self%nused
  end function nslices

  subroutine scale(self, alpha)
    class(dynamic2dMatSet), intent(inout) :: self
    real(defReal), intent(in) :: alpha
    integer :: i, tcols

    if (self%nused == 0_shortInt) return
    if (alpha == 1.0_defReal) return

    tcols = 0
    do i = 1, self%nused
      if (allocated(self%slice(i)%data)) tcols = tcols + size(self%slice(i)%data, 2)
    end do

    if (alpha == 0.0_defReal) then
      do i = 1, self%nused
        if (allocated(self%slice(i)%data)) self%slice(i)%data = 0.0_defReal
      end do
      return
    end if

    if (tcols < 64) then
      do i = 1, self%nused
        if (allocated(self%slice(i)%data)) self%slice(i)%data = alpha * self%slice(i)%data
      end do
    else
!$omp parallel do default(none) private(i) shared(self,alpha) if(self%nused>1)
      do i = 1, self%nused
        if (allocated(self%slice(i)%data)) self%slice(i)%data = alpha * self%slice(i)%data
      end do
!$omp end parallel do
    end if
  end subroutine scale

  subroutine update_replace(self, k, src, tags)
    class(dynamic2dMatSet), intent(inout) :: self
    integer(shortInt), intent(in) :: k
    real(defReal),     intent(in) :: src(:,:)
    integer(shortInt), intent(in) :: tags(:)

    integer :: ncols, j

    if (k < 1_shortInt .or. k > self%nused) then
       call fatalError("dynamic2dMatSet:update_replace", "k out of range")
    end if
    if (size(src,1) /= 3) then
       call fatalError("dynamic2dMatSet:update_replace", "first dim must be 3")
    end if

    ncols = size(src,2)
    if (size(tags) /= ncols) then
       call fatalError("dynamic2dMatSet:update_replace", "tags size mismatch")
    end if

    ! Remove old tag counts
    if (allocated(self%slice(k)%tags)) then
      do j = 1, size(self%slice(k)%tags)
        call self%gids%sub(self%slice(k)%tags(j), 1_shortInt)
      end do
    end if

    ! Resize data if needed
    if (allocated(self%slice(k)%data)) then
      if (any(shape(self%slice(k)%data) /= shape(src))) then
        deallocate(self%slice(k)%data)
        allocate(self%slice(k)%data, mold=src)
      end if
    else
      allocate(self%slice(k)%data, mold=src)
    end if
    self%slice(k)%data = src

    ! Reuse tags allocation when possible
    if (allocated(self%slice(k)%tags)) then
      if (size(self%slice(k)%tags) /= ncols) then
        deallocate(self%slice(k)%tags)
        allocate(self%slice(k)%tags(ncols))
      end if
    else
      allocate(self%slice(k)%tags(ncols))
    end if
    self%slice(k)%tags = tags

    do j = 1, ncols
      call self%gids%add(tags(j), 1_shortInt)
    end do
  end subroutine update_replace

  subroutine filter_columns(self, k, keep)
    class(dynamic2dMatSet), intent(inout) :: self
    integer(shortInt), intent(in) :: k
    logical, intent(in) :: keep(:)

    integer :: ncols, y, i, j
    real(defReal),     allocatable :: tmp(:,:)
    integer(shortInt), allocatable :: ttmp(:)

    if (k < 1_shortInt .or. k > self%nused) then
       call fatalError("dynamic2dMatSet:filter_columns", "k out of range")
    end if
    if (.not. allocated(self%slice(k)%data)) return

    ncols = size(self%slice(k)%data, 2)
    if (size(keep) /= ncols) then
       call fatalError("dynamic2dMatSet:filter_columns", "keep size mismatch")
    end if

    y = count(keep)
    if (y == ncols) return
    if (y == 0) then
      do i = 1, ncols
        call self%gids%sub(self%slice(k)%tags(i), 1_shortInt)
      end do
      call self%delete(k)
      return
    end if

    do i = 1, ncols
      if (.not. keep(i)) call self%gids%sub(self%slice(k)%tags(i), 1_shortInt)
    end do

    allocate(tmp(3, y))
    allocate(ttmp(y))
    j = 0
    do i = 1, ncols
      if (keep(i)) then
        j = j + 1
        tmp(:, j) = self%slice(k)%data(:, i)
        ttmp(j)   = self%slice(k)%tags(i)
      end if
    end do

    call move_alloc(tmp,  self%slice(k)%data)
    call move_alloc(ttmp, self%slice(k)%tags)
  end subroutine filter_columns

  subroutine delete(self, k)
    class(dynamic2dMatSet), intent(inout) :: self
    integer(shortInt), intent(in) :: k
    integer :: i, j

    if (k < 1_shortInt .or. k > self%nused) then
       call fatalError("dynamic2dMatSet:delete", "k out of range")
    end if

    if (allocated(self%slice(k)%tags)) then
      do j = 1, size(self%slice(k)%tags)
        call self%gids%sub(self%slice(k)%tags(j), 1_shortInt)
      end do
    end if

    if (allocated(self%slice(k)%data)) deallocate(self%slice(k)%data)
    if (allocated(self%slice(k)%tags)) deallocate(self%slice(k)%tags)

    do i = k, self%nused - 1_shortInt
      call move_alloc(self%slice(i+1)%data, self%slice(i)%data)
      call move_alloc(self%slice(i+1)%tags, self%slice(i)%tags)
    end do
    if (self%nused >= 1_shortInt) then
      if (allocated(self%slice(self%nused)%data)) deallocate(self%slice(self%nused)%data)
      if (allocated(self%slice(self%nused)%tags)) deallocate(self%slice(self%nused)%tags)
    end if
    self%nused = self%nused - 1_shortInt
  end subroutine delete

  subroutine delete_columns(self, k, cols)
    class(dynamic2dMatSet), intent(inout) :: self
    integer(shortInt), intent(in) :: k
    integer(shortInt), intent(in) :: cols(:)

    integer :: ncols, c
    logical, allocatable :: keep(:)

    if (k < 1_shortInt .or. k > self%nused) then
       call fatalError("dynamic2dMatSet:delete_columns", "k out of range")
    end if
    if (.not. allocated(self%slice(k)%data)) return
    if (size(cols) == 0) return

    ncols = size(self%slice(k)%data, 2)

    ! Fast path: single column delete (manual splice)
    if (size(cols) == 1) then
      c = cols(1)
      if (c < 1 .or. c > ncols) then
        call fatalError("dynamic2dMatSet:delete_columns", "column index out of range")
      end if
      call self%gids%sub(self%slice(k)%tags(c), 1_shortInt)

      if (ncols == 1) then
        ! delete whole slice (O(#slices))
        call self%delete(k)
        return
      end if

      ! splice out column c
      block
        real(defReal),     allocatable :: tmp(:,:)
        integer(shortInt), allocatable :: ttmp(:)
        allocate(tmp(3, ncols-1))
        allocate(ttmp(ncols-1))
        if (c > 1) then
          tmp(:, 1:c-1) = self%slice(k)%data(:, 1:c-1)
          ttmp( 1:c-1 ) = self%slice(k)%tags(   1:c-1 )
        end if
        if (c < ncols) then
          tmp(:, c:ncols-1) = self%slice(k)%data(:, c+1:ncols)
          ttmp( c:ncols-1 ) = self%slice(k)%tags(   c+1:ncols)
        end if
        call move_alloc(tmp,  self%slice(k)%data)
        call move_alloc(ttmp, self%slice(k)%tags)
      end block
      return
    end if

    ! General case: build keep mask and reuse filter_columns (handles gids + compaction)
    allocate(keep(ncols)); keep = .true.
    block
      integer :: i, m, last
      ! Deduplicate and validate while marking
      do i = 1, size(cols)
        c = cols(i)
        if (c < 1 .or. c > ncols) then
          call fatalError("dynamic2dMatSet:delete_columns", "column index out of range")
        end if
        if (keep(c)) keep(c) = .false.   ! ignore duplicates
      end do
    end block

    call self%filter_columns(k, keep)
    deallocate(keep)
  end subroutine delete_columns

  subroutine get_copy(self, k, out, tags)
    class(dynamic2dMatSet), intent(in)  :: self
    integer(shortInt),      intent(in)  :: k
    real(defReal),          allocatable, intent(out) :: out(:,:)
    integer(shortInt),      allocatable, intent(out) :: tags(:)

    if (k < 1_shortInt .or. k > self%nused) return
    if (.not. allocated(self%slice(k)%data)) return
    if (.not. allocated(self%slice(k)%tags)) return

    allocate(out,  mold=self%slice(k)%data);  out  = self%slice(k)%data
    allocate(tags, mold=self%slice(k)%tags);  tags = self%slice(k)%tags
  end subroutine get_copy

  pure integer(shortInt) function unique_count(self) result(n)
    class(dynamic2dMatSet), intent(in) :: self
    n = self%gids%unique_count()
  end function unique_count

  pure logical function is_singleton(self) result(tf)
    class(dynamic2dMatSet), intent(in) :: self
    tf = (self%gids%n == 1_shortInt)
  end function is_singleton

  function unique_list(self) result(out)
    class(dynamic2dMatSet), intent(in) :: self
    integer(shortInt), allocatable :: out(:)
    out = self%gids%to_list()
  end function unique_list

end module dynamic2dMatSet_class