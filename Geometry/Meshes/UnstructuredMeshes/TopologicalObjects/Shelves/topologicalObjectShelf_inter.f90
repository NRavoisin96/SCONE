module topologicalObjectShelf_inter

  use genericProcedures,       only : fatalError, numToChar
  use numPrecision
  use topologicalObject_inter, only : topologicalObject, topologicalObjectBox

  implicit none
  private

  ! Extendable procedures.
  public :: kill

  !!
  !!
  !!
  type, public, abstract :: topologicalObjectShelf
    private
    type(topologicalObjectBox), dimension(:), allocatable :: shelf
    integer(shortInt)                                     :: nObjects = 0
  contains
    procedure          :: addObject
    procedure          :: allocateShelf
    generic            :: getObjectBox => getObjectBox_shortInt, getObjectBox_shortIntArray
    procedure, private :: getObjectBox_shortInt
    procedure, private :: getObjectBox_shortIntArray
    procedure          :: getObjectsNumber
    procedure          :: getSize
    procedure          :: kill
    procedure          :: shrink
  end type topologicalObjectShelf

contains
  !!
  !!
  !!
  subroutine addObject(self, objectPtr)
    class(topologicalObjectShelf), intent(inout)          :: self
    class(topologicalObject), pointer, intent(in)         :: objectPtr
    integer(shortInt)                                     :: currentSize
    type(topologicalObjectBox), dimension(:), allocatable :: tempShelf

    ! First check if the shelf is allocated.
    if (allocated(self % shelf)) then
      ! Check if shelf is full and double its size if so.
      currentSize = size(self % shelf)
      if (self % nObjects == currentSize) then
        allocate(tempShelf(2 * currentSize))
        tempShelf(1:currentSize) = self % shelf
        call move_alloc(tempShelf, self % shelf)

      end if

    else
      ! Allocate shelf with a reasonable initial size.
      allocate(self % shelf(8))

    end if
    ! Update self % nObjects and add the new item at the next available position.
    self % nObjects = self % nObjects + 1
    self % shelf(self % nObjects) % ptr => objectPtr

  end subroutine addObject

  !!
  !!
  !!
  subroutine allocateShelf(self, nObjects)
    class(topologicalObjectShelf), intent(inout) :: self
    integer(shortInt), intent(in)                :: nObjects
    character(*), parameter                      :: here = 'allocateShelf (topologicalObjectShelf_inter.f90)'
    
    if (nObjects < 1) call fatalError(here, 'Size of shelf must be positive.')
    allocate(self % shelf(nObjects))

  end subroutine allocateShelf

  !!
  !!
  !!
  function getObjectBox_shortInt(self, idx) result(box)
    class(topologicalObjectShelf), intent(in) :: self
    integer(shortInt), intent(in)             :: idx
    type(topologicalObjectBox)                :: box

    box = self % shelf(idx)

  end function getObjectBox_shortInt

  !!
  !!
  !!
  function getObjectBox_shortIntArray(self, idxs) result(boxes)
    class(topologicalObjectShelf), intent(in)         :: self
    integer(shortInt), dimension(:), intent(in)       :: idxs
    type(topologicalObjectBox), dimension(size(idxs)) :: boxes
    integer(shortInt)                                 :: i

    do i = 1, size(idxs)
      boxes(i) = self % shelf(idxs(i))

    end do

  end function getObjectBox_shortIntArray

  !!
  !!
  !!
  elemental function getObjectsNumber(self) result(nObjects)
    class(topologicalObjectShelf), intent(in) :: self
    integer(shortInt)                         :: nObjects

    nObjects = self % nObjects

  end function getObjectsNumber

  !!
  !!
  !!
  elemental function getSize(self) result(shelfSize)
    class(topologicalObjectShelf), intent(in) :: self
    integer(shortInt)                         :: shelfSize

    shelfSize = size(self % shelf)

  end function getSize

  !!
  !!
  !!
  subroutine kill(self)
    class(topologicalObjectShelf), intent(inout) :: self
    integer(shortInt)                            :: i

    ! Local.
    if (allocated(self % shelf)) then
      do i = 1, size(self % shelf)
        call self % shelf(i) % ptr % kill()
        deallocate(self % shelf(i) % ptr)

      end do
      deallocate(self % shelf)

    end if
    self % nObjects = 0

  end subroutine kill

  !!
  !!
  !!
  subroutine shrink(self)
    class(topologicalObjectShelf), intent(inout)          :: self
    type(topologicalObjectBox), dimension(:), allocatable :: tempShelf
    character(*), parameter                               :: here = 'shrinkShelf (topologicalObjectShelf_inter.f90)'

    ! First check if shelf is unallocated.
    if (.not. allocated(self % shelf)) call fatalError(here, 'Attempting to shrink an unallocated shelf.')
    allocate(tempShelf(self % nObjects))
    tempShelf = self % shelf(1:self % nObjects)
    call move_alloc(tempShelf, self % shelf)

  end subroutine shrink

end module topologicalObjectShelf_inter