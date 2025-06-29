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
  contains
    procedure          :: addObject
    procedure          :: allocateShelf
    procedure          :: expandShelf
    generic            :: getObjectBox => getObjectBox_shortInt, getObjectBox_shortIntArray
    procedure, private :: getObjectBox_shortInt
    procedure, private :: getObjectBox_shortIntArray
    procedure          :: getSize
    procedure          :: kill
    procedure          :: shrinkShelf
  end type topologicalObjectShelf

contains
  !!
  !!
  !!
  subroutine addObject(self, idx, objectPtr)
    class(topologicalObjectShelf), intent(inout)  :: self
    integer(shortInt), intent(in)                 :: idx
    class(topologicalObject), pointer, intent(in) :: objectPtr
    character(*), parameter                       :: here = 'addObject (topologicalObjectShelf_inter.f90)'

    ! Check that idx is valid.
    if (.not. allocated(self % shelf)) call fatalError(here, 'Trying to add pointer to unallocated shelf.')
    if (idx < 1 .or. size(self % shelf) < idx) call fatalError(here, 'Invalid index: '//numToChar(idx)//'.')
    self % shelf(idx) % ptr => objectPtr

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

  !! Subroutine 'expandShelf'
  !!
  !! Basic description:
  !!   Expands the shelf by a specified number of additional vertices. Copies elements
  !!   already present. Allocates the shelf if it is not allocated yet.
  !!
  !! Arguments:
  !!   nAdditionalVertices [in] -> Number of additional vertices to be included in the shelf.
  !!
  subroutine expandShelf(self, nAdditionalObjects)
    class(topologicalObjectShelf), intent(inout)          :: self
    integer(shortInt), intent(in)                         :: nAdditionalObjects
    integer(shortInt)                                     :: oldSize
    type(topologicalObjectBox), dimension(:), allocatable :: tempShelf

    if (nAdditionalObjects < 1) return

    if (allocated(self % shelf)) then
      oldSize = size(self % shelf)
      allocate(tempShelf(oldSize + nAdditionalObjects))
      tempShelf(1:oldSize) = self % shelf
      call move_alloc(tempShelf, self % shelf)

    else
      allocate(self % shelf(nAdditionalObjects))

    end if

  end subroutine expandShelf

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
  elemental function getSize(self) result(shelfSize)
    class(topologicalObjectShelf), intent(in) :: self
    integer(shortInt)                         :: shelfSize

    shelfSize = size(self % shelf)

  end function getSize

  !!
  !!
  !!
  elemental subroutine kill(self)
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

  end subroutine kill

  !!
  !!
  !!
  subroutine shrinkShelf(self, newSize)
    class(topologicalObjectShelf), intent(inout)          :: self
    integer(shortInt), intent(in)                         :: newSize
    integer(shortInt)                                     :: oldSize
    type(topologicalObjectBox), dimension(:), allocatable :: tempShelf
    character(*), parameter                               :: here = 'shrinkShelf (topologicalObjectShelf_inter.f90)'

    if (newSize < 1) return
    if (allocated(self % shelf)) then
      oldSize = size(self % shelf)
      if (oldSize < newSize) call fatalError(here, 'Attempting to shrink shelf to a size larger than its original size.')
      allocate(tempShelf(newSize))
      tempShelf = self % shelf(1:newSize)
      call move_alloc(tempShelf, self % shelf)

    else
      allocate(self % shelf(newSize))

    end if

  end subroutine shrinkShelf

end module topologicalObjectShelf_inter