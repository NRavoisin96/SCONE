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
    generic            :: getObjectBoundingBoxBounds => getObjectBoundingBoxBounds_shortInt, &
                                                        getObjectBoundingBoxBounds_shortIntArray
    procedure, private :: getObjectBoundingBoxBounds_shortInt
    procedure, private :: getObjectBoundingBoxBounds_shortIntArray
    generic            :: getObjectCentroid => getObjectCentroid_shortInt, getObjectCentroid_shortIntArray
    procedure, private :: getObjectCentroid_shortInt
    procedure, private :: getObjectCentroid_shortIntArray
    generic            :: getObjectElements => getObjectElements_shortInt, getObjectElements_shortIntArray
    procedure, private :: getObjectElements_shortInt
    procedure, private :: getObjectElements_shortIntArray
    procedure          :: getObjectsNumber
    procedure          :: getShelf
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
  function getObjectBoundingBoxBounds_shortInt(self, idx) result(bounds)
    class(topologicalObjectShelf), intent(in) :: self
    integer(shortInt), intent(in)             :: idx
    real(defReal), dimension(3, 2)            :: bounds

    bounds = self % shelf(idx) % ptr % getBoundingBoxBounds()

  end function getObjectBoundingBoxBounds_shortInt

  !!
  !!
  !!
  function getObjectBoundingBoxBounds_shortIntArray(self, idxs) result(bounds)
    class(topologicalObjectShelf), intent(in)   :: self
    integer(shortInt), dimension(:), intent(in) :: idxs
    real(defReal), dimension(3, 2 * size(idxs)) :: bounds
    integer(shortInt)                           :: i

    do i = 1, size(idxs)
      bounds(:, 2 * (i - 1) + 1:2 * i) = self % shelf(idxs(i)) % ptr % getBoundingBoxBounds()

    end do

  end function getObjectBoundingBoxBounds_shortIntArray

  !!
  !!
  !!
  function getObjectCentroid_shortInt(self, idx) result(centroid)
    class(topologicalObjectShelf), intent(in) :: self
    integer(shortInt), intent(in)             :: idx
    real(defReal), dimension(3)               :: centroid

    centroid = self % shelf(idx) % ptr % getCentroid()

  end function getObjectCentroid_shortInt

  !!
  !!
  !!
  function getObjectCentroid_shortIntArray(self, idxs) result(centroids)
    class(topologicalObjectShelf), intent(in)   :: self
    integer(shortInt), dimension(:), intent(in) :: idxs
    real(defReal), dimension(3, size(idxs))     :: centroids
    integer(shortInt)                           :: i

    do i = 1, size(idxs)
      centroids(:, i) = self % shelf(idxs(i)) % ptr % getCentroid()

    end do

  end function getObjectCentroid_shortIntArray

  !!
  !!
  !!
  function getShelf(self) result(shelf)
    class(topologicalObjectShelf), intent(in)             :: self
    type(topologicalObjectBox), dimension(:), allocatable :: shelf

    if (allocated(self % shelf)) then
      shelf = self % shelf

    else
      allocate(shelf(0))

    end if

  end function getShelf

  !!
  !!
  !!
  function getObjectElements_shortInt(self, idx) result(elements)
    class(topologicalObjectShelf), intent(in)             :: self
    integer(shortInt), intent(in)                         :: idx
    type(topologicalObjectBox), dimension(:), allocatable :: elements

    elements = self % shelf(idx) % ptr % getElements()

  end function getObjectElements_shortInt

  !!
  !!
  !!
  function getObjectElements_shortIntArray(self, idxs) result(elements)
    class(topologicalObjectShelf), intent(in)             :: self
    integer(shortInt), dimension(:), intent(in)           :: idxs
    type(topologicalObjectBox), dimension(:), allocatable :: elements, objectElements, tempElements
    type(topologicalObjectBox), dimension(size(idxs))     :: boxes
    integer(shortInt)                                     :: i, j, k, nElements
    logical(defBool)                                      :: alreadyFound

    boxes = self % getObjectBox(idxs)
    do i = 1, size(idxs)
      objectElements = boxes(i) % ptr % getElements()
      do j = 1, size(objectElements)
        alreadyFound = .false.
        if (allocated(elements)) then
          do k = 1, size(elements)
            if (associated(objectElements(j) % ptr, elements(k) % ptr)) then
              alreadyFound = .true.
              exit

            end if

          end do

        end if

        if (.not. alreadyFound) then
          if (allocated(elements)) then
            nElements = size(elements)
            allocate(tempElements(nElements + 1))
            tempElements(1:nElements) = elements
            tempElements(nElements + 1) = objectElements(j)
            call move_alloc(tempElements, elements)

          else
            allocate(elements(1))
            elements(1) = objectElements(j)

          end if

        end if

      end do

    end do

  end function getObjectElements_shortIntArray

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