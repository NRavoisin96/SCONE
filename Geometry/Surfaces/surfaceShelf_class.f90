module surfaceShelf_class

  use numPrecision
  use genericProcedures,   only : fatalError, numToChar
  use dictionary_class,    only : dictionary
  use intMap_class,        only : intMap
  use surface_inter,       only : surface
  use surfaceFactory_func, only : new_surface_ptr

  implicit none
  private

  !!
  !! Small, local container to store polymorphic surfaces in a single array
  !!
  !! Public members:
  !!   name -> Name of the surface
  !!   ptr  -> Pointer to the surface
  !!
  type :: surfaceBox
    character(:), allocatable :: name
    class(surface), pointer   :: ptr => null()
  end type surfaceBox

  !!
  !! Storage space for surfaces defined in the geometry
  !!
  !! Sample dictionary input:
  !!   surfaces {
  !!    surf1 { <surface definition> }
  !!    surf2 { <surface definition> }
  !!    surf3 { <surface definition> }
  !!    ...
  !!    }
  !!
  !! Private Members:
  !!   surfaces -> Array to store pointers to polymorphic surfaces
  !!   idMap -> Map between surface ID and corresponding index
  !!
  !! Interface:
  !!   init    -> Initialise from a dictionary
  !!   getPtr  -> Return pointer to a surface given its index
  !!   getIdx  -> Return index of a surface given its id
  !!   getID   -> Return id of a surface given its idx
  !!   getSize -> Return the number of surfaces (max surfIdx)
  !!   kill -> Return to uninitialised state
  !!
  !! NOTE: Becouse surfaces are stored as pointers, calling `kill` is crutial to prevent
  !!   memory leaks. TODO: Add `final` procedure here ?
  !!
  type, public :: surfaceShelf
    private
    type(surfaceBox), dimension(:), allocatable :: surfaces
    type(intMap)                                :: idMap

  contains
    procedure :: init
    procedure :: getPtr
    procedure :: getIdx
    procedure :: getId
    procedure :: getSize
    procedure :: kill
  end type surfaceShelf

contains

  !!
  !! Loads surfaces into shelf
  !!
  !! Args:
  !!   dict [in] -> Dictionary with subdictionaries that contain surface definitions.
  !!
  !! Errors:
  !!   fatalError if there are clashes in surface ids.
  !!
  subroutine init(self, dict)
    class(surfaceShelf), intent(inout)            :: self
    class(dictionary), intent(in)                 :: dict
    character(nameLen), dimension(:), allocatable :: names
    character(:), allocatable                     :: trimmedName
    integer(shortInt)                             :: nSurfaces, i, id, idx
    integer(shortInt), parameter                  :: NOT_PRESENT = -7
    character(100), parameter                     :: Here = 'init (surfaceShelf_class.f90)'

    ! Get all keys for subdictionaries and compute number of surfaces to allocate.
    call dict % keys(names, 'dict')
    nSurfaces = size(names)

    ! Allocate space and build surfaces.
    allocate (self % surfaces(nSurfaces))
    do i = 1, nSurfaces
      trimmedName = trim(names(i))
      self % surfaces(i) % name = trimmedName
      self % surfaces(i) % ptr => new_surface_ptr(dict % getDictPtr(names(i)))
      id = self % surfaces(i) % ptr % getId()

      ! Add id to the map detecting any conflicts.
      idx = self % idMap % getOrDefault(id, NOT_PRESENT)
      if (idx /= NOT_PRESENT) call fatalError(Here,'Surfaces '//trimmedName// ' & '//&
                                              self % surfaces(idx) % name//&
                                              ' have the same Id: '//numToChar(id)//'.')

      call self % idMap % add(id, i)

    end do

  end subroutine init

  !!
  !! Returns pointer to the surface indicated by index
  !!
  !! Args:
  !!   idx [in] -> Index of the surface
  !!
  !! Result:
  !!   Pointer to a surface.
  !!
  !! Error:
  !!   fatalError if idx does not correspond to a surface (is out-of-bounds).
  !!
  function getPtr(self, idx) result (ptr)
    class(surfaceShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    integer(shortInt)               :: nSurfaces
    class(surface), pointer         :: ptr
    character(100), parameter       :: Here = 'getPtr (surfaceShelf_class.f90)'

    ! Catch invalid idx.
    nSurfaces = size(self % surfaces)
    if (idx < 1 .or. idx > nSurfaces) call fatalError(Here, 'Requested index: '//numToChar(idx)//' is not valid. Must be between &
                                                      &1 and '//numToChar(nSurfaces)//'.')

    ! Return pointer
    ptr => self % surfaces(idx) % ptr

  end function getPtr

  !!
  !! Returns index of a surface with id.
  !!
  !! Args:
  !!   id [in] -> Id of the surface
  !!
  !! Result:
  !!   Index of the surface.
  !!
  !! Error:
  !!   fatalError if there is no surface with id.
  !!
  function getIdx(self, id) result(idx)
    class(surfaceShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: id
    integer(shortInt)               :: idx
    integer(shortInt), parameter    :: NOT_PRESENT = -7
    character(100), parameter       :: Here = 'getIdx (surfaceShelf_class.f90)'

    idx = self % idMap % getOrDefault(id, NOT_PRESENT)
    if (idx == NOT_PRESENT) call fatalError(Here, 'There is no surface with Id: '//numToChar(id)//'.')

  end function getIdx

  !!
  !! Returns the id of a surface with index = idx.
  !!
  !! Args:
  !!   idx [in] -> Index of the surface.
  !!
  !! Result:
  !!   Id of the surface.
  !!
  !! Error:
  !!   fatalError if idx does not correspond to a surface (is out-of-bounds).
  !!
  function getId(self, idx) result(id)
    class(surfaceShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    integer(shortInt)               :: id, nSurfaces
    character(100), parameter       :: Here = 'getId (surfaceShelf_class.f90)'

    ! Catch invalid idx.
    nSurfaces = size(self % surfaces)
    if (idx < 1 .or. idx > nSurfaces) call fatalError(Here, 'Requested index: '//numToChar(idx)//' is not valid. Must be between &
                                                      &1 and '//numToChar(nSurfaces)//'.')

    id = self % surfaces(idx) % ptr % getId()

  end function getId

  !!
  !! Returns size of the shelf
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Number of surfaces in the shelf.
  !!
  elemental function getSize(self) result(nSurfaces)
    class(surfaceShelf), intent(in) :: self
    integer(shortInt)               :: nSurfaces

    nSurfaces = size(self % surfaces)

  end function getSize

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(surfaceShelf), intent(inout) :: self
    integer(shortInt)                  :: i

    if (allocated(self % surfaces)) then
      do i = 1, size(self % surfaces)
        call self % surfaces(i) % ptr % kill()
        if (allocated(self % surfaces(i) % name)) deallocate(self % surfaces(i) % name)

      end do
      deallocate(self % surfaces)
    end if

    call self % idMap % kill()

  end subroutine kill

end module surfaceShelf_class
