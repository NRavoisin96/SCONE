module universeShelf_class

  use numPrecision
  use universalVariables,   only : NOT_PRESENT
  use genericProcedures,    only : fatalError, numToChar
  use dictionary_class,     only : dictionary
  use intMap_class,         only : intMap
  use charMap_class,        only : charMap
  use surfaceShelf_class,   only : surfaceShelf
  use cellShelf_class,      only : cellShelf
  use meshShelf_class,      only : meshShelf

  ! Universe objects
  use universe_inter,       only : universe
  use universeFactory_func, only : new_universe_ptr
  use uniFills_class,       only : uniFills

  implicit none
  private

  !!
  !! Small, local container to store polymorphic universes in an array
  !!
  !! Public members:
  !!   name -> Name of the universe
  !!   ptr  -> Pointer to the universe
  !!
  type :: uniBox
    character(:), allocatable :: name
    class(universe), pointer  :: ptr => null()
  end type uniBox

  !!
  !! Storage space for universes defined in the geometry
  !!
  !! Sample Dictionary Input:
  !!   universes {
  !!     uni1 { <universe definition>}
  !!     uni2 { <universe definition>}
  !!     ...
  !!   }
  !!
  !! Private Members:
  !!   unis -> Array with pointers to different universes
  !!   idMap -> Map between uniId and uniIdx
  !!
  !! Interface:
  !!   init    -> Initialise and build uniFills
  !!   getPtr  -> Get pointer to a universe given by its index
  !!   getPtr_fast -> Get pointer to a universe without bounds checking. Should be used in
  !!     speed-critical parts.
  !!   getIdx  -> Get uniIdx of a universe given by uniId
  !!   getId   -> Get uniId of a universe given by uniIdx
  !!   getSize -> Return the number of universes (max uniIdx)
  !!   kill    -> Return to uninitialised state
  !!
  !! NOTE: Because universes are stored as pointers, calling `kill` is crucial
  !!   to prevent memory leaks. TODO: Add `final` procedure here?
  !!
  type, public :: universeShelf
    private
    type(uniBox), dimension(:), allocatable :: unis
    type(intMap)                            :: idMap
  contains
    procedure :: init
    procedure :: getPtr
    procedure :: getPtr_fast
    procedure :: getIdx
    procedure :: getId
    procedure :: getSize
    procedure :: kill
  end type universeShelf

contains

  !!
  !! Initialise universeShelf
  !!
  !! Builds all definitions and load fill information into uniFills
  !!
  !! Args:
  !!   fills [out]   -> `uniFills` that contains fill information for every universe
  !!   dict [in]     -> Dictionary with universe definitions
  !!   cells [inout] -> `cellShelf` with used defined cells
  !!   surfs [inout] -> `surfaceShelf` with user defined surfaces
  !!   mats [in]     -> Map of material names to matIdx
  !!
  !! Errors:
  !!   fatalError if there are multiple universes with the same id
  !!
  subroutine init(self, dict, mats, fills, cells, surfs, meshes)
    class(universeShelf), intent(inout)           :: self
    class(dictionary), intent(in)                 :: dict
    type(charMap), intent(in)                     :: mats
    type(uniFills), intent(out)                   :: fills
    type(cellShelf), intent(inout)                :: cells
    type(surfaceShelf), intent(inout)             :: surfs
    type(meshShelf), intent(inout)                :: meshes
    character(nameLen), dimension(:), allocatable :: names
    integer(shortInt)                             :: nUniverses, i, id, idx
    character(:), allocatable                     :: name
    integer(shortInt), dimension(:), allocatable  :: fillInfo
    character(100), parameter                     :: Here = 'init (universeShelf_class.f90)'

    ! Get all universe names.
    call dict % keys(names, 'dict')

    ! Allocate space and uniFills size.
    nUniverses = size(names)
    allocate (self % unis(nUniverses))
    call fills % init(nUniverses)

    ! Build universes.
    do i = 1, nUniverses
      ! Retrieve name of current universe and build it.
      name = trim(names(i))
      self % unis(i) % name = name

      ! Build new interface
      call new_universe_ptr(dict % getDictPtr(name), mats, cells, surfs, meshes, fillInfo, self % unis(i) % ptr)

      ! Add ID to map detecting any conflicts
      id = self % unis(i) % ptr % id()
      idx = self % idMap % getOrDefault(id, NOT_PRESENT)
      if (idx /= NOT_PRESENT) call fatalError(Here, 'Universes '//name// ' & '//trim(self % unis(idx) % name)//&
                                              ' have the same id: '//numToChar(id)//'.')
      
      call self % idMap % add(id, i)

      ! Store content info into fills and set new universe index.
      call fills % addUniverse(i, id, fillInfo)
      call self % unis(i) % ptr % setIdx(i)

    end do

    ! Finish build.
    call fills % finishBuild(self % idMap)

  end subroutine init

  !!
  !! Return pointer to the universe indicated by idx (FAST) without bound check
  !!
  !! Args:
  !!   idx [in] -> Index of the universe
  !!
  !! Result:
  !!   Pointer to the universe under idx.
  !!
  !! Errors:
  !!   Pointer will have undefined status if the idx is not valid (will point to some place
  !!   it is unlikely it will be null so associated procedure may return true!)
  !!
  function getPtr_fast(self, idx) result(ptr)
    class(universeShelf), intent(in) :: self
    integer(shortInt), intent(in)    :: idx
    class(universe), pointer         :: ptr

    ptr => self % unis(idx) % ptr

  end function getPtr_fast

  !!
  !! Return pointer to the universe indicated by idx (SLOW) with bound checking
  !!
  !! Slow access with bounds checking. To be used in non-performance
  !! critical areas
  !!
  !! Args:
  !!   idx [in] -> Index of the universe
  !!
  !! Result:
  !!   Pointer to the universe under index idx
  !!
  !! Errors:
  !!   fatalError is idx does not correspond to a universe (is out-of-bounds)
  !!
  function getPtr(self, idx) result (ptr)
    class(universeShelf), intent(in) :: self
    integer(shortInt), intent(in)    :: idx
    class(universe), pointer         :: ptr
    integer(shortInt)                :: nUniverses
    character(100), parameter        :: Here = 'getPtr (universeShelf_class.f90)'

    ! Return pointer
    nUniverses = size(self % unis)
    if (idx < 1 .or. idx > nUniverses) call fatalError(Here, 'Requested index: '//numToChar(idx)//' is not valid. Must be between &
                                                       &1 and '//numToChar(nUniverses)//'.')
    
    ptr => self % unis(idx) % ptr

  end function getPtr

  !!
  !! Return IDX of a universe with the ID
  !!
  !! Args:
  !!   id [in] -> Id of the universe
  !!
  !! Result:
  !!   Index of the universe with the given Id
  !!
  !! Error:
  !!   fatalError if there is no universe with the given Id
  !!
  function getIdx(self, id) result(idx)
    class(universeShelf), intent(in) :: self
    integer(shortInt), intent(in)    :: id
    integer(shortInt)                :: idx
    character(100), parameter        :: Here = 'getIdx (universeShelf_class.f90)'

    idx = self % idMap % getOrDefault(id, NOT_PRESENT)

    if (idx == NOT_PRESENT) call fatalError(Here, 'There is no universe with id: '//numToChar(id)//'.')

  end function getIdx

  !!
  !! Return ID of the universe given by index
  !!
  !! Args:
  !!   idx [in] -> Index of the universe
  !!
  !! Result:
  !!   ID of the universe under the index
  !!
  !! Error:
  !!   fatalError is idx does not correspond to a universe (is out-of-bounds)
  !!
  function getId(self, idx) result(id)
    class(universeShelf), intent(in) :: self
    integer(shortInt), intent(in)    :: idx
    integer(shortInt)                :: nUniverses, id
    character(100), parameter        :: Here = 'getId (universeShelf_class.f90)'

    ! Catch invalid idx
    nUniverses = size(self % unis)
    if (idx < 1 .or. idx > nUniverses) call fatalError(Here, 'Requested index: '//numToChar(idx)//' is not valid. Must be between &
                                                       &1 and '//numToChar(nUniverses)//'.')

    id = self % unis(idx) % ptr % id()

  end function getId

  !!
  !! Return size of the shelf
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Number of universes on the shelf
  !!
  elemental function getSize(self) result(N)
    class(universeShelf), intent(in) :: self
    integer(shortInt)                :: N

    N = size(self % unis)

  end function getSize

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(universeShelf), intent(inout) :: self
    integer(shortInt)                   :: i

    if (allocated(self % unis)) then
        do i = 1, size(self % unis)
            call self % unis(i) % ptr % kill()
            if (allocated(self % unis(i) % name)) deallocate(self % unis(i) % name)
        
        end do
        deallocate(self % unis)
    end if
    call self % idMap % kill()

  end subroutine kill

end module universeShelf_class
