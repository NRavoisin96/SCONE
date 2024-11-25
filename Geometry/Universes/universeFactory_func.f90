module universeFactory_func

  use numPrecision
  use genericProcedures,  only : fatalError
  use dictionary_class,   only : dictionary
  use charMap_class,      only : charMap
  use surfaceShelf_class, only : surfaceShelf
  use cellShelf_class,    only : cellShelf
  use meshShelf_class,    only : meshShelf

  ! Universe interface
  use universe_inter,     only : universe

  ! Universes
  use rootUniverse_class, only : rootUniverse
  use cellUniverse_class, only : cellUniverse
  use pinUniverse_class,  only : pinUniverse
  use latUniverse_class,  only : latUniverse
  use meshUniverse_class, only : meshUniverse
  
  implicit none
  private

  ! List contains acceptable types of universe
  ! NOTE: It is necessary to adjust trailing blanks so all entries have the same length
  character(nameLen), dimension(*), parameter :: AVAILABLE_UNIVERSES = ['cellUniverse',&
                                                                        'latUniverse ',&
                                                                        'meshUniverse',&
                                                                        'pinUniverse ',&
                                                                        'rootUniverse']

  ! Public Interface
  public :: new_universe_ptr
  public :: new_universe

contains
  !!
  !! Point a pointer to a new instance of an allocated universe
  !!
  !! Args:
  !!  dict [in]      -> Dictionary with universe definition
  !!  mats [in]      -> Map of material names to matIdx
  !!  cells [inout]  -> Shelf with user-defined cells
  !!  surfs [inout]  -> Shelf with user-defined surfaces
  !!  meshes [inout] -> Shelf with user-defined meshes
  !!  fills [out]    -> Allocatable integer array with filling of different uniqueIds
  !!  ptr [out]      -> Pointer to the new universe
  !!
  !! Errors:
  !!   fatalError if type of universe is unknown
  !!
  subroutine new_universe_ptr(dict, mats, cells, surfs, meshes, fills, ptr)
    class(dictionary), intent(in)                             :: dict
    type(charMap), intent(in)                                 :: mats
    type(cellShelf), intent(inout)                            :: cells
    type(surfaceShelf), intent(inout)                         :: surfs
    type(meshShelf), intent(inout)                            :: meshes
    integer(shortInt), dimension(:), allocatable, intent(out) :: fills
    class(universe), pointer, intent(out)                     :: ptr
    character(nameLen)                                        :: type
    character(100), parameter                                 :: Here = 'new_universe_ptr (universeFactory_func.f90)'

    ! Obtain type of the universe
    call dict % get(type, 'type')

    ! Allocate appropriate universe
    select case (type)
      case ('rootUniverse')
        allocate(rootUniverse :: ptr)

      case ('cellUniverse')
        allocate(cellUniverse :: ptr)

      case ('pinUniverse')
        allocate(pinUniverse :: ptr)

      case ('latUniverse')
        allocate(latUniverse :: ptr)

      case ('meshUniverse')
        allocate(meshUniverse :: ptr)

      case default
        print '(A)', 'AVAILABLE UNIVERSES: '
        print '(A)', AVAILABLE_UNIVERSES
        call fatalError(Here, 'Unrecognised type of universe: '//trim(type)//'.')

    end select

    ! Initialise universe.
    call ptr % init(dict, mats, fills, cells, surfs, meshes)

  end subroutine new_universe_ptr

  !!
  !! Allocate an allocatable universe
  !!
  !! Args:
  !!  dict [in]      -> Dictionary with universe definition
  !!  mats [in]      -> Map of material names to matIdx
  !!  cells [inout]  -> Shelf with user-defined cells
  !!  surfs [inout]  -> Shelf with user-defined surfaces
  !!  meshes [inout] -> Shelf with user-defined meshes
  !!  fills [out]    -> Allocatable integer array with filling of different uniqueID
  !!  new [out]      -> Universe to be allocated
  !!
  !! Errors:
  !!   fatalError if type of universe is unknown
  !!
  subroutine new_universe(dict, mats, cells, surfs, meshes, fills, new)
    class(dictionary), intent(in)                             :: dict
    type(charMap), intent(in)                                 :: mats
    type(cellShelf), intent(inout)                            :: cells
    type(surfaceShelf), intent(inout)                         :: surfs
    type(meshShelf), intent(inout)                            :: meshes
    integer(shortInt), dimension(:), allocatable, intent(out) :: fills
    class(universe), allocatable, intent(out)                 :: new
    class(universe), pointer                                  :: temp

    ! Allocate temporary
    call new_universe_ptr(dict, mats, cells, surfs, meshes, fills, temp)
    allocate(new, source = temp)

    ! Clean up
    call temp % kill()
    deallocate(temp)

  end subroutine new_universe

end module universeFactory_func