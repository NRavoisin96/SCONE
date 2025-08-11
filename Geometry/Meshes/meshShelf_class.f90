module meshShelf_class
  
  use dictionary_class,  only : dictionary
  use genericProcedures, only : fatalError, linFind, numToChar
  use intMap_class,      only : intMap
  use meshFactory_func,  only : new_mesh_ptr
  use mesh_inter,        only : mesh
  use numPrecision
  use publicObjects,     only : coordData
  use universalVariables
  
  implicit none
  private
  
  !!
  !! Small, local container to store polymorphic meshes in a single array.
  !!
  !! Public members:
  !!   name -> Name of the mesh.
  !!   ptr  -> Pointer to the mesh.
  !!
  type, public :: meshBox
    character(:), allocatable :: name
    class(mesh), pointer      :: ptr => null()
  end type meshBox
  
  !!
  !! Stores space for multiple mesh geometries.
  !!
  !! Private members:
  !!   meshes -> Array to store polymorphic meshes.
  !!   idMap  -> intMap to store mesh ids.
  !!
  type, public :: meshShelf
    private
    type(meshBox), dimension(:), allocatable :: meshes
    type(intMap)                             :: idMap
  contains
    procedure :: findHostMesh
    procedure :: getMeshBox
    procedure :: getMeshId
    procedure :: getMeshIdx
    procedure :: getMeshPtr
    procedure :: getOverallBoundingBoxBounds
    procedure :: getSize
    procedure :: init
    procedure :: kill
  end type meshShelf

contains
  
  !! Subroutine 'init'
  !!
  !! Basic description:
  !!   Initialises the meshShelf.
  !!
  !! Arguments:
  !!   dict [in] -> Dictionary with all mesh definitions.
  !!
  !! Errors:
  !!   fatalError if multiple meshes have the same Id.
  !!
  subroutine init(self, dict)
    class(meshShelf), intent(inout)               :: self
    class(dictionary), intent(in)                 :: dict
    character(nameLen), dimension(:), allocatable :: names
    integer(shortInt)                             :: i, id, idx, nMeshes
    character(*), parameter                       :: Here = 'init (meshShelf_class.f90)'
    
    ! Get all keys to subdictionaries.
    call dict % keys(names, 'dict')
    ! Allocate space.
    nMeshes = size(names)
    allocate(self % meshes(nMeshes))
    ! Build mesh geometries.
    do i = 1, nMeshes
      self % meshes(i) % name = names(i)
      self % meshes(i) % ptr => new_mesh_ptr(dict % getDictPtr(names(i)))
      id = self % meshes(i) % ptr % getId()
      ! Add Id to the map detecting any conflicts.
      idx = self % idMap % getOrDefault(id, NOT_PRESENT)
      if (idx /= NOT_PRESENT) call fatalError(Here,'Mesh geometries '//trim(names(i))//' and '//&
                                              trim(self % meshes(idx) % name)//' have the same id: '&
                                              //numToChar(id)//'.')
      call self % idMap % add(id, i)

    end do

  end subroutine init
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  subroutine kill(self)
    class(meshShelf), intent(inout) :: self
    integer(shortInt)               :: i
    
    ! Clear everything.
    if (allocated(self % meshes)) then
        do i = 1, size(self % meshes)
            call self % meshes(i) % ptr % kill()
            deallocate(self % meshes(i) % ptr)
            if (allocated(self % meshes(i) % name)) deallocate(self % meshes(i) % name)
            
        end do
        deallocate(self % meshes)
        
    end if
    call self % idMap % kill()

  end subroutine kill

  !!
  !!
  !!
  subroutine findHostMesh(self, data)
    class(meshShelf), intent(in)   :: self
    type(coordData), intent(inout) :: data
    integer(shortInt)              :: i
    character(*), parameter        :: here = 'findHostMesh (meshShelf_class.f90)'

    ! Check if shelf is allocated.
    if (.not. allocated(self % meshes)) call fatalError(here, 'Unallocated shelf.')
    do i = 1, size(self % meshes)
      call self % meshes(i) % ptr % findHostElement(data)
      if (0 < data % elementIdx) then
        data % meshIdx = i
        return

      end if

    end do

  end subroutine findHostMesh

  !!
  !!
  !!
  function getMeshBox(self, idx) result(box)
    class(meshShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    type(meshBox)                 :: box

    box = self % meshes(idx)

  end function getMeshBox
  
  !! Function 'getId'
  !!
  !! Basic description:
  !!   Returns the Id of the mesh with index idx.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the mesh.
  !!
  !! Result:
  !!   id -> Id of the mesh.
  !!
  !! Errors:
  !!   fatalError if idx is invalid.
  !!
  function getMeshId(self, idx) result(id)
    class(meshShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    integer(shortInt)             :: id, size
    character(*), parameter       :: Here = 'getId (meshShelf_class.f90)'
    
    ! Catch invalid idx.
    size = self % getSize()
    if (idx < 1 .or. idx > size) call fatalError(Here, 'Requested index: '//numToChar(idx)//' is not valid. Must be between &
                                                 &1 and '//numToChar(size))
    id = self % meshes(idx) % ptr % getId()

  end function getMeshId
  
  !! Function 'getIdx'
  !!
  !! Basic description:
  !!   Returns the index of the mesh with Id id.
  !!
  !! Arguments:
  !!   id [in] -> Id of the mesh.
  !!
  !! Result:
  !!   idx -> Index of the mesh.
  !!
  !! Errors:
  !!   fatalError if id is invalid.
  !!
  function getMeshIdx(self, id) result(idx)
    class(meshShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: id
    integer(shortInt)             :: idx
    integer(shortInt), parameter  :: NOT_PRESENT = -7
    character(*), parameter       :: Here = 'getIdx (meshShelf_class.f90)'
    
    idx = self % idMap % getOrDefault(id, NOT_PRESENT)
    if (idx == NOT_PRESENT) call fatalError(Here, 'There is no mesh with Id: '//numToChar(id)//'.')

  end function getMeshIdx
  
  !! Function 'getPtr'
  !!
  !! Basic description:
  !!   Returns a pointer to a mesh with index idx.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the mesh.
  !!
  !! Result:
  !!   ptr -> Pointer to the mesh.
  !!
  !! Errors:
  !!   fatalError if idx is invalid.
  !!
  function getMeshPtr(self, idx) result(ptr)
    class(meshShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    class(mesh), pointer          :: ptr
    integer(shortInt)             :: nMeshes
    character(*), parameter       :: Here = 'getPtr (meshShelf_class.f90)'
    
    ! Catch invalid idx.
    nMeshes = self % getSize()
    if (idx < 1 .or. nMeshes < idx) &
    call fatalError(Here, 'Index: '//numToChar(idx)//' is not valid. Must be between 1 and '//numToChar(nMeshes)//'.')
    
    ! Return pointer.
    ptr => self % meshes(idx) % ptr

  end function getMeshPtr
  
  !! Function 'getSize'
  !!
  !! Basic description:
  !!   Returns the size of the shelf.
  !!
  !! Result:
  !!   size -> Size of the shelf.
  !!
  !! Errors:
  !!   fatalError if shelf is unallocated.
  !!
  function getSize(self) result(nMeshes)
    class(meshShelf), intent(in) :: self
    integer(shortInt)            :: nMeshes
    character(*), parameter      :: Here = 'getSize (meshShelf_class.f90)'

    ! Check allocation and return size.
    if (.not. allocated(self % meshes)) call fatalError(Here, 'Requested size of unallocated shelf.')
    nMeshes = size(self % meshes)

  end function getSize

  !!
  !!
  !!
  function getOverallBoundingBoxBounds(self) result(overallBoundingBoxBounds)
    class(meshShelf), intent(in)   :: self
    integer(shortInt)              :: i
    real(defReal), dimension(3, 2) :: meshBoundingBoxBounds, overallBoundingBoxBounds
    character(*), parameter        :: here = 'getOverallBoundingBoxBounds (meshShelf_class.f90)'

    overallBoundingBoxBounds = reshape([INF, INF, INF, -INF, -INF, -INF], [3, 2])
    if (.not. allocated(self % meshes)) call fatalError(here, 'Shelf is unallocated.')
    do i = 1, size(self % meshes)
      meshBoundingBoxBounds = self % meshes(i) % ptr % getBoundingBoxBounds()
      overallBoundingBoxBounds(:, 1) = min(overallBoundingBoxBounds(:, 1), meshBoundingBoxBounds(:, 1))
      overallBoundingBoxBounds(:, 2) = max(overallBoundingBoxBounds(:, 2), meshBoundingBoxBounds(:, 2))

    end do

  end function getOverallBoundingBoxBounds

end module meshShelf_class
