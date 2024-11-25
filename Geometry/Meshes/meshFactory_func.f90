module meshFactory_func
  
  use numPrecision
  use genericProcedures,  only : fatalError
  use dictionary_class,   only : dictionary
  
  ! Interface mesh.
  use mesh_inter,         only : mesh
  
  ! Meshes.
  use OpenFOAMMesh_class, only : OpenFOAMMesh
  
  implicit none
  private
  
  ! ** ADD NAME OF NEW MESH TO THE LIST **!
  ! List that contains acceptable types of meshes
  ! NOTE: It is necessary to adjust trailing blanks so all entries have the same length
  character(nameLen), dimension(*), parameter :: AVAILABLE_MESHES = ['OpenFOAMMesh']
  ! Public interface.
  public :: new_mesh, new_mesh_ptr

  contains
  
  !! Function 'new_mesh_ptr'
  !!
  !! Basic description:
  !!   Returns a pointer to a new instance of an allocated mesh.
  !!
  !! Arguments:
  !!   dict [in] -> Dictionary with mesh definition.
  !!   name [in] -> Name of the mesh.
  !!
  !! Result:
  !!   ptr -> Pointer to the allocated mesh.
  !!
  !! Errors:
  !!   - fatalError if mesh folder path does not exist;
  !!   - fatalError if type of mesh is unknown.
  !!
  function new_mesh_ptr(dict) result(new)
    class(dictionary), intent(in) :: dict
    class(mesh), pointer          :: new
    character(nameLen)            :: type
    character(pathLen)            :: path
    character(:), allocatable     :: trimmedPath
    logical(defBool)              :: pathExists
    character(100), parameter     :: Here = 'new_mesh_ptr (meshFactory_func.f90)'
    ! Retrieve type of the mesh.
    call dict % get(type, 'type')
    ! Retrieve path to the mesh folder.
    call dict % get(path, 'path')
    ! Check that the folder corresponding to path exists.
    trimmedPath = trim(path)
    inquire(file = trimmedPath, exist = pathExists)
    ! If the provided mesh folder path does not exist call fatalError.
    if (.not. pathExists) call fatalError(Here, 'The provided mesh folder path does not exist.')
    ! Allocate appropriate mesh.
    ! ** FOR NEW MESH ADD CASE STATEMENT HERE ** !
    select case(type)
      case('OpenFOAMMesh')
        allocate(OpenFOAMMesh :: new)
      case default
        print '(A)', 'AVAILABLE MESHES: '
        print '(A)', AVAILABLE_MESHES
        call fatalError(Here, 'Unrecognised mesh type: '//trim(type)//'.')
    end select
    ! Initialise the mesh geometry.
    call new % init(trimmedPath, dict)
  end function new_mesh_ptr
  
  !! Subroutine 'new_mesh'
  !!
  !! Basic description:
  !!   Allocates an allocatable mesh.
  !!
  !! Arguments:
  !!   new [out] -> Mesh to be allocated.
  !!   dict [in] -> Dictionary with the mesh definition.
  !!
  subroutine new_mesh(new, dict)
    class(mesh), allocatable, intent(out) :: new
    class(dictionary), intent(in)         :: dict
    class(mesh), pointer                  :: temp
    ! Obtain pointer to temporary mesh.
    temp => new_mesh_ptr(dict)
    ! Now copy the temporary mesh into the new mesh.
    allocate(new, source = temp)
    ! Clear the temporary mesh and deallocate pointer to prevent memory leaks.
    call temp % kill()
    deallocate(temp)
  end subroutine new_mesh
end module meshFactory_func