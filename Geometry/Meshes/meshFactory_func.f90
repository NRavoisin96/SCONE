module meshFactory_func
  
  use numPrecision
  use genericProcedures,     only : fatalError
  use dictionary_class,      only : dictionary
  
  ! Interface mesh.
  use mesh_inter,            only : mesh
  
  ! Meshes.
  use OpenFOAMMesh_class,    only : OpenFOAMMesh
  
  implicit none
  private
  
  ! ** ADD NAME OF NEW MESH TO THE LIST **!
  ! List that contains acceptable types of meshes
  ! NOTE: It is necessary to adjust trailing blanks so all entries have the same length
  character(nameLen), dimension(*), parameter :: AVAILABLE_MESHES = ['OpenFOAMMesh']
  ! Public interface.
  public :: new_mesh_ptr

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

end module meshFactory_func