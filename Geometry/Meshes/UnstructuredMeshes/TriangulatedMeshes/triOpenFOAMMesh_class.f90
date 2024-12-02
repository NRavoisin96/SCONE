module triOpenFOAMMesh_class

  use numPrecision
  use dictionary_class,          only : dictionary
  use OpenFOAMFunctions,         only : importMesh
  use triUnstructuredMesh_inter, only : triUnstructuredMesh

  implicit none
  private

  type, public, extends(triUnstructuredMesh) :: triOpenFOAMMesh
    private
  contains
    procedure, non_overridable               :: init
  end type triOpenFOAMMesh

contains

  !! Subroutine 'init'
  !!
  !! Basic description:
  !!   Imports an OpenFOAM mesh from the path of the folder containing the mesh files and sets the 
  !!   mesh Id from the dictionary.
  !!
  !! Detailed description:
  !!   'init' first sets the mesh Id from the supplied dictionary. It then checks the existence and 
  !!   consistency of files located in the appropriate mesh folder and allocates memory to the 
  !!   'vertices', 'faces' and 'elements' structures of the 'mesh' structure. It then imports data 
  !!   from the 'points' file into the 'vertices' structures and builds a kd-tree for the mesh from 
  !!   the various vertices' coordinates. The subroutine then proceeds to import data from the 
  !!   'faces', 'owner' and 'neighbour' files and stores it into the appropriate 'faces' and 
  !!   'elements' structure. 'init' also computes the area, centroid and normal vector for each 
  !!   face, as well as the centroid and volume of each element. From the 'faces', 'owner' and 
  !!   'neighbour' files, mesh connectivity information is also assigned to the various structures. 
  !!   The subroutine then imports data from the 'cellZones' file (if it exists) into the 
  !!   'cellZones' structures.
  !!
  !! Arguments:
  !!   folderPath [in] -> Path of the folder containing the files for the mesh geometry.
  !!   dict [in]       -> Input dictionary.
  !!
  !! Errors:
  !!   - fatalError if id < 1.
  !!   - fatalError if the mesh contains concave elements.
  !!
  subroutine init(self, folderPath, dict)
    class(triOpenFOAMMesh), intent(inout)        :: self
    character(*), intent(in)                     :: folderPath
    class(dictionary), intent(in)                :: dict
    integer(shortInt)                            :: freeVertexIdx
    
    ! Import OpenFOAM mesh and split non-tetrahedral mesh elements into tetrahedra.
    call importMesh(self, folderPath, dict)
    call self % split(freeVertexIdx)
    
    ! Collapse the 'vertices' structure to save memory and initialise kd-tree for the mesh.
    call self % vertices % collapse(freeVertexIdx)
    call self % tree % init(self % vertices % getAllCoordinates(), .true.)

  end subroutine init

end module triOpenFOAMMesh_class