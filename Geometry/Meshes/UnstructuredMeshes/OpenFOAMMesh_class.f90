module OpenFOAMMesh_class

  use numPrecision
  use dictionary_class,       only : dictionary
  use OpenFOAMFunctions,      only : importMesh
  use unstructuredMesh_inter, only : unstructuredMesh

  implicit none
  private

  type, public, extends(unstructuredMesh) :: OpenFOAMMesh
    private
  contains
    procedure, non_overridable            :: init
  end type OpenFOAMMesh

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
  subroutine init(self, folderPath, dict)
    class(OpenFOAMMesh), intent(inout) :: self
    character(*), intent(in)           :: folderPath
    class(dictionary), intent(in)      :: dict
    
    ! Import OpenFOAM mesh and initialise kd-tree.
    call importMesh(self, folderPath, dict)
    call self % tree % init(self % vertices % getAllCoordinates(), .true.)

  end subroutine init

end module OpenFOAMMesh_class