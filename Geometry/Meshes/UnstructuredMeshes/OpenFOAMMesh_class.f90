module OpenFOAMMesh_class

  use numPrecision
  use coord_class,            only : coord
  use dictionary_class,       only : dictionary
  use OpenFOAMFunctions,      only : importMesh
  use unstructuredMesh_inter, only : unstructuredMesh, &
                                     distanceToBoundaryFace_super => distanceToBoundaryFace, &
                                     distanceToNextFace_super => distanceToNextFace, &
                                     findElementAndParentIdxs_super => findElementAndParentIdxs

  implicit none
  private

  type, public, extends(unstructuredMesh) :: OpenFOAMMesh
    private
  contains
    ! Superclass procedures.
    procedure                             :: distanceToBoundaryFace
    procedure                             :: distanceToNextFace
    procedure                             :: findElementAndParentIdxs
    procedure                             :: init
  end type OpenFOAMMesh

contains

  !! Subroutine 'distanceToBoundaryFace'
  !!
  !! Basic description:
  !!   Returns the distance to the mesh boundary face intersected by a particle's path. Also returns the index
  !!   of the parent element containing the intersected boundary face.
  !!
  !! See unstructuredMesh_inter for details.
  !!
  pure subroutine distanceToBoundaryFace(self, d, coords, parentIdx)
    class(OpenFOAMMesh), intent(in) :: self
    real(defReal), intent(out)      :: d
    type(coord), intent(inout)      :: coords
    integer(shortInt), intent(out)  :: parentIdx

    call distanceToBoundaryFace_super(self, d, coords, parentIdx)

  end subroutine distanceToBoundaryFace

  !! Subroutine 'distanceToNextFace'
  !!
  !! Basic description:
  !!   Returns the distance to the next face intersected by the particle's path.
  !!
  !! See unstructuredMesh_inter for details.
  !!
  pure subroutine distanceToNextFace(self, d, coords)
    class(OpenFOAMMesh), intent(in) :: self
    real(defReal), intent(out)      :: d
    type(coord), intent(inout)      :: coords

    call distanceToNextFace_super(self, d, coords)

  end subroutine distanceToNextFace

  !! Subroutine 'findElementAndParentIdxs'
  !!
  !! Basic description:
  !!   Returns the index of the mesh element occupied by a particle. Also returns the index of the parent mesh
  !!   element containing the occupied element.
  !!
  !! See unstructuredMesh_inter for details.
  !!
  pure subroutine findElementAndParentIdxs(self, r, u, elementIdx, parentIdx)
    class(OpenFOAMMesh), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r, u
    integer(shortInt), intent(out)          :: elementIdx, parentIdx

    call findElementAndParentIdxs_super(self, r, u, elementIdx, parentIdx)

  end subroutine findElementAndParentIdxs

  !! Subroutine 'init'
  !!
  !! Basic description:
  !!   Imports an OpenFOAM mesh from the path of the folder containing the mesh files.
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