module triOpenFOAMMesh_class

  use numPrecision
  use coord_class,               only : coord
  use dictionary_class,          only : dictionary
  use OpenFOAMFunctions,         only : importMesh
  use triUnstructuredMesh_inter, only : triUnstructuredMesh, &
                                        distanceToBoundaryFace_super => distanceToBoundaryFace, &
                                        distanceToNextFace_super => distanceToNextFace, &
                                        findElementAndParentIdxs_super => findElementAndParentIdxs, &
                                        kill_super => kill

  implicit none
  private

  type, public, extends(triUnstructuredMesh) :: triOpenFOAMMesh
    private
  contains
    ! Superclass procedures.
    procedure                                :: distanceToBoundaryFace
    procedure                                :: distanceToNextFace
    procedure                                :: findElementAndParentIdxs
    procedure                                :: init
    ! Local procedures.
    procedure                                :: kill
  end type triOpenFOAMMesh

contains

  !! Subroutine 'distanceToBoundaryFace'
  !!
  !! Basic description:
  !!   Returns the distance to the mesh boundary face intersected by a particle's path. Also returns the index
  !!   of the parent element containing the intersected boundary face.
  !!
  !! See unstructuredMesh_inter for details.
  !!
  elemental subroutine distanceToBoundaryFace(self, d, coords, parentIdx)
    class(triOpenFOAMMesh), intent(in) :: self
    real(defReal), intent(out)         :: d
    type(coord), intent(inout)         :: coords
    integer(shortInt), intent(out)     :: parentIdx

    call distanceToBoundaryFace_super(self, d, coords, parentIdx)

  end subroutine distanceToBoundaryFace

  !! Subroutine 'distanceToNextFace'
  !!
  !! Basic description:
  !!   Returns the distance to the next face intersected by the particle's path.
  !!
  !! See unstructuredMesh_inter for details.
  !!
  elemental subroutine distanceToNextFace(self, d, coords)
    class(triOpenFOAMMesh), intent(in) :: self
    real(defReal), intent(out)         :: d
    type(coord), intent(inout)         :: coords

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
    class(triOpenFOAMMesh), intent(in)      :: self
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
    class(triOpenFOAMMesh), intent(inout) :: self
    character(*), intent(in)              :: folderPath
    class(dictionary), intent(in)         :: dict
    
    ! Import OpenFOAM mesh, split non-tetrahedral mesh elements into tetrahedra and 
    ! initialise kd-tree for the mesh.
    call importMesh(self, folderPath, dict)
    call self % split()
    call self % tree % init(self % getAllVertexCoordinates(), .true.)

  end subroutine init

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an unitialised state.
  !!
  elemental subroutine kill(self)
    class(triOpenFOAMMesh), intent(inout) :: self

    ! Call triUnstructuredMesh procedure.
    call kill_super(self)

  end subroutine kill

end module triOpenFOAMMesh_class