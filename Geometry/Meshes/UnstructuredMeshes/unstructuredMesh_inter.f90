module unstructuredMesh_inter

  use numPrecision
  use universalVariables
  use genericProcedures,   only : findDifferent, numToChar
  use coord_class,         only : coord
  use dictionary_class,    only : dictionary
  use mesh_inter,          only : mesh, kill_super => kill
  use element_class,       only : element
  use elementShelf_class,  only : elementShelf
  use face_class,          only : face
  use faceShelf_class,     only : faceShelf
  use vertexShelf_class,   only : vertexShelf
  use kdTree_class,        only : kdTree

  implicit none
  private

  ! Extendable procedures.
  public :: kill, distanceToNextFace, distanceToBoundaryFace, findElementAndParentIdxs
  
  !! Abstract interface to group all unstructured meshes. An unstructured mesh uses a vertex -> face 
  !! -> element representation of space. Each element is composed by a set of faces which are themselves 
  !! composed by a number of vertices. Elements can be grouped together into zones. This is useful to 
  !! assign material filling to mesh elements. Local ids are assigned in the order of the cell zone 
  !! definition.
  !!
  !! Public members:
  !!   cellZones                -> Shelf that stores cell zones.
  !!   edges                    -> Shelf that stores edges.
  !!   elements                 -> Shelf that stores elements.
  !!   faces                    -> Shelf that stores faces.
  !!   vertices                 -> Shelf that stores vertices.
  !!   nVertices                -> Number of vertices in the mesh.
  !!   nFaces                   -> Number of faces in the mesh.
  !!   nEdges                   -> Number of edges in the mesh.
  !!   nElements                -> Number of elements in the mesh.
  !!   nInternalFaces           -> Number of internal faces in the mesh.
  !!   tree                     -> kd-tree used for nearest-neighbour searches and entry checks.
  !!
  !! Interface:
  !!   kill                     -> Returns to an unitialised state.
  !!   printComposition         -> Displays mesh composition to the user.
  !!   distanceToBoundaryFace   -> Checks if a particle enters the mesh and returns distance to entry 
  !!                               intersection.
  !!   distanceToNextFace       -> Returns the distance to the next mesh face.
  !!   findElementAndParentIdxs -> Returns the index of the mesh element occupied by a particle. Also
  !!                               returns the index of the parent mesh element containing the occupied
  !!                               element.
  !!
  type, public, abstract, extends(mesh) :: unstructuredMesh
    private
    integer(shortInt), public           :: nVertices = 0, nFaces = 0, nEdges = 0, &
                                           nElements = 0, nInternalFaces = 0
!    type(edgeShelf), public            :: edges
    type(elementShelf), public          :: elements
    type(faceShelf), public             :: faces
    type(vertexShelf), public           :: vertices
    type(kdTree), public                :: tree
  contains
    ! Build procedures.
    procedure                           :: kill
    procedure, non_overridable          :: printComposition
    ! Runtime procedures.
    procedure                           :: distanceToBoundaryFace
    procedure                           :: distanceToNextFace
    procedure                           :: findElementAndParentIdxs
  end type unstructuredMesh

contains

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an unitialised state.
  !!
  elemental subroutine kill(self)
    class(unstructuredMesh), intent(inout) :: self

    ! Superclass.
    call kill_super(self)
    
    ! Local.
    self % nVertices = 0
    self % nFaces = 0
    self % nInternalFaces = 0
    self % nElements = 0
    self % nEdges = 0
    call self % vertices % kill()
    call self % faces % kill()
    call self % elements % kill()
    call self % tree % kill()

  end subroutine kill

  !! Subroutine 'printComposition'
  !!
  !! Basic description:
  !!   Prints the initial polyhedral composition of the mesh.
  !!
  !! Arguments:
  !!   nTetrahedra [out] -> Number of tetrahedra in the mesh.
  !!
  subroutine printComposition(self, nTetrahedra)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt), intent(out)      :: nTetrahedra
    integer(shortInt)                   :: nFaces, nPentahedra, nHexahedra, nOthers, i
    
    ! Initialise the numbers of various polyhedra to zero.
    nTetrahedra = 0
    nPentahedra = 0
    nHexahedra = 0
    nOthers = 0
    
    ! Loop over all elements in the mesh.
    do i = 1, self % nElements
      ! Retrieve the number of faces in the current element and increment specific polyhedra
      ! accordingly.
      nFaces = size(self % elements % shelf(i) % getFaces())
      select case (nFaces)
        case (4)
          nTetrahedra = nTetrahedra + 1
        case (5)
          nPentahedra = nPentahedra + 1
        case (6)
          nHexahedra = nHexahedra + 1
        case default
          nOthers = nOthers + 1

      end select

    end do
    
    ! Print to screen.
    print *, 'Displaying unstructured mesh composition:'
    print *, '  Number of tetrahedra     : '//numToChar(nTetrahedra)//'.'
    print *, '  Number of pentahedra     : '//numToChar(nPentahedra)//'.'
    print *, '  Number of hexahedra      : '//numToChar(nHexahedra)//'.'
    print *, '  Number of other polyhedra: '//numToChar(nOthers)//'.'

  end subroutine printComposition

  !! Subroutine 'distanceToBoundaryFace'
  !!
  !! Basic description:
  !!   Returns the distance to the mesh boundary face intersected by a particle's path. Also returns the index
  !!   of the parent element containing the intersected boundary face.
  !!
  !! See mesh_inter for details.
  !! TODO: finish this.
  pure subroutine distanceToBoundaryFace(self, d, coords, parentIdx)
    class(unstructuredMesh), intent(in) :: self
    real(defReal), intent(out)          :: d
    type(coord), intent(inout)          :: coords
    integer(shortInt), intent(out)      :: parentIdx
    
    ! TODO: find intersected face from tree.

  end subroutine distanceToBoundaryFace

  !! Subroutine 'distanceToNextFace'
  !!
  !! Basic description:
  !!   Returns the distance to the next face intersected by the particle's path. Returns INF if the particle
  !!   does not intersect any face (i.e., if its path is entirely contained in the element the particle 
  !!   currently is). Algorithm adapted from Macpherson, et al. (2009). DOI: 10.1002/cnm.1128.
  !!
  !! See mesh_inter for details.
  !! TODO: finish this.
  !!
  pure subroutine distanceToNextFace(self, d, coords)
    class(unstructuredMesh), intent(in)          :: self
    real(defReal), intent(out)                   :: d
    type(coord), intent(inout)                   :: coords
    real(defReal), dimension(3)                  :: r, rEnd
    type(element)                                :: currentElement
    type(face)                                   :: intersectedFace
    
    ! Initialise d = INF, retrieve the tetrahedron currently occupied by the particle and compute
    ! potential triangle intersections.
    d = INF
    currentElement = self % elements % shelf(coords % elementIdx)
    rEnd = coords % rEnd

  end subroutine distanceToNextFace

  !! Subroutine 'findElementAndParentIdxs'
  !!
  !! Basic description:
  !!   Returns the index of the mesh element occupied by a particle. Also returns the index of the parent mesh
  !!   element containing the occupied element.
  !!
  !! See mesh_inter for details.
  !!
  pure subroutine findElementAndParentIdxs(self, r, u, elementIdx, parentIdx)
    class(unstructuredMesh), intent(in)          :: self
    real(defReal), dimension(3), intent(in)      :: r, u
    integer(shortInt), intent(out)               :: elementIdx, parentIdx
    real(defReal), dimension(6)                  :: boundingBox
    integer(shortInt), dimension(:), allocatable :: potentialElements, faceToElements, zeroDotProductFaceIdxs
    type(element)                                :: currentElement
    integer(shortInt)                            :: i, nearestVertexIdx, failedFaceIdx
    type(face)                                   :: failedFace
    
    ! Initialise elementIdx = 0 and parentIdx = 0. Retrieve the mesh's bounding box. If the particle is outside the
    ! bounding box we can return early.
    elementIdx = 0
    parentIdx = 0
    boundingBox = self % getBoundingBox()
    do i = 1, 3
      if (r(i) <= boundingBox(i) .or. boundingBox(i + 3) <= r(i)) return

    end do

    ! If the point is inside the bounding box then we need to determine if the point is inside a
    ! pseudoCell. First find the vertex which is nearest to the coordinates.
    nearestVertexIdx = self % tree % findNearestVertex(r)
    
    ! Retrieve potential elements occupied by the particle from mesh connectivity information.
    potentialElements = self % vertices % shelf(nearestVertexIdx) % getVertexToElements()
    
    ! Initialise the search to the first element in potentialElements.
    currentElement = self % elements % shelf(potentialElements(1))
    searchLoop: do
      ! Check if the particle is inside the current element.
      call currentElement % testForInclusion(self % faces, r, failedFaceIdx, zeroDotProductFaceIdxs)
      
      ! If the current element is not occupied by the particle, retrieve the element sharing the face
      ! for which the search failed and update the search.
      if (failedFaceIdx > 0) then
        failedFace = self % faces % shelf(failedFaceIdx)
        
        ! If the failed face is a boundary face the particle is outside the mesh and we can return.
        if (failedFace % isBoundary()) return
        
        ! Retrieve the elements sharing the face from mesh connectivity, update element to be searched
        ! and cycle.
        faceToElements = failedFace % getFaceToElements()
        elementIdx = findDifferent(faceToElements, currentElement % getIdx())
        currentElement = self % elements % shelf(elementIdx)
        cycle searchLoop

      end if
      
      ! If there are faces on which the particle lies we need to employ some more specific
      ! procedures to correctly determine which element is actually occupied.
      if (allocated(zeroDotProductFaceIdxs)) then
        select case (size(zeroDotProductFaceIdxs))
          ! Particle is on a element's face.
          case(1)
          ! Particle is on a element's edge.
          case(2)
          ! Particle is on a element's vertex.
          case(3)

        end select
        exit
      
      end if
      ! If reached this point the particle is in the current element. Update elementIdx and parentIdx and return.
      elementIdx = currentElement % getIdx()
      parentIdx = elementIdx
      return

    end do searchLoop

  end subroutine findElementAndParentIdxs

end module unstructuredMesh_inter