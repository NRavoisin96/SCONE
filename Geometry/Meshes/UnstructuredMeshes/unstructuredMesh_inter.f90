module unstructuredMesh_inter

  use numPrecision
  use universalVariables
  use genericProcedures,   only : findDifferent, numToChar
  use coord_class,         only : coord
  use dictionary_class,    only : dictionary
  use mesh_inter,          only : mesh, kill_super => kill
  use cellZoneShelf_class, only : cellZoneShelf
  use element_class,       only : element
  use elementShelf_class,  only : elementShelf
  use face_class,          only : face
  use faceShelf_class,     only : faceShelf
  use vertexShelf_class,   only : vertexShelf
  use kdTree_class,        only : kdTree

  implicit none
  private

  ! Extendable methods.
  public :: kill
  
  !! Abstract interface to group all unstructured meshes. An unstructured mesh uses a vertex -> face 
  !! -> element representation of space. Each element is composed by a set of faces which are themselves 
  !! composed by a number of vertices. Elements can be grouped together into zones. This is useful to 
  !! assign material filling to mesh elements. Local ids are assigned in the order of the cell zone 
  !! definition.
  !!
  !! Public members:
  !!   cellZones              -> Shelf that stores cell zones.
  !!   edges                  -> Shelf that stores edges.
  !!   elements               -> Shelf that stores elements.
  !!   faces                  -> Shelf that stores faces.
  !!   vertices               -> Shelf that stores vertices.
  !!
  !! Private members:
  !!   nVertices              -> Number of vertices in the mesh.
  !!   nFaces                 -> Number of faces in the mesh.
  !!   nEdges                 -> Number of edges in the mesh.
  !!   nElements              -> Number of elements in the mesh.
  !!   nInternalFaces         -> Number of internal faces in the mesh.
  !!   cellZonesFile          -> .true. if a file is present to detail the different cell zones.
  !!   tree                   -> kd-tree used for nearest-neighbour searches and entry checks.
  !!
  !! Interface:
  !!   checkFiles             -> Checks the existence of required files to import the mesh.
  !!   getMeshInfo            -> Retrieves preleminary information about mesh composition.
  !!   kill                   -> Returns to an unitialised state.
  !!   printComposition       -> Displays mesh composition to the user.
  !!   checkForEntry          -> Checks if a particle enters the mesh and returns distance to entry 
  !!                             intersection.
  !!   distanceToNextFace     -> Returns the distance to the next mesh face.
  !!   distance               -> Returns the distance travelled by a particle within the mesh.
  !!   findElement            -> Returns the index of the mesh element occupied by a particle.
  !!   getVerticesNumber      -> Returns the number of vertices in the mesh.
  !!   getFacesNumber         -> Returns the number of faces in the mesh.
  !!   getElementsNumber      -> Returns the number of elements in the mesh.
  !!   getInternalFacesNumber -> Returns the number of internal faces in the mesh.
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
    procedure                           :: checkForEntry
    procedure                           :: distanceToNextFace
    procedure                           :: distance
    procedure                           :: findElement
    procedure, non_overridable          :: getVerticesNumber
    procedure, non_overridable          :: getFacesNumber
    procedure, non_overridable          :: getElementsNumber
    procedure, non_overridable          :: getInternalFacesNumber
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
    self % nElementZones = 0
    self % cellZonesFile = .false.
    if (allocated(self % vertices % shelf)) call self % vertices % kill()
    if (allocated(self % faces % shelf)) call self % faces % kill()
    if (allocated(self % elements % shelf)) call self % elements % kill()
    if (allocated(self % cellZones % shelf)) call self % cellZones % kill()
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

  !! Subroutine 'checkForEntry'
  !!
  !! Basic description:
  !!   Checks whether a given particle enters from the CSG cell into the mesh. Returns the distance
  !!   to the face of intersection in case the particle enters the mesh. Returns INF otherwise.
  !!
  !! Arguments:
  !!   d [out]        -> Distance to the intersected triangle.
  !!   coords [inout] -> Particle's coordinates.
  !!
  !! TODO: finish this.
  pure subroutine checkForEntry(self, d, coords)
    class(unstructuredMesh), intent(in) :: self
    real(defReal), intent(out)          :: d
    type(coord), intent(inout)          :: coords
    real(defReal)                       :: rComponent, rEndComponent
    real(defReal), dimension(6)         :: boundingBox
    integer(shortInt)                   :: i
    
    ! Initialise d = INF and check that the particle's path intersects the mesh's bounding box.
    d = INF
    boundingBox = self % getBoundingBox()
    do i = 1, 3
      ! If both the particle's initial and end locations are on the same side of the plane of the
      ! bounding box return early.
      rComponent = coords % r(i)
      rEndComponent = coords % rEnd(i)
      if ((rComponent < boundingBox(i) .and. rEndComponent < boundingBox(i)) .or. &
          (rComponent > boundingBox(i + 3) .and. rEndComponent > boundingBox(i + 3))) return

    end do
    
    ! TODO: find intersected face from tree.

  end subroutine checkForEntry

  !! Subroutine 'distanceToNextFace'
  !!
  !! Basic description:
  !!   Returns the distance to the next triangle intersected by the particle's path. Returns
  !!   INF if the particle does not intersect any triangle (i.e., if its path is entirely
  !!   contained in the tetrahedron the particle currently is). Algorithm adapted from Macpherson, 
  !!   et al. (2009). DOI: 10.1002/cnm.1128.
  !!
  !! Arguments:
  !!   d [out]        -> Distance to the next intersected triangle.
  !!   coords [inout] -> Particle's coordinates.
  !!
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

  !! Subroutine 'distance'
  !!
  !! Basic description:
  !!   Computes the distance of intersection between a particle and the mesh. Also updates the element of the mesh
  !!   occupied by the particle and the associated local id.
  !!
  !! Arguments:
  !!   d [out]        -> Distance to next intersection.
  !!   coords [inout] -> Particle's coordinates.
  !!   isInside [out] -> .true. if particle intersects a face of the mesh.
  !!
  pure subroutine distance(self, d, coords, isInside)
    class(unstructuredMesh), intent(in) :: self
    real(defReal), intent(out)          :: d
    type(coord), intent(inout)          :: coords
    logical(defBool), intent(out)       :: isInside

    ! Initialise isInside and compute the particle's end position.
    isInside = .true.
    
    ! If particle is already inside a tetrahedron, simply compute the distance to the next mesh triangle and return.
    if (coords % elementIdx > 0) then
      call self % distanceToNextFace(d, coords)
      return

    end if

    ! If not, we need to check if the particle enters the mesh.
    call self % checkForEntry(d, coords)

    ! If the particle enters the mesh retrieve the localId based on the element associated with the entry tetrahedron and return.
    if (coords % elementIdx > 0) then
      coords % localId = self % cellZones % findCellZone(coords % elementIdx)
      return

    end if
      
    ! If reached here, the particle does not enter the mesh and CSG tracking resumes.
    isInside = .false.

  end subroutine distance

  !! Subroutine 'findElement'
  !!
  !! Basic description:
  !!   Returns the index of the mesh element occupied by a particle as well as the associated local id.
  !!
  !! Arguments:
  !!   r [in]           -> Particle's location.
  !!   u [in]           -> Particle's direction.
  !!   elementIdx [out] -> Index of the mesh element occupied by the particle. 0 if particle is outside the mesh.
  !!   localId [out]    -> Local id of the element zone containing the element occupied by the particle.
  !!
  pure subroutine findElement(self, r, u, elementIdx, localId)
    class(unstructuredMesh), intent(in)          :: self
    real(defReal), dimension(3), intent(in)      :: r, u
    integer(shortInt), intent(out)               :: elementIdx, localId
    real(defReal), dimension(6)                  :: boundingBox
    integer(shortInt), dimension(:), allocatable :: potentialElements, faceToElements, zeroDotProductFaceIdxs
    type(element)                                :: currentElement
    integer(shortInt)                            :: i, nearestVertexIdx, failedFaceIdx
    type(face)                                   :: failedFace
    
    ! Initialise localId = 1 (corresponds to the particle being in the CSG cell) and elementIdx = 0.
    localId = 1
    elementIdx = 0

    ! Retrieve the mesh's bounding box. If the particle is outside the bounding box we can return early.
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
      ! If reached this point the particle is in the current element. Retrieve its index and exit.
      elementIdx = currentElement % getIdx()
      exit

    end do searchLoop

    ! If elementIdx = 0 the particle is in the CSG cell. Return early in this case.
    if (elementIdx == 0) return

    ! Update localId based on the element associated with the element.
    localId = self % cellZones % findCellZone(elementIdx)

  end subroutine findElement

  !! Function 'getVeticesNumber'
  !!
  !! Basic description:
  !!   Returns the number of vertices in the mesh.
  !!
  elemental function getVerticesNumber(self) result(nVertices)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt)                   :: nVertices

    nVertices = self % nVertices

  end function getVerticesNumber

  !! Function 'getFacesNumber'
  !!
  !! Basic description:
  !!   Returns the number of faces in the mesh.
  !!
  elemental function getFacesNumber(self) result(nFaces)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt)                   :: nFaces

    nFaces = self % nFaces

  end function getFacesNumber

  !! Function 'getElementsNumber'
  !!
  !! Basic description:
  !!   Returns the number of elements in the mesh.
  !!
  elemental function getElementsNumber(self) result(nElements)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt)                   :: nElements

    nElements = self % nElements

  end function getElementsNumber

  !! Function 'getInternalFacesNumber'
  !!
  !! Basic description:
  !!   Returns the number of internal faces in the mesh.
  !!
  elemental function getInternalFacesNumber(self) result(nInternalFaces)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt)                   :: nInternalFaces

    nInternalFaces = self % nInternalFaces

  end function getInternalFacesNumber

end module unstructuredMesh_inter