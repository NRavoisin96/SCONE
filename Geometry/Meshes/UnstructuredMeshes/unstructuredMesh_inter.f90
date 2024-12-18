module unstructuredMesh_inter

  use numPrecision
  use universalVariables
  use genericProcedures,   only : findDifferent, numToChar
  use coord_class,         only : coord
  use dictionary_class,    only : dictionary
  use mesh_inter,          only : mesh, kill_super => kill
  use edgeShelf_class,     only : edgeShelf
  use elementShelf_class,  only : elementShelf
  use faceShelf_class,     only : faceShelf
  use vertexShelf_class,   only : vertexShelf
  use kdTree_class,        only : kdTree

  implicit none
  private

  ! Extendable procedures.
  public :: distanceToBoundaryFace, distanceToNextFace, findElementAndParentIdxs, kill
  
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
    type(edgeShelf), public             :: edges
    type(elementShelf), public          :: elements
    type(faceShelf), public             :: faces
    type(vertexShelf), public           :: vertices
    type(kdTree), public                :: tree
  contains
    ! Build procedures.
    procedure                           :: kill
    procedure, non_overridable          :: printComposition
    procedure                           :: setEdgeShelf
    procedure                           :: setElementShelf
    procedure                           :: setFaceShelf
    procedure                           :: setVertexShelf
    ! Runtime procedures.
    procedure                           :: distanceToBoundaryFace
    procedure                           :: distanceToNextFace
    procedure                           :: findElementAndParentIdxs
    procedure                           :: findElementFromDirection
    procedure                           :: getAllVertexCoordinates
  end type unstructuredMesh

contains

  !! Subroutine 'distanceToBoundaryFace'
  !!
  !! Basic description:
  !!   Returns the distance to the mesh boundary face intersected by a particle's path. Also returns the index
  !!   of the parent element containing the intersected boundary face.
  !!
  !! See mesh_inter for details.
  !!
  elemental subroutine distanceToBoundaryFace(self, d, coords, parentIdx)
    class(unstructuredMesh), intent(in)          :: self
    real(defReal), intent(out)                   :: d
    type(coord), intent(inout)                   :: coords
    integer(shortInt), intent(out)               :: parentIdx
    integer(shortInt)                            :: edgeIdx, vertexIdx
    integer(shortInt), dimension(:), allocatable :: testFaceIdxs, elementIdxs
    
    ! Initialise parentIdx = 0, edgeIdx = 0 and vertexIdx = 0 then search the tree for the intersected boundary face.
    parentIdx = 0
    edgeIdx = 0
    vertexIdx = 0
    call self % tree % findIntersectedFace(self % vertices, self % faces, d, coords, edgeIdx, vertexIdx)

    ! If edgeIdx > 0 or vertexIdx > 0 then the particle intersects a boundary face through and edge or a vertex. In this
    ! case assign element occupation from particle direction.
    if (edgeIdx > 0 .or. vertexIdx > 0) then
      if (edgeIdx > 0) then
        testFaceIdxs = self % edges % getEdgeFaceIdxs(edgeIdx)
        elementIdxs = self % edges % getEdgeElementIdxs(edgeIdx)

      elseif (vertexIdx > 0) then
        testFaceIdxs = self % vertices % getVertexFaceIdxs(vertexIdx)
        elementIdxs = self % vertices % getVertexElementIdxs(vertexIdx)

      end if
      coords % elementIdx = self % findElementFromDirection(coords % dir, elementIdxs, testFaceIdxs)
      
      ! If the particle's elementIdx is 0 (e.g., particle enters the mesh via a boundary vertex but still points outside)
      ! then update d = INF and return.
      if (coords % elementIdx == 0) then
        d = INF
        return

      end if

    end if

    ! Update parentIdx only if particle enters the mesh.
    if (coords % elementIdx > 0) parentIdx = coords % elementIdx

  end subroutine distanceToBoundaryFace

  !! Subroutine 'distanceToNextFace'
  !!
  !! Basic description:
  !!   Returns the distance to the next face intersected by the particle's path. Returns INF if the particle
  !!   does not intersect any face (i.e., if its path is entirely contained in the element the particle 
  !!   currently is). Algorithm adapted from Macpherson, et al. (2009). DOI: 10.1002/cnm.1128.
  !!
  !! See mesh_inter for details.
  !!
  elemental subroutine distanceToNextFace(self, d, coords)
    class(unstructuredMesh), intent(in)          :: self
    real(defReal), intent(out)                   :: d
    type(coord), intent(inout)                   :: coords
    real(defReal), dimension(3)                  :: r, rEnd
    integer(shortInt), dimension(:), allocatable :: potentialFaces, faceToElements
    integer(shortInt)                            :: elementIdx, intersectedFaceIdx
    real(defReal)                                :: lambda
    
    ! Initialise d = INF, retrieve the element currently occupied by the particle and compute potential 
    ! face intersections.
    d = INF
    elementIdx = coords % elementIdx
    rEnd = coords % rEnd
    potentialFaces = self % elements % computePotentialFaceIdxs(elementIdx, rEnd, self % faces)

    ! If no potential intersections are detected return early.
    if (size(potentialFaces) == 0) return
    
    ! If reached here, compute which face is actually intersected and update d.
    r = coords % r
    call self % elements % computeFaceIntersection(elementIdx, r, rEnd, potentialFaces, self % faces, &
                                                   intersectedFaceIdx, lambda)
    d = norm2(min(ONE, max(ZERO, lambda)) * (rEnd - r))
    
    ! If the intersected face is a boundary face then the particle is leaving the mesh.
    if (self % faces % getFaceIsBoundary(intersectedFaceIdx)) then
      coords % elementIdx = 0
      coords % localId = 1
      return

    end if

    ! Else, retrieve the elements sharing the intersected face from mesh connectivity then
    ! update elementIdx and localId.
    faceToElements = self % faces % getFaceElementIdxs(intersectedFaceIdx)
    coords % elementIdx = findDifferent(faceToElements, elementIdx)
    coords % localId = self % findElementZoneIdx(coords % elementIdx)

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
    integer(shortInt), dimension(:), allocatable :: potentialElements, faceToElements, zeroDotProductFaceIdxs, &
                                                    testFaceIdxs
    integer(shortInt)                            :: i, nearestVertexIdx, failedFaceIdx, commonEdgeIdx, commonVertexIdx
    
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
    potentialElements = self % vertices % getVertexElementIdxs(nearestVertexIdx)
    
    ! Initialise the search to the first element in potentialElements.
    elementIdx = potentialElements(1)
    searchLoop: do
      ! Check if the particle is inside the current element.
      call self % elements % testForInclusion(elementIdx, r, self % faces, failedFaceIdx, zeroDotProductFaceIdxs)
      
      ! If the current element is not occupied by the particle, retrieve the element sharing the face
      ! for which the search failed and update the search.
      if (failedFaceIdx > 0) then
        ! If the failed face is a boundary face the particle is outside the mesh and we can return.
        if (self % faces % getFaceIsBoundary(failedFaceIdx)) return
        
        ! Retrieve the elements sharing the face from mesh connectivity, update element to be searched
        ! and cycle.
        faceToElements = self % faces % getFaceElementIdxs(failedFaceIdx)
        elementIdx = findDifferent(faceToElements, elementIdx)
        cycle searchLoop

      end if
      
      ! If there are faces on which the particle lies we need to employ some more specific
      ! procedures to correctly determine which element is actually occupied.
      if (allocated(zeroDotProductFaceIdxs)) then
        select case (size(zeroDotProductFaceIdxs))
          ! Particle is on a element's face.
          case(1)
            potentialElements = self % faces % getFaceElementIdxs(zeroDotProductFaceIdxs(1))
            testFaceIdxs = zeroDotProductFaceIdxs
          ! Particle is on a element's edge.
          case(2)
            ! Find the common edge and retrieve elements sharing this edge.
            commonEdgeIdx = self % faces % findCommonEdgeIdx(zeroDotProductFaceIdxs(1), zeroDotProductFaceIdxs(2))
            potentialElements = self % edges % getEdgeElementIdxs(commonEdgeIdx)
            testFaceIdxs = self % edges % getEdgeFaceIdxs(commonEdgeIdx)
          ! Particle is on a element's vertex.
          case default
            ! Find the common vertex and retrieve all elements sharing this vertex.
            commonVertexIdx = self % faces % findCommonVertexIdx(zeroDotProductFaceIdxs)
            potentialElements = self % vertices % getVertexElementIdxs(commonVertexIdx)
            testFaceIdxs = self % vertices % getVertexFaceIdxs(commonVertexIdx)

        end select
        elementIdx = self % findElementFromDirection(u, potentialElements, testFaceIdxs)
        parentIdx = elementIdx
        return
      
      end if
      ! If reached this point the particle is in the current element. Update parentIdx and return.
      parentIdx = elementIdx
      return

    end do searchLoop

  end subroutine findElementAndParentIdxs

  !! Function 'findElementFromDirection'
  !!
  !! Basic description:
  !!   Returns the index of the element occupied by a particle in case the particle lies on one or more
  !!   element face(s).
  !!
  !! Detailed description:
  !!   When the particle lies on one or more element face(s), it is technically in more than one element
  !!   at once. In this case element occupation can be assigned by using the particle's direction vector.
  !!   The idea is to take all the faces on which the particle lies and perform the dot product of the
  !!   particle's direction with each of the face's (signed) normal vector. The element actually occupied
  !!   by the partice is then that for which the number of negative dot products is the greatest.
  !!
  !! Arguments:
  !!   u [in]            -> Particle's direction.
  !!   elementIdxs [in]  -> Indices of the elements to be tested.
  !!   testFaceIdxs [in] -> Indices of the faces to be tested.
  !!
  !! Result:
  !!   elementIdx        -> Index of the element occupied by the particle.
  !!
  pure function findElementFromDirection(self, u, elementIdxs, testFaceIdxs) result(elementIdx)
    class(unstructuredMesh), intent(in)          :: self
    real(defReal), dimension(3), intent(in)      :: u
    integer(shortInt), dimension(:), intent(in)  :: elementIdxs, testFaceIdxs
    integer(shortInt)                            :: elementIdx, i, j, count, maxCount, maxIdx
    integer(shortInt), dimension(:), allocatable :: faceIdxs, absFaceIdxs
    real(defReal), dimension(3)                  :: normal

    ! Initialise elementIdx = 0 and maxCount = 0 then loop over all elements.
    elementIdx = 0
    maxCount = 0
    do i = 1, size(elementIdxs)
      ! Retrieve the faces making the current element and convert them to absolute indices.
      faceIdxs = self % elements % getElementFaceIdxs(elementIdxs(i))
      absFaceIdxs = abs(faceIdxs)
      
      ! Initialise count = 0 then loop through all faces in the current element.
      count = 0
      do j = 1, size(faceIdxs)

        ! Cycle to the next face if it is not a face on which the particle is.
        if (.not. any(testFaceIdxs == absFaceIdxs(j))) cycle

        ! Retrieve the current face's signed normal vector.
        normal = self % faces % getFaceNormal(faceIdxs(j))

        ! If the dot product is negative increment count and cycle to the next face.
        if (dot_product(u, normal) <= ZERO) then
          count = count + 1
          cycle

        end if
        
        ! If the test fails and the current face is a boundary face, the particle is outside the mesh
        ! and we can return early.
        if (self % faces % getFaceIsBoundary(absFaceIdxs(j))) return

      end do

      ! If count > maxCount, update maxCount and maxIdx.
      if (count > maxCount) then
        maxCount = count
        maxIdx = i

      end if

    end do

    ! If reached here, update elementIdx.
    elementIdx = elementIdxs(maxIdx)

  end function findElementFromDirection

  !! Function 'getAllVertexCoordinates'
  !!
  !! Basic description:
  !!   Returns the 3-D coordinates of all the vertices in the mesh.
  !!
  !! Result:
  !!   coords -> 3-D coordinates of all the vertices in the mesh.
  !!
  pure function getAllVertexCoordinates(self) result(coords)
    class(unstructuredMesh), intent(in)           :: self
    real(defReal), dimension(self % nVertices, 3) :: coords

    coords = self % vertices % getAllCoordinates()

  end function getAllVertexCoordinates

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
    call self % edges % kill()
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
      nFaces = size(self % elements % getElementFaceIdxs(i))
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

  !! Subroutine 'setEdgeShelf'
  !!
  !! Basic description:
  !!   Sets the edgeShelf of the unstructured mesh.
  !!
  !! Arguments:
  !!   edges [in] -> An edgeShelf.
  !!
  elemental subroutine setEdgeShelf(self, edges)
    class(unstructuredMesh), intent(inout) :: self
    type(edgeShelf), intent(in)            :: edges

    self % edges = edges

  end subroutine setEdgeShelf

  !! Subroutine 'setElementShelf'
  !!
  !! Basic description:
  !!   Sets the elementShelf of the unstructured mesh.
  !!
  !! Arguments:
  !!   elements [in] -> An elementShelf.
  !!
  elemental subroutine setElementShelf(self, elements)
    class(unstructuredMesh), intent(inout) :: self
    type(elementShelf), intent(in)         :: elements

    self % elements = elements

  end subroutine setElementShelf

  !! Subroutine 'setFaceShelf'
  !!
  !! Basic description:
  !!   Sets the faceShelf of the unstructured mesh.
  !!
  !! Arguments:
  !!   faces [in] -> A faceShelf.
  !!
  elemental subroutine setFaceShelf(self, faces)
    class(unstructuredMesh), intent(inout) :: self
    type(faceShelf), intent(in)            :: faces

    self % faces = faces

  end subroutine setFaceShelf

  !! Subroutine 'setVertexShelf'
  !!
  !! Basic description:
  !!   Sets the vertexShelf of the unstructured mesh.
  !!
  !! Arguments:
  !!   vertices [in] -> A vertexShelf.
  !!
  elemental subroutine setVertexShelf(self, vertices)
    class(unstructuredMesh), intent(inout) :: self
    type(vertexShelf), intent(in)          :: vertices

    self % vertices = vertices

  end subroutine setVertexShelf

end module unstructuredMesh_inter