module unstructuredMesh_inter

  use numPrecision
  use universalVariables
  use genericProcedures,   only : append, findCommon, findDifferent, numToChar, removeDuplicates
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
    procedure                           :: findElementFromEdge
    procedure                           :: findElementFromFace
    procedure                           :: findElementFromVertex
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
  !!
  pure subroutine distanceToBoundaryFace(self, d, coords, parentIdx)
    class(unstructuredMesh), intent(in) :: self
    real(defReal), intent(out)          :: d
    type(coord), intent(inout)          :: coords
    integer(shortInt), intent(out)      :: parentIdx
    
    ! Initialise parentIdx = 0 then search the tree for the intersected boundary face. Update parentIdx only
    ! if boundary intersection has been found.
    parentIdx = 0
    call self % tree % findIntersectedFace(d, coords, self % vertices, self % faces)
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
  pure subroutine distanceToNextFace(self, d, coords)
    class(unstructuredMesh), intent(in)          :: self
    real(defReal), intent(out)                   :: d
    type(coord), intent(inout)                   :: coords
    real(defReal), dimension(3)                  :: r, rEnd
    type(element)                                :: currentElement
    integer(shortInt), dimension(:), allocatable :: potentialFaces, faceToElements
    integer(shortInt)                            :: intersectedFaceIdx
    real(defReal)                                :: lambda
    type(face)                                   :: intersectedFace
    
    ! Initialise d = INF, retrieve the element currently occupied by the particle and compute potential 
    ! face intersections.
    d = INF
    currentElement = self % elements % shelf(coords % elementIdx)
    rEnd = coords % rEnd
    call currentElement % computePotentialFaces(rEnd, self % faces, potentialFaces)

    ! If no potential intersections are detected return early.
    if (size(potentialFaces) == 0) return
    
    ! If reached here, compute which face is actually intersected and update d.
    r = coords % r
    call currentElement % computeIntersectedFace(r, rEnd, potentialFaces, intersectedFaceIdx, lambda, &
                                                 self % faces)
    d = norm2(min(ONE, max(ZERO, lambda)) * (rEnd - r))
    
    ! If the intersected face is a boundary face then the particle is leaving the mesh.
    intersectedFace = self % faces % shelf(intersectedFaceIdx)
    if (intersectedFace % getIsBoundary()) then
      coords % elementIdx = 0
      coords % localId = 1
      return

    end if

    ! Else, retrieve the elements sharing the intersected face from mesh connectivity then
    ! update elementIdx and localId.
    faceToElements = intersectedFace % getFaceToElements()
    coords % elementIdx = findDifferent(faceToElements, currentElement % getIdx())
    coords % localId = self % cellZones % findCellZone(coords % elementIdx)

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
        if (failedFace % getIsBoundary()) return
        
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
            elementIdx = self % findElementFromFace(u, currentElement % getIdx(), zeroDotProductFaceIdxs)
          ! Particle is on a element's edge.
          case(2)
            elementIdx = self % findElementFromEdge(u, zeroDotProductFaceIdxs)
          ! Particle is on a element's vertex.
          case default
            elementIdx = self % findElementFromVertex(u, zeroDotProductFaceIdxs)

        end select
        parentIdx = elementIdx
        return
      
      end if
      ! If reached this point the particle is in the current element. Update elementIdx and parentIdx and return.
      elementIdx = currentElement % getIdx()
      parentIdx = elementIdx
      return

    end do searchLoop

  end subroutine findElementAndParentIdxs

  !! Function 'findElementFromEdge'
  !!
  !! Basic description:
  !!   Returns the index of the element occupied by a particle in case the particle lies on an
  !!   element edge.
  !!
  !! Detailed description:
  !!   Generalisation of 'findTetrahedronFromFace' applied to an edge. In this case the element
  !!   assigned to the particle is that for which the number of negative dot products between the
  !!   element's faces containing the edge and the particle's direction is the greatest.
  !!
  !! Arguments:
  !!   u [in]                   -> Particle's direction.
  !!   zeroDotProductFaces [in] -> Indices of the faces on which the particle lies.
  !!
  !! Result:
  !!   elementIdx               -> Index of the tetrahedron occupied by the particle.
  !!
  !! Notes:
  !!   TODO: Rework this. Import edges during mesh construction.
  !!
  pure function findElementFromEdge(self, u, zeroDotProductFaces) result(elementIdx)
    class(unstructuredMesh), intent(in)                      :: self
    real(defReal), dimension(3), intent(in)                  :: u
    integer(shortInt), dimension(:), allocatable, intent(in) :: zeroDotProductFaces
    integer(shortInt)                                        :: i, j, elementIdx, count, maxCount, maxIdx
    integer(shortInt), dimension(:), allocatable             :: testFaceIdxs, edgeVertices, &
                                                                potentialElements, &
                                                                absFaceIdxs, faceIdxs
    type(element)                                            :: currentElement
    type(face)                                               :: currentFace
    real(defReal), dimension(3)                              :: normal
    
    ! Initialise elementIdx = 0 and retrieve the vertices in the common edge.
    elementIdx = 0
    testFaceIdxs = abs(zeroDotProductFaces)
    edgeVertices = findCommon(self % faces % shelf(testFaceIdxs(1)) % getVertices(), &
                              self % faces % shelf(testFaceIdxs(2)) % getVertices())
    
    ! Loop through all faces and retrieve the indices of all faces sharing the common edge.
    do i = 1, self % nFaces
      ! If the index of the 'do' loop corresponds to one of the original faces sharing the common edge
      ! there is no need to search so cycle to the next face.
      if (any(testFaceIdxs == i)) cycle
      
      ! Check if the current face contains the common edge and if so append the index to testFaceIdxs.
      if (size(findCommon(self % faces % shelf(i) % getVertices(), edgeVertices)) == 2) call append(testFaceIdxs, i)

    end do
    
    ! Loop through all the faces containing the common edge and retrieve their associated elements.
    ! We might introduce duplicates in the list of potential elements so remove them.
    do i = 1, size(testFaceIdxs)
      currentFace = self % faces % shelf(testFaceIdxs(i))
      call append(potentialElements, currentFace % getFaceToElements()) 

    end do
    potentialElements = removeDuplicates(potentialElements)
    
    ! Initialise maxCount = 0 then loop through all potential elements.
    maxCount = 0
    do i = 1, size(potentialElements)
      ! Retrieve the faces making the current element and convert them to absolute indices.
      currentElement = self % elements % shelf(potentialElements(i))
      faceIdxs = currentElement % getFaces()
      absFaceIdxs = abs(faceIdxs)
      
      ! Initialise count = 0 then loop through all faces in the current element.
      count = 0
      do j = 1, size(faceIdxs)

        ! If the current face does not contain the common edge cycle to the next face.
        if (.not. any(testFaceIdxs == absFaceIdxs(j))) cycle

        ! Retrieve the current face's signed normal vector and 
        currentFace = self % faces % shelf(absFaceIdxs(j))
        normal = currentFace % getNormal(faceIdxs(j))

        ! If the dot product is negative increment count and cycle to the next face.
        if (dot_product(u, normal) <= ZERO) then
          count = count + 1
          cycle

        end if
        
        ! If the test fails and the current face is a boundary face, the particle is outside the mesh
        ! and we can return early.
        if (currentFace % getIsBoundary()) return

      end do

      ! If count > maxCount, update maxCount and maxIdx.
      if (count > maxCount) then
        maxCount = count
        maxIdx = i

      end if

    end do

    ! If reached here, update elementIdx.
    elementIdx = potentialElements(maxIdx)

  end function findElementFromEdge

  !! Function 'findElementFromFace'
  !!
  !! Basic description:
  !!   Returns the index of the element occupied by the particle in case the particle lies on an
  !!   element face.
  !!
  !! Detailed description:
  !!   When a particle is on an element face it is not as straightforward to assign an element to
  !!   a particle (it is actually in both elements at the same time). The workaround here is to
  !!   use the particle's direction to determine into which element the particle's direction 
  !!   points. This is given by the dot product of the direction and the face's normal being
  !!   negative.
  !!
  !! Arguments:
  !!   u [in]                   -> Particle's direction.
  !!   currentElementIdx [in]   -> Index of the current element being searched.
  !!   zeroDotProductFaces [in] -> Indices of the triangles on which the particle lies.
  !!
  !! Result:
  !!   elementIdx               -> Index of the element occupied by the particle.
  !!
  pure function findElementFromFace(self, u, currentElementIdx, zeroDotProductFaces) result(elementIdx)
    class(unstructuredMesh), intent(in)                      :: self
    real(defReal), dimension(3), intent(in)                  :: u
    integer(shortInt), intent(in)                            :: currentElementIdx
    integer(shortInt), dimension(:), allocatable, intent(in) :: zeroDotProductFaces
    integer(shortInt)                                        :: absZeroDotProductFace, elementIdx
    integer(shortInt), dimension(:), allocatable             :: faceToElements
    real(defReal), dimension(3)                              :: normal
    type(face)                                               :: currentFace
    
    ! Initialise elementIdx = currentElementIdx then create absolute indices and retrieve the 
    ! face's normal vector.
    elementIdx = currentElementIdx
    absZeroDotProductFace = abs(zeroDotProductFaces(1))
    currentFace = self % faces % shelf(absZeroDotProductFace)
    normal = currentFace % getNormal(zeroDotProductFaces(1))
    
    ! Check the sign of the dot product between the particle's direction and the face's normal.
    ! If negative it is in the current tetrahedron and we can return early.
    if (dot_product(u, normal) <= ZERO) return

    ! If not, the particle's in the neighbouring element. Update elementIdx = 0.
    elementIdx = 0

    ! If the current face is a boundary face return. Else, retrieve the elements sharing the face
    ! from mesh connectivity information and update elementIdx to the neighbouring element.
    if (currentFace % getIsBoundary()) return
    faceToElements = currentFace % getFaceToElements()
    elementIdx = findDifferent(faceToElements, currentElementIdx)

  end function findElementFromFace

  !! Function 'findElementFromVertex'
  !!
  !! Basic description:
  !!   Returns the index of the element occupied by a particle in case the particle lies on an
  !!   element's vertex.
  !!
  !! Detailed description:
  !!   Generalisation of 'findElementFromFace' applied to a vertex. In this case the element
  !!   assigned to the particle is that for which the number of negative dot products between the
  !!   element's faces containing the vertex and the particle's direction is the greatest.
  !!
  !! Arguments:
  !!   direction [in]           -> Particle's direction.
  !!   zeroDotProductFaces [in] -> Indices of the faces on which the particle lies.
  !!
  !! Result:
  !!   elementIdx               -> Index of the element occupied by the particle.
  !!
  pure function findElementFromVertex(self, u, zeroDotProductFaces) result(elementIdx)
    class(unstructuredMesh), intent(in)                      :: self
    real(defReal), dimension(3), intent(in)                  :: u
    integer(shortInt), dimension(:), allocatable, intent(in) :: zeroDotProductFaces
    integer(shortInt)                                        :: commonVertex, i, j, count, maxCount, maxIdx, elementIdx
    integer(shortInt), dimension(:), allocatable             :: absFaces, absZeroDotProductFaces, &
                                                                commonVertices, potentialElements, faces
    type(face)                                               :: currentFace
    real(defReal), dimension(3)                              :: normal
    
    ! Initialise elementIdx = 0. Find the common vertex between the faces and retrieve all 
    ! elements sharing this vertex.
    elementIdx = 0
    absZeroDotProductFaces = abs(zeroDotProductFaces)
    commonVertices = self % faces % shelf(absZeroDotProductFaces(1)) % getVertices()
    do i = 2, size(zeroDotProductFaces)
      commonVertices = findCommon(commonVertices, self % faces % shelf(absZeroDotProductFaces(i)) % getVertices())

    end do
    commonVertex = commonVertices(1)
    potentialElements = self % vertices % shelf(commonVertex) % getVertexToElements()
    
    ! Initialise maxCount = 0 then loop through all potential elements.
    maxCount = 0
    do i = 1, size(potentialElements)
      ! Retrieve the faces in the current element and create absolute indices.
      faces = self % elements % shelf(potentialElements(i)) % getFaces()
      absFaces = abs(faces)
      
      ! Initialise count = 0 then loop through all the faces in the current element.
      count = 0
      do j = 1, size(faces)
        currentFace = self % faces % shelf(absFaces(j))
        ! Cycle if the current face does not contain the common vertex.
        if (.not. any(currentFace % getVertices() == commonVertex)) cycle
        
        ! Retrieve the current face's signed normal vector.
        normal = currentFace % getNormal(faces(j))
        
        ! If the dot product is negative increment the number of faces for the current 
        ! element and cycle to the next face.
        if (dot_product(u, normal) <= ZERO) then
          count = count + 1
          cycle
  
        end if
        
        ! If the test fails and the current face is a boundary face, the particle is outside the mesh
        ! and we can return early.
        if (currentFace % getIsBoundary()) return

      end do

      ! If count > maxCount, update maxCount and maxIdx.
      if (count > maxCount) then
        maxCount = count
        maxIdx = i

      end if

    end do

    ! If reached here, update elementIdx.
    elementIdx = potentialElements(maxIdx)
    
  end function findElementFromVertex

end module unstructuredMesh_inter