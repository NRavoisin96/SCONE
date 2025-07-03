module unstructuredMesh_inter

  use accelerationStructure_inter, only : accelerationStructure
  use coord_class,                 only : coord
  use dictionary_class,            only : dictionary
  use edge_class,                  only : edgeBox
  use edgeShelf_class,             only : edgeShelf
  use element_class,               only : buildElementInfo, elementBox, inclusionTestResult
  use elementShelf_class,          only : elementShelf
  use face_class,                  only : buildFaceInfo, faceBox
  use faceShelf_class,             only : faceShelf
  use genericProcedures,           only : append, findDifferent, numToChar, removeDuplicates
  use mesh_inter,                  only : mesh, kill_super => kill
  use numPrecision
  use octreeAcceleration_class,    only : octreeAcceleration
  use publicObjects,               only : basicEdgeInfo, basicElementInfo, basicFaceInfo, buildEdgeInfo, meshLocalIdInfo
  use universalVariables
  use vertex_class,                only : vertexBox
  use vertexShelf_class,           only : vertexShelf

  implicit none
  private

  ! Extendable procedures.
  public :: distanceToBoundaryFace, distanceToNextFace, findHostElement, kill
  
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
  type, public, abstract, extends(mesh)       :: unstructuredMesh
    private
    integer(shortInt)                         :: nVertices = 0, nFaces = 0, nEdges = 0, &
                                                 nElements = 0, nInternalFaces = 0
    class(accelerationStructure), allocatable :: acceleration
    type(edgeShelf)                           :: edges
    type(elementShelf)                        :: elements
    type(faceShelf)                           :: faces
    type(vertexShelf)                         :: vertices
  contains
    ! Build procedures.
    procedure                       :: assignLocalIds
    procedure(importMesh), deferred :: importMesh
    procedure                       :: init
    procedure                       :: initEdgeShelf
    procedure                       :: initElementShelf
    procedure                       :: initFaceShelf
    procedure                       :: initVertexShelf
    procedure                       :: kill
    procedure, non_overridable      :: printComposition
    procedure                       :: split
    procedure                       :: splitElements
    procedure                       :: splitFaces
    ! Runtime procedures.
    procedure                       :: distanceToBoundaryFace
    procedure                       :: distanceToNextFace
    procedure                       :: findHostElement
    procedure                       :: getAllVertexCoordinates
    procedure                       :: getEdgesNumber
    procedure                       :: getFacesNumber
    procedure                       :: getElementsNumber
    procedure                       :: getInternalFacesNumber
    procedure                       :: getVerticesNumber
  end type unstructuredMesh

  abstract interface

    !! Subroutine 'distanceToNextFace'
    !!
    !! Basic description:
    !!   Returns the distance to the next intersected mesh face.
    !!
    !! Arguments:
    !!   d [out]        -> Distance to the next intersected face.
    !!   coords [inout] -> Particle's coordinates.
    !!
    subroutine importMesh(self, folderPath)
      import                                 :: unstructuredMesh
      class(unstructuredMesh), intent(inout) :: self
      character(*), intent(in)               :: folderPath
    end subroutine importMesh

  end interface

contains
  !!
  !!
  !!
  subroutine assignLocalIds(self, localIdInfos)
    class(unstructuredMesh), intent(inout)          :: self
    type(meshLocalIdInfo), dimension(:), intent(in) :: localIdInfos
    integer(shortInt)                               :: i, nLocalIds

    nLocalIds = size(localIdInfos)
    call self % setLocalIdsNumber(nLocalIds)
    do i = 1, nLocalIds
      call self % elements % setElementLocalId(localIdInfos(i) % elementIdxs, localIdInfos(i) % localId)

    end do

  end subroutine assignLocalIds

  !! Subroutine 'distanceToBoundaryFace'
  !!
  !! Basic description:
  !!   Returns the distance to the mesh boundary face intersected by a particle's path. Also returns the index
  !!   of the parent element containing the intersected boundary face.
  !!
  !! See mesh_inter for details.
  !!
  subroutine distanceToBoundaryFace(self, d, coords)
    class(unstructuredMesh), intent(in)          :: self
    real(defReal), intent(out)                   :: d
    type(coord), intent(inout)                   :: coords
    integer(shortInt)                            :: i, boundaryFaceIdx
    real(defReal)                                :: update
    integer(shortInt), dimension(:), allocatable :: elementIdxs
    type(inclusionTestResult)                    :: testResult
    
    ! Initialise parentIdx = 0, edgeIdx = 0 and vertexIdx = 0 then search the tree for the intersected boundary face.
    d = INF
    call coords % setParentElementIdx(0)
    boundaryFaceIdx = 0
    do i = 1, self % nFaces
      if (.not. self % faces % getFaceIsBoundary(i)) cycle

      ! Compute distance to boundary face.
      call self % faces % computeFaceIntersection(i, coords, update)
      if (update < d) then
        d = update
        boundaryFaceIdx = i

      end if

    end do

    ! Retrieve the element associated with the boundary face.
    if (boundaryFaceIdx > 0) then
      elementIdxs = self % faces % getFaceElementIdxs(boundaryFaceIdx)
      call coords % setElementIdx(elementIdxs(1))
      call coords % setParentElementIdx(self % elements % getElementParentIdx(elementIdxs(1)))

      ! Set coords % endPosition to the minimum computed distance plus a slight forward nudge.
      call coords % setEndPosition(coords % getPosition() + (d + NUDGE) * coords % getDirection())

      ! If the element associated with the intersected face does not contain the end position, begin rescue.
      testResult = self % elements % isPointInside(elementIdxs(1), coords % getEndPosition())
      if (.not. testResult % status == INSIDE_ELEMENT) then
        call coords % setNudgeEndPosition(.true.)
        call self % findHostElement(coords)
        ! Update distance.
        if (coords % getElementIdx() == 0) then
          d = INF

        else
          d = norm2(coords % getEndPosition() - coords % getPosition())

        end if
        call coords % setNudgeEndPosition(.false.)

      end if

    end if

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
  subroutine distanceToNextFace(self, d, coords)
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
    elementIdx = coords % getElementIdx()
    rEnd = coords % getEndPosition()
    potentialFaces = self % elements % computePotentialFaceIdxs(elementIdx, rEnd)

    ! If no potential intersections are detected return early.
    if (size(potentialFaces) == 0) return
    
    ! If reached here, compute which face is actually intersected and update d.
    r = coords % getPosition()
    call self % elements % computeFaceIntersection(elementIdx, r, rEnd, potentialFaces, intersectedFaceIdx, lambda)
    d = norm2(min(ONE, max(ZERO, lambda)) * (rEnd - r))
    
    ! If the intersected face is a boundary face then the particle is leaving the mesh.
    if (self % faces % getFaceIsBoundary(intersectedFaceIdx)) then
      call coords % setElementIdx(0)
      call coords % setLocalId(1)
      return

    end if

    ! Else, retrieve the elements sharing the intersected face from mesh connectivity then
    ! update elementIdx and localId.
    faceToElements = self % faces % getFaceElementIdxs(intersectedFaceIdx)
    call coords % setElementIdx(findDifferent(faceToElements, elementIdx))
    call coords % setLocalId(self % elements % getElementLocalId(coords % getElementIdx()))

  end subroutine distanceToNextFace

  !! Subroutine 'findElementAndParentIdxs'
  !!
  !! Basic description:
  !!   Returns the index of the mesh element occupied by a particle. Also returns the index of the parent mesh
  !!   element containing the occupied element.
  !!
  !! See mesh_inter for details.
  !!
  subroutine findHostElement(self, coords)
    class(unstructuredMesh), intent(in)          :: self
    type(coord), intent(inout)                   :: coords
    integer(shortInt), dimension(:), allocatable :: potentialElementsIdxs
    integer(shortInt)                            :: i, nPotentialElements, potentialElementIdx
    real(defReal), dimension(3)                  :: r
    type(inclusionTestResult)                    :: testResult
    
    ! Initialise parentIdx = 0. Retrieve the mesh's bounding box. If the particle is outside the bounding box we can return early.
    call coords % setElementIdx(0)
    call coords % setParentElementIdx(0)
    if (allocated(self % acceleration)) then
      call self % acceleration % findHostElement(self % elements, coords)

    else
      ! Perform brute-force search.
      searchLoop: do
        do i = 1, self % nElements
          testResult = self % elements % isPointInside(i, coords % getPositionToNudge())
          if (testResult % status == INSIDE_ELEMENT) then
            call coords % setElementIdx(i)
            call coords % setParentElementIdx(self % elements % getElementParentIdx(i))
            call coords % setLocalId(self % elements % getElementLocalId(i))
            return

          elseif (testResult % status == ON_BOUNDARY_ELEMENT) then
            ! If coordinates are on the element boundary (very rare), we need to push them off.
            do while (testResult % status == ON_BOUNDARY_ELEMENT)
              call self % elements % pushFromElementBoundary(i, coords)

              ! Perform containment test again.
              testResult = self % elements % isPointInside(i, coords % getPositionToNudge())

            end do

            ! Now the coordinates are not on the boundary of the element anymore.
            if (testResult % status == INSIDE_ELEMENT) then
              ! If coordinates are now well inside the element, we have found our element.
              call coords % setElementIdx(i)
              call coords % setParentElementIdx(self % elements % getElementParentIdx(i))
              call coords % setLocalId(self % elements % getElementLocalId(i))
              return

            elseif (testResult % status == OUTSIDE_ELEMENT) then
              ! If the nudge has resulted in an overshoot, we cycle searchLoop and begin the entire process again.
              cycle searchLoop

            end if

          end if

        end do
        return

      end do searchLoop

    end if

  end subroutine findHostElement

  !! Function 'getAllVertexCoordinates'
  !!
  !! Basic description:
  !!   Returns the 3-D coordinates of all the vertices in the mesh.
  !!
  !! Result:
  !!   coords -> 3-D coordinates of all the vertices in the mesh.
  !!
  function getAllVertexCoordinates(self) result(coords)
    class(unstructuredMesh), intent(in)           :: self
    real(defReal), dimension(3, self % nVertices) :: coords

    coords = self % vertices % getAllCoordinates()

  end function getAllVertexCoordinates

  !!
  !!
  !!
  elemental function getEdgesNumber(self) result(nEdges)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt)                   :: nEdges

    nEdges = self % nEdges

  end function getEdgesNumber

  !!
  !!
  !!
  elemental function getElementsNumber(self) result(nElements)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt)                   :: nElements

    nElements = self % nElements

  end function getElementsNumber

  !!
  !!
  !!
  elemental function getFacesNumber(self) result(nFaces)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt)                   :: nFaces

    nFaces = self % nFaces

  end function getFacesNumber

  !!
  !!
  !!
  elemental function getInternalFacesNumber(self) result(nInternalFaces)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt)                   :: nInternalFaces

    nInternalFaces = self % nInternalFaces

  end function getInternalFacesNumber

  !!
  !!
  !!
  elemental function getVerticesNumber(self) result(nVertices)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt)                   :: nVertices

    nVertices = self % nVertices

  end function getVerticesNumber

  !!
  !!
  !!
  subroutine init(self, folderPath, dict)
    class(unstructuredMesh), intent(inout) :: self
    character(*), intent(in)               :: folderPath
    class(dictionary), intent(in)          :: dict
    logical(defBool)                       :: triangulate
    character(nameLen)                     :: acceleration

    ! Set up base components.
    call self % setupBase(dict)
    
    ! Import mesh from files.
    call self % importMesh(folderPath)

    ! Check if triangulation was requested.
    call dict % getOrDefault(triangulate, 'triangulate', .false.)

    ! Check if acceleration structure was required by user and initialise it if applicable.
    call dict % getOrDefault(acceleration, 'accelerationMethod', 'none')
    if (acceleration /= 'none') then
      if (acceleration == 'octree') allocate(octreeAcceleration :: self % acceleration)
      call self % acceleration % init(self % vertices, self % faces, self % elements)

    end if

  end subroutine init

  !!
  !!
  !!
  subroutine initEdgeShelf(self, edgeInfos)
    class(unstructuredMesh), intent(inout)          :: self
    type(basicEdgeInfo), dimension(:), intent(in)   :: edgeInfos
    type(buildEdgeInfo), dimension(size(edgeInfos)) :: buildInfos
    integer(shortInt)                               :: i, nEdges

    nEdges = size(edgeInfos)
    self % nEdges = nEdges
    do i = 1, nEdges
      buildInfos(i) % idx = edgeInfos(i) % idx
      buildInfos(i) % vertices = self % vertices % getVertexBox(edgeInfos(i) % vertexIdxs)
      
      ! Update connectivity.
      call self % vertices % addEdgeIdxToVertex(edgeInfos(i) % vertexIdxs, edgeInfos(i) % idx)

    end do
    call self % edges % init(buildInfos)

  end subroutine initEdgeShelf

  !!
  !!
  !!
  subroutine initElementShelf(self, elementInfos)
    class(unstructuredMesh), intent(inout)                :: self
    type(basicElementInfo), dimension(:), intent(inout)   :: elementInfos
    type(buildElementInfo), dimension(size(elementInfos)) :: buildInfos
    integer(shortInt)                                     :: absFaceIdx, i, j, nElements, nFaces

    nElements = size(elementInfos)
    self % nElements = nElements
    do i = 1, nElements
      buildInfos(i) % idx = elementInfos(i) % idx
      buildInfos(i) % parentIdx = elementInfos(i) % parentIdx
      nFaces = size(elementInfos(i) % faceIdxs)
      allocate(buildInfos(i) % orientatedFaces(nFaces))
      do j = 1, nFaces
        ! Create absFaceIdx then retrieve face properties.
        absFaceIdx = abs(elementInfos(i) % faceIdxs(j))
        buildInfos(i) % orientatedFaces(j) % face = self % faces % getFaceBox(absFaceIdx)
        buildInfos(i) % orientatedFaces(j) % outwardNormal = self % faces % getFaceNormal(elementInfos(i) % faceIdxs(j))

        ! Update connectivity.
        call self % faces % addElementIdxToFace(absFaceIdx, elementInfos(i) % idx)

      end do

      ! Check if we need to construct edges.
      if (.not. allocated(elementInfos(i) % edgeIdxs)) call createEdgeIdxs(elementInfos(i))
      buildInfos(i) % edges = self % edges % getEdgeBox(elementInfos(i) % edgeIdxs)

      ! Update connectivity.
      call self % edges % addElementIdxToEdge(elementInfos(i) % edgeIdxs, elementInfos(i) % idx)

      ! Check if we need to construct vertices.
      if (.not. allocated(elementInfos(i) % vertexIdxs)) call createVertexIdxs(elementInfos(i))
      buildInfos(i) % vertices = self % vertices % getVertexBox(elementInfos(i) % vertexIdxs)

      ! Update connectivity.
      call self % vertices % addElementIdxToVertex(elementInfos(i) % vertexIdxs, elementInfos(i) % idx)

    end do
    call self % elements % init(buildInfos)
  
  contains
    !!
    !!
    !!
    subroutine createEdgeIdxs(info)
      type(basicElementInfo), intent(inout)                 :: info
      logical(defBool), dimension(self % edges % getSize()) :: isPresent
      integer(shortInt)                                     :: currentSize, idx, k, l, nEdges
      type(edgeBox), dimension(:), allocatable              :: faceEdges
      integer(shortInt), dimension(:), allocatable          :: tempIdxs

      ! Initialise isPresent = .false. and allocate info % edgeIdxs to an appropriate initial size.
      isPresent = .false.
      allocate(info % edgeIdxs(6))
      nEdges = 0
      do k = 1, size(info % faceIdxs)
        faceEdges = self % faces % getFaceEdges(abs(info % faceIdxs(k)))
        do l = 1, size(faceEdges)
          idx = faceEdges(l) % ptr % getIdx()
          if (.not. isPresent(idx)) then
            nEdges = nEdges + 1
            currentSize = size(info % edgeIdxs)
            if (currentSize < nEdges) then
              allocate(tempIdxs(2 * currentSize))
              tempIdxs(1:currentSize) = info % edgeIdxs
              call move_alloc(tempIdxs, info % edgeIdxs)

            end if
            info % edgeIdxs(nEdges) = idx
            isPresent(idx) = .true.

          end if

        end do

      end do

      ! Resize info % edgeIdxs to correct size if necessary.
      if (nEdges < size(info % edgeIdxs)) then
        allocate(tempIdxs(nEdges))
        tempIdxs = info % edgeIdxs(1:nEdges)
        call move_alloc(tempIdxs, info % edgeIdxs)

      end if

    end subroutine createEdgeIdxs

    !!
    !!
    !!
    subroutine createVertexIdxs(info)
      type(basicElementInfo), intent(inout)                    :: info
      logical(defBool), dimension(self % vertices % getSize()) :: isPresent
      integer(shortInt)                                        :: currentSize, idx, k, l, nVertices
      type(vertexBox), dimension(:), allocatable               :: faceVertices
      integer(shortInt), dimension(:), allocatable             :: tempIdxs

      ! Initialise isPresent = .false. and allocate info % vertexIdxs to an appropriate initial size.
      isPresent = .false.
      allocate(info % vertexIdxs(4))
      nVertices = 0
      do k = 1, size(info % faceIdxs)
        faceVertices = self % faces % getFaceVertices(abs(info % faceIdxs(k)))
        do l = 1, size(faceVertices)
          idx = faceVertices(l) % ptr % getIdx()
          if (.not. isPresent(idx)) then
            nVertices = nVertices + 1
            currentSize = size(info % vertexIdxs)
            if (currentSize < nVertices) then
              allocate(tempIdxs(2 * currentSize))
              tempIdxs(1:currentSize) = info % vertexIdxs
              call move_alloc(tempIdxs, info % vertexIdxs)

            end if
            info % vertexIdxs(nVertices) = idx
            isPresent(idx) = .true.

          end if

        end do

      end do

      ! Resize info % vertexIdxs to correct size if necessary.
      if (nVertices < size(info % vertexIdxs)) then
        allocate(tempIdxs(nVertices))
        tempIdxs = info % vertexIdxs(1:nVertices)
        call move_alloc(tempIdxs, info % vertexIdxs)

      end if

    end subroutine createVertexIdxs

  end subroutine initElementShelf

  !!
  !!
  !!
  subroutine initFaceShelf(self, faceInfos)
    class(unstructuredMesh), intent(inout)          :: self
    type(basicFaceInfo), dimension(:), intent(in)   :: faceInfos
    type(buildFaceInfo), dimension(size(faceInfos)) :: buildInfos
    integer(shortInt)                               :: i, nFaces

    nFaces = size(faceInfos)
    self % nFaces = nFaces
    do i = 1, nFaces
      buildInfos(i) % idx = faceInfos(i) % idx
      buildInfos(i) % parentIdx = faceInfos(i) % parentIdx
      buildInfos(i) % isBoundary = faceInfos(i) % isBoundary
      buildInfos(i) % vertices = self % vertices % getVertexBox(faceInfos(i) % vertexIdxs)
      buildInfos(i) % edges = self % edges % getEdgeBox(faceInfos(i) % edgeIdxs)

      ! Update connectivity.
      call self % vertices % addFaceIdxToVertex(faceInfos(i) % vertexIdxs, faceInfos(i) % idx)
      call self % edges % addFaceIdxToEdge(faceInfos(i) % edgeIdxs, faceInfos(i) % idx)

    end do
    call self % faces % init(buildInfos)

  end subroutine initFaceShelf

  !!
  !!
  !!
  subroutine initVertexShelf(self, coords)
    class(unstructuredMesh), intent(inout)     :: self
    real(defReal), dimension(:, :), intent(in) :: coords

    self % nVertices = size(coords, 2)
    call self % vertices % init(coords)
    call self % initBoundingBox(self % vertices % getExtremalCoordinates())

  end subroutine

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an unitialised state.
  !!
  subroutine kill(self)
    class(unstructuredMesh), intent(inout) :: self

    ! Superclass.
    call kill_super(self)
    
    ! Local.
    self % nVertices = 0
    self % nFaces = 0
    self % nInternalFaces = 0
    self % nElements = 0
    self % nEdges = 0
    call self % elements % kill()
    call self % faces % kill()
    call self % edges % kill()
    call self % vertices % kill()
    if (allocated(self % acceleration)) then
      call self % acceleration % kill()
      deallocate(self % acceleration)

    end if

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
      nFaces = size(self % elements % getElementFaces(i))
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

  !! Subroutine 'split'
  !!
  !! Basic description:
  !!   Splits a mesh into tetrahedral elements. If a given element is already a tetrahedron it is
  !!   not split but simply added to the shelf of tetrahedra in the mesh.
  !!
  !! Detailed description:
  !!   'split' starts by computing the number of pyramids, tetrahedra and triangles that will be
  !!   generated in the resulting mesh. Then, each element is split into a set of pyramids, whose
  !!   bases are each of the element's face and whose (common) apex is the element's centroid. This
  !!   apex is also appended to the list of vertices in the mesh in the process. Once this is done,
  !!   each face in the original mesh is subdivided into triangles. Lastly, each pyramid previously
  !!   created is further split into tetrahedra.
  !!
  !! Arguments:
  !!   lastVertexIdx [out] -> Index of the last vertex in the resulting mesh.
  !!
  subroutine split(self, edges, elements, faces, vertices, newEdges, newElements, newFaces, newVertices)
    class(unstructuredMesh), intent(inout)    :: self
    type(edgeShelf), intent(in)               :: edges
    type(elementShelf), intent(inout)         :: elements
    type(faceShelf), intent(inout)            :: faces
    type(vertexShelf), intent(in)             :: vertices
    type(edgeShelf), intent(out)              :: newEdges
    type(elementShelf), intent(out)           :: newElements
    type(faceShelf), intent(out)              :: newFaces
    type(vertexShelf), intent(out)            :: newVertices
    integer(shortInt)                         :: i, j, nEdges, nInternalTriangles, nNewEdges, nTetrahedra, nTriangles, &
                                                 nVertices, nNewVertices, lastEdgeIdx, lastFaceIdx, lastElementIdx, lastVertexIdx
    integer(shortInt), dimension(:), allocatable :: edgeIdxs
    type(elementBox), dimension(:), allocatable  :: tetrahedra
    type(faceBox), dimension(:), allocatable     :: triangles

  end subroutine split

  !! Subroutine 'splitElements'
  !!
  !! Basic description:
  !!   Splits all elements in the original mesh into pyramids. If a given element is already a
  !!   tetrahedron it is not split but simply appended to the list of existing tetrahedra.
  !!
  !! Arguments:
  !!   lastEdgeIdx [inout]        -> Index of the last edge in the mesh.
  !!   lastPyramidIdx [inout]     -> Index of the last pyramid in the mesh.
  !!   lastTetrahedronIdx [inout] -> Index of the last tetrahedron in the mesh.
  !!   lastTriangleIdx [inout]    -> Index of the last triangle in the mesh.
  !!   lastVertexIdx [inout]      -> Index of the last vertex in the mesh.
  !!
  subroutine splitElements(self, elements, faces, lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx, &
                           newEdges, newElements, newFaces, newVertices, tetrahedra, triangles)
    class(unstructuredMesh), intent(inout)        :: self
    type(elementShelf), intent(inout)             :: elements, newElements
    type(faceShelf), intent(inout)                :: faces, newFaces
    integer(shortInt), intent(inout)              :: lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx
    type(edgeShelf), intent(inout)                :: newEdges
    type(vertexShelf), intent(inout)              :: newVertices
    type(elementBox), dimension(:), intent(inout) :: tetrahedra
    type(faceBox), dimension(:), intent(inout)    :: triangles
    integer(shortInt)                             :: i, initialElementIdx, j

  end subroutine splitElements

  !!
  !!
  !!
  subroutine splitFaces(self, faces, newEdges, newFaces, newVertices, lastNewEdgeIdx, lastNewFaceIdx, triangles)
    class(unstructuredMesh), intent(inout)     :: self
    type(faceShelf), intent(inout)             :: faces, newFaces
    type(edgeShelf), intent(inout)             :: newEdges
    type(vertexShelf), intent(inout)           :: newVertices
    integer(shortInt), intent(inout)           :: lastNewEdgeIdx, lastNewFaceIdx
    type(faceBox), dimension(:), intent(inout) :: triangles
    integer(shortInt)                          :: i, initialFaceIdx, j

  end subroutine splitFaces

end module unstructuredMesh_inter