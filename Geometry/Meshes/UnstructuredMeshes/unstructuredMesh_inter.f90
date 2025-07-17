module unstructuredMesh_inter

  use accelerationStructure_inter,       only : accelerationStructure
  use accelerationStructureFactory_func, only : newAccelerationStructurePtr
  use axisAlignedBoundingBox_class,      only : axisAlignedBoundingBox
  use coord_class,                       only : coord
  use dictionary_class,                  only : dictionary
  use edge_class,                        only : edgeBox
  use element_class,                     only : buildElementPayload, element, elementBox, inclusionTestResult, &
                                                rayIntersectionTestResult
  use extentTopologicalObject_inter,     only : buildExtentTopologicalObjectPayload
  use face_class,                        only : buildFacePayload, face, faceBox
  use genericProcedures,                 only : append, fatalError, numToChar
  use mesh_inter,                        only : mesh, distanceToBoundary_super => distanceToBoundary, &
                                                findHostElement_super => findHostElement, kill_super => kill
  use numPrecision
  use publicObjects,                     only : basicEdgeInfo, basicElementInfo, basicFaceInfo, basicVertexInfo, &
                                                intersectionTestResult, meshLocalIdInfo
  use topologicalObject_inter,           only : topologicalObjectBox
  use topologicalObjectShelf_class,      only : topologicalObjectShelf
  use triangulationFactory_func,         only : newTriangulationPtr
  use triangulationMethod_inter,         only : triangulationMethod
  use universalVariables
  use vertex_class,                      only : buildVertexPayload, vertexBox

  implicit none
  private

  ! Extendable procedures.
  public :: distanceToBoundary, distanceToNextFace, findHostElement, kill
  
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
  type, public, abstract, extends(mesh)   :: unstructuredMesh
    private
    integer(shortInt)                     :: nVertices = 0, nFaces = 0, nEdges = 0, &
                                             nElements = 0, nInternalFaces = 0
    class(accelerationStructure), pointer :: acceleration
    type(topologicalObjectShelf)          :: edges, elements, faces, vertices
    class(triangulationMethod), pointer   :: triangulation
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
    procedure                       :: setEdgesNumber
    procedure                       :: setElementsNumber
    procedure                       :: setFacesNumber
    procedure                       :: setInternalFacesNumber
    procedure                       :: setVerticesNumber
    ! Runtime procedures.
    procedure                       :: distanceToBoundary
    procedure                       :: distanceToNextFace
    procedure                       :: findHostElement
    procedure                       :: getEdgesNumber
    procedure                       :: getElementsNumber
    procedure                       :: getFacesNumber
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
    integer(shortInt)                               :: i, j, nLocalIds
    type(elementBox)                                :: element

    nLocalIds = size(localIdInfos)
    call self % setLocalIdsNumber(nLocalIds)
    do i = 1, nLocalIds
      do j = 1, size(localIdInfos(i) % elementIdxs)
        element = self % elements % getElementBox(localIdInfos(i) % elementIdxs(j))
        call element % ptr % setLocalId(localIdInfos(i) % localId)

      end do

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
  subroutine distanceToBoundary(self, d, coords)
    class(unstructuredMesh), intent(in)                   :: self
    real(defReal), intent(out)                            :: d
    type(coord), intent(inout)                            :: coords
    type(axisAlignedBoundingBox), pointer                 :: boundingBoxPtr
    type(intersectionTestResult)                          :: boundingBoxIntersectionResult
    type(faceBox)                                         :: boundaryFace
    type(topologicalObjectBox), dimension(:), allocatable :: faceElements
    type(inclusionTestResult)                             :: testResult
    character(*), parameter                               :: here = 'distanceToBoundary (unstructuredMesh_inter.f90)'
    
    ! Initialise parentIdx = 0, edgeIdx = 0 and vertexIdx = 0 then search the tree for the intersected boundary face.
    call distanceToBoundary_super(self, d, coords)
    boundingBoxPtr => self % getBoundingBoxPtr()
    boundingBoxIntersectionResult = boundingBoxPtr % intersects(coords % getPosition(), coords % getDirection())
    if (.not. boundingBoxIntersectionResult % intersects) return

    call self % acceleration % findEntranceBoundaryFace(self % faces, coords, d, boundaryFace)
    if (.not. associated(boundaryFace % ptr)) return

    ! Retrieve the element associated with the boundary face.
    faceElements = boundaryFace % ptr % getSharingElements()
    ! Downcast elements to correct type.
    select type(ptr => faceElements(1) % ptr)
      type is(element)
        call coords % setElementIdx(ptr % getIdx())
        call coords % setParentElementIdx(ptr % getParentIdx())
        call coords % setLocalId(ptr % getLocalId())

        ! Set coords % endPosition to the minimum computed distance plus a slight forward nudge.
        call coords % setEndPosition(coords % getPosition() + (d + NUDGE) * coords % getDirection())

        ! If the element associated with the intersected face does not contain the end position, begin rescue.
        testResult = ptr % isPointInside(coords % getEndPosition())
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

      class default
        call fatalError(here, 'Element with index: '//numToChar(ptr % getIdx())//' is not an element.')

    end select

  end subroutine distanceToBoundary

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
    class(unstructuredMesh), intent(in)                   :: self
    real(defReal), intent(out)                            :: d
    type(coord), intent(inout)                            :: coords
    type(elementBox)                                      :: currentElement
    type(rayIntersectionTestResult)                       :: intersectionResults
    integer(shortInt)                                     :: i, nElements
    type(topologicalObjectBox), dimension(:), allocatable :: faceElements
    character(*), parameter                               :: here = 'distanceToNextFace (unstructuredMesh_inter.f90)'
    
    ! Retrieve the element currently occupied by the particle and compute potential 
    ! face intersections.
    d = INF
    currentElement = self % elements % getElementBox(coords % getElementIdx())
    intersectionResults = currentElement % ptr % intersects_Ray(coords % getPosition(), coords % getEndPosition())

    if (.not. intersectionResults % intersects) return
    d = intersectionResults % d
    
    ! If the intersected face is a boundary face then the particle is leaving the mesh.
    if (intersectionResults % intersectedFace % ptr % getIsBoundary()) then
      call coords % setElementIdx(0)
      call coords % setParentElementIdx(0)
      call coords % setLocalId(1)

    else
      ! Else, retrieve the elements sharing the intersected face from mesh connectivity then
      ! update elementIdx and localId.
      faceElements = intersectionResults % intersectedFace % ptr % getSharingElements()
      nElements = size(faceElements)
      if (nElements /= 2) &
      call fatalError(here, 'Internal face: '//numToChar(intersectionResults % intersectedFace % ptr % getIdx())// &
                            ' is not associated to the correct number of elements.')

      do i = 1, 2
        ! Downcast element to correct type.
        select type(ptr => faceElements(i) % ptr)
          type is(element)
            if (.not. associated(currentElement % ptr, ptr)) then
              ! We have found our new element.
              call coords % setElementIdx(ptr % getIdx())
              call coords % setParentElementIdx(ptr % getParentIdx())
              call coords % setLocalId(ptr % getLocalId())

            end if

          class default
            call fatalError(here, 'Element with index: '//numToChar(ptr % getIdx())//' is not an element.')

        end select

      end do

    end if

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
    class(unstructuredMesh), intent(in)   :: self
    type(coord), intent(inout)            :: coords
    type(axisAlignedBoundingBox), pointer :: boundingBoxPtr
    logical(defBool)                      :: stopSearch
    
    ! Set elementIdx = 0, parentElementIdx = 0 and localId = 1.
    call findHostElement_super(self, coords)
    boundingBoxPtr => self % getBoundingBoxPtr()
    
    searchLoop: do
      if (.not. boundingBoxPtr % contains(coords % getPositionToNudge())) return
      call self % acceleration % findHostElement(self % elements, coords, stopSearch)
      if (stopSearch) return

    end do searchLoop

  end subroutine findHostElement

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
    integer(shortInt)                      :: i
    type(elementBox)                       :: element
    type(faceBox)                          :: face

    ! Set up base components.
    call self % setupBase(dict)
    
    ! Import mesh from files.
    call self % importMesh(folderPath)

    ! Initialise triangulation method from dictionary then triangulate mesh.
    call newTriangulationPtr(dict, self % triangulation)
    call self % triangulation % triangulate(self % edges, self % elements, self % faces, self % vertices)

    ! Shrink shelves to their correct size after triangulation.
    call self % edges % shrink()
    call self % elements % shrink()
    call self % faces % shrink()
    call self % vertices % shrink()

    ! Update number of edges, elements, faces, internal faces, and vertices.
    self % nEdges = self % edges % getObjectsNumber()
    self % nVertices = self % vertices % getObjectsNumber()
    self % nElements = 0
    do i = 1, self % elements % getObjectsNumber()
      element = self % elements % getElementBox(i)
      if (element % ptr % getIsActive()) self % nElements = self % nElements + 1

    end do

    self % nFaces = 0
    self % nInternalFaces = 0
    do i = 1, self % faces % getObjectsNumber()
      face = self % faces % getFaceBox(i)
      if (face % ptr % getIsActive()) then
        self % nFaces = self % nFaces + 1
        if (.not. face % ptr % getIsBoundary()) self % nInternalFaces = self % nInternalFaces + 1

      end if

    end do

    ! Initialise acceleration method from dictionary.
    call newAccelerationStructurePtr(dict, self % edges, self % elements, self % faces, self % vertices, self % acceleration)

  end subroutine init

  !!
  !!
  !!
  subroutine initEdgeShelf(self, edgeInfos)
    class(unstructuredMesh), intent(inout)                                :: self
    type(basicEdgeInfo), dimension(:), intent(in)                         :: edgeInfos
    type(buildExtentTopologicalObjectPayload), dimension(size(edgeInfos)) :: payloads
    integer(shortInt)                                                     :: i, nEdges

    nEdges = size(edgeInfos)
    self % nEdges = nEdges
    do i = 1, nEdges
      payloads(i) % idx = edgeInfos(i) % idx
      payloads(i) % vertices = self % vertices % getVertexBox(edgeInfos(i) % vertexIdxs)

    end do
    call self % edges % init(payloads)

  end subroutine initEdgeShelf

  !!
  !!
  !!
  subroutine initElementShelf(self, elementInfos)
    class(unstructuredMesh), intent(inout)                   :: self
    type(basicElementInfo), dimension(:), intent(inout)      :: elementInfos
    type(buildElementPayload), dimension(size(elementInfos)) :: payloads
    integer(shortInt)                                        :: i, j, nElements, nFaces
    type(faceBox)                                            :: face

    nElements = size(elementInfos)
    self % nElements = nElements
    do i = 1, nElements
      payloads(i) % idx = elementInfos(i) % idx
      payloads(i) % parentIdx = elementInfos(i) % parentIdx
      nFaces = size(elementInfos(i) % faceIdxs)
      allocate(payloads(i) % orientatedFaces(nFaces))
      do j = 1, nFaces
        face = self % faces % getFaceBox(abs(elementInfos(i) % faceIdxs(j)))
        payloads(i) % orientatedFaces(j) % face = face
        if (0 < elementInfos(i) % faceIdxs(j)) then
          payloads(i) % orientatedFaces(j) % isOwner = .true.
          payloads(i) % orientatedFaces(j) % outwardNormal = face % ptr % getNormal()

        else
          payloads(i) % orientatedFaces(j) % outwardNormal = -face % ptr % getNormal()

        end if

      end do

      ! Check if we need to construct edges.
      if (.not. allocated(elementInfos(i) % edgeIdxs)) call createEdgeIdxs(elementInfos(i))
      payloads(i) % edges = self % edges % getEdgeBox(elementInfos(i) % edgeIdxs)

      ! Check if we need to construct vertices.
      if (.not. allocated(elementInfos(i) % vertexIdxs)) call createVertexIdxs(elementInfos(i))
      payloads(i) % vertices = self % vertices % getVertexBox(elementInfos(i) % vertexIdxs)

    end do
    call self % elements % init(payloads)
  
  contains
    !!
    !!
    !!
    subroutine createEdgeIdxs(info)
      type(basicElementInfo), intent(inout)                 :: info
      logical(defBool), dimension(self % edges % getSize()) :: isPresent
      integer(shortInt)                                     :: currentSize, idx, k, l, nEdges
      type(faceBox)                                         :: fBox
      type(edgeBox), dimension(:), allocatable              :: faceEdges
      integer(shortInt), dimension(:), allocatable          :: tempIdxs

      ! Initialise isPresent = .false. and allocate info % edgeIdxs to an appropriate initial size.
      isPresent = .false.
      allocate(info % edgeIdxs(6))
      nEdges = 0
      do k = 1, size(info % faceIdxs)
        fBox = self % faces % getFaceBox(abs(info % faceIdxs(k)))
        faceEdges = fBox % ptr % getEdges()
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
      type(faceBox)                                            :: fBox
      type(vertexBox), dimension(:), allocatable               :: faceVertices
      integer(shortInt), dimension(:), allocatable             :: tempIdxs

      ! Initialise isPresent = .false. and allocate info % vertexIdxs to an appropriate initial size.
      isPresent = .false.
      allocate(info % vertexIdxs(4))
      nVertices = 0
      do k = 1, size(info % faceIdxs)
        fBox = self % faces % getFaceBox(abs(info % faceIdxs(k)))
        faceVertices = fBox % ptr % getVertices()
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
    class(unstructuredMesh), intent(inout)             :: self
    type(basicFaceInfo), dimension(:), intent(in)      :: faceInfos
    type(buildFacePayload), dimension(size(faceInfos)) :: payloads
    integer(shortInt)                                  :: i, nFaces

    nFaces = size(faceInfos)
    self % nFaces = nFaces
    do i = 1, nFaces
      payloads(i) % idx = faceInfos(i) % idx
      payloads(i) % parentIdx = faceInfos(i) % parentIdx
      payloads(i) % isBoundary = faceInfos(i) % isBoundary
      payloads(i) % vertices = self % vertices % getVertexBox(faceInfos(i) % vertexIdxs)
      payloads(i) % edges = self % edges % getEdgeBox(faceInfos(i) % edgeIdxs)

    end do
    call self % faces % init(payloads)

  end subroutine initFaceShelf

  !!
  !!
  !!
  subroutine initVertexShelf(self, vertexInfos)
    class(unstructuredMesh), intent(inout)                 :: self
    type(basicVertexInfo), dimension(:), intent(in)        :: vertexInfos
    type(buildVertexPayload), dimension(size(vertexInfos)) :: payloads
    integer(shortInt)                                      :: i, nVertices
    real(defReal), dimension(3)                            :: coords
    real(defReal), dimension(:, :), allocatable            :: allCoords

    nVertices = size(vertexInfos)
    self % nVertices = nVertices
    allocate(allCoords(3, nVertices))

    do i = 1, nVertices
      payloads(i) % idx = vertexInfos(i) % idx
      coords = vertexInfos(i) % coordinates
      payloads(i) % coordinates = coords
      allCoords(:, i) = coords

    end do
    call self % vertices % init(payloads)
    call self % initBoundingBox(allCoords)

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
    call self % acceleration % kill()
    deallocate(self % acceleration)
    call self % elements % kill()
    call self % faces % kill()
    call self % edges % kill()
    call self % vertices % kill()
    deallocate(self % triangulation)

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
    type(elementBox)                    :: element
    
    ! Initialise the numbers of various polyhedra to zero.
    nTetrahedra = 0
    nPentahedra = 0
    nHexahedra = 0
    nOthers = 0
    
    ! Loop over all elements in the mesh.
    do i = 1, self % nElements
      ! Retrieve the number of faces in the current element and increment specific polyhedra
      ! accordingly.
      element = self % elements % getElementBox(i)
      nFaces = size(element % ptr % getOrientatedFaces())
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

  !!
  !!
  !!
  elemental subroutine setEdgesNumber(self, nEdges)
    class(unstructuredMesh), intent(inout) :: self
    integer(shortInt), intent(in)          :: nEdges

    self % nEdges = nEdges

  end subroutine setEdgesNumber

  !!
  !!
  !!
  elemental subroutine setElementsNumber(self, nElements)
    class(unstructuredMesh), intent(inout) :: self
    integer(shortInt), intent(in)          :: nElements

    self % nElements = nElements

  end subroutine setElementsNumber

  !!
  !!
  !!
  elemental subroutine setFacesNumber(self, nFaces)
    class(unstructuredMesh), intent(inout) :: self
    integer(shortInt), intent(in)          :: nFaces

    self % nFaces = nFaces

  end subroutine setFacesNumber

  !!
  !!
  !!
  elemental subroutine setInternalFacesNumber(self, nInternalFaces)
    class(unstructuredMesh), intent(inout) :: self
    integer(shortInt), intent(in)          :: nInternalFaces

    self % nInternalFaces = nInternalFaces

  end subroutine setInternalFacesNumber

  !!
  !!
  !!
  elemental subroutine setVerticesNumber(self, nVertices)
    class(unstructuredMesh), intent(inout) :: self
    integer(shortInt), intent(in)          :: nVertices

    self % nVertices = nVertices

  end subroutine setVerticesNumber

end module unstructuredMesh_inter