module unstructuredMesh_inter

  use accelerationStructure_inter,   only : accelerationStructure
  use cellZoneShelf_class,           only : cellZoneShelf
  use coord_class,                   only : coord
  use dictionary_class,              only : dictionary
  use edgeShelf_class,               only : edgeShelf
  use element_inter,                 only : elementBox, inclusionTestResult
  use elementShelf_class,            only : elementShelf
  use face_inter,                    only : faceBox
  use faceShelf_class,               only : faceShelf
  use genericProcedures,             only : append, findDifferent, numToChar
  use mesh_inter,                    only : mesh, kill_super => kill
  use numPrecision
  use octreeAcceleration_class,      only : octreeAcceleration
  use patchSingleAcceleration_class, only : patchSingleAcceleration
  use patchMultiAcceleration_class,  only : patchMultiAcceleration
  use universalVariables
  use vertexShelf_class,             only : vertexShelf
  use errors_mod,                    only : fatalError !!!

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
    integer(shortInt), public                 :: nVertices = 0, nFaces = 0, nEdges = 0, &
                                                 nElements = 0, nInternalFaces = 0
    class(accelerationStructure), allocatable :: acceleration
    type(edgeShelf), public                   :: edges
    type(elementShelf), public                :: elements
    type(faceShelf), public                   :: faces
    type(vertexShelf), public                 :: vertices
  contains
    ! Build procedures.
    procedure                       :: computePrimitives
    procedure(importMesh), deferred :: importMesh
    procedure                       :: init
    procedure                       :: kill
    procedure, non_overridable      :: printComposition
    procedure                       :: setEdgeShelf
    procedure                       :: setElementShelf
    procedure                       :: setFaceShelf
    procedure                       :: setVertexShelf
    procedure                       :: split
    procedure                       :: splitElements
    procedure                       :: splitFaces
    ! Runtime procedures.
    procedure                       :: distanceToBoundaryFace
    procedure                       :: distanceToNextFace
    procedure                       :: findHostElement
    procedure                       :: getAllVertexCoordinates
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
    subroutine importMesh(self, folderPath, edges, elements, elementZones, faces, vertices)
      import                                 :: unstructuredMesh, cellZoneShelf, edgeShelf, elementShelf, faceShelf, vertexShelf
      class(unstructuredMesh), intent(inout) :: self
      character(*), intent(in)               :: folderPath
      type(edgeShelf), intent(out)           :: edges
      type(elementShelf), intent(out)        :: elements
      type(cellZoneShelf), intent(out)       :: elementZones
      type(faceShelf), intent(out)           :: faces
      type(vertexShelf), intent(out)         :: vertices

    end subroutine importMesh

  end interface

contains

  !! Subroutine 'computePrimitives'
  !!
  !! Basic description:
  !!   Computes the number of pyramids, triangles and tetrahedra to be created during the mesh
  !!   splitting process.
  !!
  !! Detailed description:
  !!   The number of pyramids is simply given by the sum of the number of faces in each element in
  !!   the original element. The number of triangles is more complex: each pyramid created during
  !!   the splitting process also creates a number of triangles equal to the number of edges (or
  !!   vertices) in the current face. However, since all these triangles are internal they are
  !!   always shared between two pyramids; therefore, the number of triangles to be generated during
  !!   the pyramid creation process is, for a given element, equal to the sum of the number of
  !!   vertices in each of the element's face divided by two. Triangles are also created during the
  !!   splitting of the original mesh's faces: for a given face, the number of triangles to be
  !!   created is simply equal to the number of vertices in the face less two. Lastly, during the
  !!   splitting of pyramids into tetrahedra, additional internal triangles are created, given by
  !!   the number of vertices in a given pyramid's base less three. The number of tetrahedra to be
  !!   generated simply is, for a given face, the number of triangles it is decomposed into.
  !!
  !! Arguments:
  !!   nEdges [out]      -> Number of edges to be generated.
  !!   nTriangles [out]  -> Number of triangles to be generated.
  !!   nTetrahedra [out] -> Number of tetrahedra to be generated.
  !!   nVertices [out]   -> Number of vertices to be generated.
  !!
  elemental subroutine computePrimitives(self, elements, faces, nEdges, nInternalTriangles, nTetrahedra, nTriangles, nVertices)
    class(unstructuredMesh), intent(in)          :: self
    type(elementShelf), intent(in)               :: elements
    type(faceShelf), intent(in)                  :: faces
    integer(shortInt), intent(out)               :: nEdges, nInternalTriangles, nTetrahedra, nTriangles, nVertices
    integer(shortInt)                            :: i, j, nVerticesInElement, nFaces, nVerticesInFace, &
                                                    absFaceIdx
    integer(shortInt), dimension(:), allocatable :: faceIdxs

    ! Initialise nEdges = 0, nInternalTriangles = 0, nTetrahedra = 0, nTriangles = 0 and nVertices = 0.
    nEdges = 0
    nInternalTriangles = 0
    nTetrahedra = 0
    nTriangles = 0
    nVertices = 0

    ! Loop through all elements.
    do i = 1, self % nElements
      ! Retrieve the number of vertices and indices of the faces in the current element.
      nVerticesInElement = size(elements % getElementVertexIdxs(i))
      faceIdxs = elements % getElementFaceIdxs(i)
      
      ! Check if the current element is already a tetrahedron. If yes, increment nTetrahedra by 1
      ! and nTriangles by the number of triangles owned by the tetrahedron then cycle.
      if (nVerticesInElement == 4) then
        nTetrahedra = nTetrahedra + 1
        nTriangles = nTriangles + count(faceIdxs > 0)

        do j = 1, 4
          if (faceIdxs(j) > 0) then
            if (.not. faces % getFaceIsBoundary(faceIdxs(j))) nInternalTriangles = nInternalTriangles + 1

          end if

        end do
        cycle
      
      end if

      ! If the current element is not a tetrahedron it will be split from its centroid so we need
      ! to add the current element's centroid to the list of vertices.
      nVertices = nVertices + 1

      ! Increment nEdges.
      nEdges = nEdges + nVerticesInElement

      ! Compute the number of faces in the current element.
      nFaces = size(faceIdxs)
      
      ! Initialise nVertices and loop through all faces.
      nVerticesInElement = 0
      do j = 1, nFaces
        absFaceIdx = abs(faceIdxs(j))
        ! Retrieve the number of vertices in the current face and increase the total 
        ! number of vertices by the number of vertices in the current face.
        nVerticesInFace = size(faces % getFaceVertexIdxs(absFaceIdx))
        nVerticesInElement = nVerticesInElement + nVerticesInFace
        
        ! Increase the number of triangles corresponding to new internal faces by nVerticesInFace - 3.
        nTriangles = nTriangles + nVerticesInFace - 3
        nInternalTriangles = nInternalTriangles + nVerticesInFace - 3

        ! If the element owns the current face, increase the number of triangles by nVerticesInFace - 2.
        if (faceIdxs(j) > 0) then
          nEdges = nEdges + nVerticesInFace - 3
          nTriangles = nTriangles + nVerticesInFace - 2
          if (.not. faces % getFaceIsBoundary(faceIdxs(j))) nInternalTriangles = nInternalTriangles + nVerticesInFace - 2

        end if
        
        ! There will be as many tetrahedra as the number of triangles in each face, which is given
        ! by nVerticesInFace - 2.
        nTetrahedra = nTetrahedra + nVerticesInFace - 2
      
      end do
      
      ! The number of pyramids' faces is given by half the total number of vertices.
      nTriangles = nTriangles + nVerticesInElement / 2
      nInternalTriangles = nInternalTriangles + nVerticesInElement / 2

    end do

  end subroutine computePrimitives

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
      call self % faces % computeFaceIntersection(i, coords, self % vertices, update)
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
      testResult = self % elements % isPointInside(elementIdxs(1), coords % getEndPosition(), self % faces)
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
    elementIdx = coords % getElementIdx()
    rEnd = coords % getEndPosition()
    potentialFaces = self % elements % computePotentialFaceIdxs(elementIdx, rEnd, self % faces)

    ! If no potential intersections are detected return early.
    if (size(potentialFaces) == 0) return
    
    ! If reached here, compute which face is actually intersected and update d.
    r = coords % getPosition()
    call self % elements % computeFaceIntersection(elementIdx, r, rEnd, potentialFaces, self % faces, &
                                                   intersectedFaceIdx, lambda)
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
    call coords % setLocalId(self % findElementZoneIdx(self % elements % getElementParentIdx(coords % getElementIdx())))

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
    integer(shortInt)           :: coordPatch, coordBrute, coordOctree !!!
    
    ! Initialise parentIdx = 0. Retrieve the mesh's bounding box. If the particle is outside the bounding box we can return early.
    call coords % setElementIdx(0)
    call coords % setParentElementIdx(0)
    if (allocated(self % acceleration)) then
      call self % acceleration % findHostElement(self % vertices, self % edges, self % faces, self % elements, coords)
      
    !!!
    coordPatch = coords % getElementIdx()
    !coordPatch = coords % getParentElementIdx()
    end if
    !!!

    !else !!!
      ! Perform brute-force search.
      searchLoop: do
        do i = 1, self % nElements
          testResult = self % elements % isPointInside(i, coords % getPositionToNudge(), self % faces)
          if (testResult % status == INSIDE_ELEMENT) then
            call coords % setElementIdx(i)
            call coords % setParentElementIdx(self % elements % getElementParentIdx(i))

            !!!
            coordBrute = coords % getElementIdx()
            !coordBrute = coords % getParentElementIdx()
            ! if (coordBrute /= 0) then
            ! print*, coordBrute
            ! end if
            ! print*, "INSIDE"
            ! print*, coordBrute, coordPatch
            if (coordBrute /= coordPatch) then
              print*, "Brute:", coordBrute
              print*, "Patch:", coordPatch
              call fatalError("INSIDE_ELEMENT", "Element indices not matching between the two methods")
            end if
            !!!

            return

          elseif (testResult % status == ON_BOUNDARY_ELEMENT) then
            ! If coordinates are on the element boundary (very rare), we need to push them off.
            do while (testResult % status == ON_BOUNDARY_ELEMENT)
              call self % elements % pushFromElementBoundary(i, self % faces, coords)

              ! Perform containment test again.
              testResult = self % elements % isPointInside(i, coords % getPositionToNudge(), self % faces)

            end do

            ! Now the coordinates are not on the boundary of the element anymore.
            if (testResult % status == INSIDE_ELEMENT) then
              ! If coordinates are now well inside the element, we have found our element.
              call coords % setElementIdx(i)
              call coords % setParentElementIdx(self % elements % getElementParentIdx(i))
              
              !!!
              coordBrute = coords % getElementIdx()
              !coordBrute = coords % getParentElementIdx()
              ! if (coordBrute /= 0) then
              ! print*, coordBrute
              ! end if
              ! print*, "INSIDE Element"
              ! print*, coordBrute, coordPatch
              if (coordBrute /= coordPatch) then
                print*, "Brute:", coordBrute
                print*, "Patch:", coordPatch
                call fatalError("ON_BOUNDARY_ELEMENT", "Element indices not matching between the two methods")
              end if
              !!!
              
              
              
              
              return

            elseif (testResult % status == OUTSIDE_ELEMENT) then
              ! If the nudge has resulted in an overshoot, we cycle searchLoop and begin the entire process again.
              cycle searchLoop

            end if

          end if

        end do


        !!!
        coordBrute = coords % getElementIdx()
        !coordBrute = coords % getParentElementIdx()
        ! if (coordBrute /= 0) then
        ! print*, coordBrute
        ! end if
        ! print*, "OUTSIDE Element"
        ! print*, coordBrute, coordPatch
        if (coordBrute /= coordPatch) then
          print*, "Brute:", coordBrute
          print*, "Patch:", coordPatch
          call fatalError("OUTSIDE_ELEMENT", "Element indices not matching between the two methods")
        end if
        !!!

        return

      end do searchLoop


    !end if !!!

  end subroutine findHostElement

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
    real(defReal), dimension(3, self % nVertices) :: coords

    coords = self % vertices % getAllCoordinates()

  end function getAllVertexCoordinates

  !!
  !!
  !!
  subroutine init(self, folderPath, dict)
    class(unstructuredMesh), intent(inout) :: self
    character(*), intent(in)               :: folderPath
    class(dictionary), intent(in)          :: dict
    type(edgeShelf)                        :: edges, newEdges
    type(elementShelf)                     :: elements, newElements
    type(cellZoneShelf)                    :: elementZones
    type(faceShelf)                        :: faces, newFaces
    type(vertexShelf)                      :: newVertices, vertices
    logical(defBool)                       :: triangulate
    character(nameLen)                     :: acceleration

    ! Set up base components.
    call self % setupBase(dict)
    
    ! Import mesh from files.
    call self % importMesh(folderPath, edges, elements, elementZones, faces, vertices)

    ! Check if triangulation was requested.
    call dict % getOrDefault(triangulate, 'triangulate', .false.)
    if (triangulate) then
      call self % split(edges, elements, faces, vertices, newEdges, newElements, newFaces, newVertices)
      call self % setEdgeShelf(newEdges)
      call self % setElementShelf(newElements)
      call self % setFaceShelf(newFaces)
      call self % setVertexShelf(newVertices)

    else
      call self % setEdgeShelf(edges)
      call self % setElementShelf(elements)
      call self % setFaceShelf(faces)
      call self % setVertexShelf(vertices)

    end if

    ! Set elements zones.
    call self % setElementZones(elementZones)

    ! Check if acceleration structure was required by user and initialise it if applicable.
    call dict % getOrDefault(acceleration, 'accelerationMethod', 'none')
    if (acceleration /= 'none') then
      if (acceleration == 'octree') allocate(octreeAcceleration :: self % acceleration)
      if (acceleration == 'patchSingle') allocate(patchSingleAcceleration :: self % acceleration)
      if (acceleration == 'patchMulti') allocate(patchMultiAcceleration :: self % acceleration)
      call self % acceleration % init(self % vertices, self % edges, self % faces, self % elements)
    end if

  end subroutine init

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
    call self % edges % kill()
    call self % elements % kill()
    call self % faces % kill()
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
    
    ! Retrieve sizes of the original shelves.
    nEdges = self % nEdges
    nVertices = self % nVertices
    
    ! Compute the number of edges, pyramids, triangles and tetrahedra to be created and
    ! allocate memory to the corresponding structures.
    call self % computePrimitives(elements, faces, nNewEdges, nInternalTriangles, nTetrahedra, nTriangles, nNewVertices)
    
    ! Allocate memory in the new shelves.
    call newEdges % allocateShelf(nEdges + nNewEdges)
    call newElements % allocateShelf(nTetrahedra)
    call newFaces % allocateShelf(nTriangles)
    call newVertices % allocateShelf(nVertices + nNewVertices)

    ! Copy original edges and vertices into the new shelves.
    do i = 1, nEdges
      call newEdges % initEdge(i, edges % getEdgeVertexIdxs(i))

    end do

    call newVertices % setExtremalCoordinates(vertices % getExtremalCoordinates())
    call newVertices % setOffset(vertices % getOffset())
    do i = 1, nVertices
      call newVertices % initVertex(i, vertices % getVertexCoordinates(i))
      edgeIdxs = vertices % getVertexEdgeIdxs(i)

      do j = 1, size(edgeIdxs)
        call newVertices % addEdgeIdxToVertex(i, edgeIdxs(j))

      end do

    end do

    ! Initialise new triangles and tetrahedra to be generated.
    allocate(triangles(nTriangles))
    allocate(tetrahedra(nTetrahedra))
    
    ! Initialise lastVertexIdx, lastPyramidIdx, lastTetrahedronIdx and lastTriangleIdx then
    ! split all elements into pyramids and all pyramids into tetrahedra.
    lastEdgeIdx = nEdges
    lastElementIdx = 0
    lastFaceIdx = 0
    lastVertexIdx = nVertices
    call self % splitFaces(faces, newEdges, newFaces, newVertices, lastEdgeIdx, lastFaceIdx, triangles)
    call self % splitElements(elements, faces, lastEdgeIdx, lastElementIdx, lastFaceIdx, lastVertexIdx, &
                              newEdges, newElements, newFaces, newVertices, tetrahedra, triangles)

    ! Update the number of edges, faces, elements and vertices in the mesh.
    self % nEdges = nEdges + nNewEdges
    self % nElements = nTetrahedra
    self % nFaces = nTriangles
    self % nInternalFaces = nInternalTriangles
    self % nVertices = nVertices + nNewVertices

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
                                        
    ! Loop through all original elements and split them into tetrahedra.
    do i = 1, self % nElements
      ! Initialise initialElementIdx then split the current element.
      initialElementIdx = lastNewElementIdx + 1
      call elements % splitElement(i, faces, lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx, &
                                   newEdges, newFaces, newVertices, tetrahedra, triangles)

      ! Set all new tetrahedra.
      do j = initialElementIdx, lastNewElementIdx
        call newElements % addElement(j, tetrahedra(j))

      end do
      
    end do

  end subroutine splitElements

  subroutine splitFaces(self, faces, newEdges, newFaces, newVertices, lastNewEdgeIdx, lastNewFaceIdx, triangles)
    class(unstructuredMesh), intent(inout)     :: self
    type(faceShelf), intent(inout)             :: faces, newFaces
    type(edgeShelf), intent(inout)             :: newEdges
    type(vertexShelf), intent(inout)           :: newVertices
    integer(shortInt), intent(inout)           :: lastNewEdgeIdx, lastNewFaceIdx
    type(faceBox), dimension(:), intent(inout) :: triangles
    integer(shortInt)                          :: i, initialFaceIdx, j

    ! Loop through all original faces in the mesh and split them into triangles.
    do i = 1, self % nFaces
      ! Initialise initialFaceIdx then split the current face.
      initialFaceIdx = lastNewFaceIdx + 1
      call faces % splitFace(i, newEdges, newVertices, lastNewEdgeIdx, lastNewFaceIdx, triangles)

      ! Set all new triangles.
      do j = initialFaceIdx, lastNewFaceIdx
        call newFaces % addFace(j, triangles(j))

      end do

    end do

  end subroutine splitFaces

end module unstructuredMesh_inter