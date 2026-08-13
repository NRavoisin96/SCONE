module centroidTriangulationMethod_class

  use element_class,                only : buildElementPayload, elementBox
  use face_class,                   only : buildFacePayload, faceBox, orientatedFaceBox
  use numPrecision
  use topologicalObjectShelf_class, only : topologicalObjectShelf
  use triangulationMethod_inter,    only : triangulationMethod
  use vertex_class,                 only : buildVertexPayload, vertexBox

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(triangulationMethod) :: centroidTriangulationMethod
    private
  contains
    procedure :: decomposeFaces
    procedure :: generateTetrahedra
    procedure :: triangulate
  end type centroidTriangulationMethod

contains
  !!
  !!
  !!
  subroutine decomposeFaces(self, edges, faces)
    class(centroidTriangulationMethod), intent(in)    :: self
    type(topologicalObjectShelf), intent(inout)       :: edges, faces
    integer(shortInt)                                 :: faceIdx, i, idx, infoIdx, j, minFaceVertexIdx, minFaceVertexIdxLoc, &
                                                         newEdgeIdx, nFaces, newFaceIdx, nTriangles, nVertices
    integer(shortInt), dimension(:), allocatable      :: nTrianglesInFace
    type(faceBox)                                     :: face
    type(buildFacePayload), dimension(:), allocatable :: facePayloads
    type(vertexBox), dimension(:), allocatable        :: faceVertices

    ! Compute number of triangles to be created first.
    nFaces = faces % getObjectsNumber()
    allocate(nTrianglesInFace(nFaces))
    nTrianglesInFace = 0
    do i = 1, nFaces
      face = faces % getFaceBox(i)
      nVertices = size(face % ptr % getVertices())
      if (nVertices == 3) cycle
      nTrianglesInFace(i) = nVertices - 2

    end do

    ! Allocate number of new triangles to be generated and populate infos.
    nTriangles = sum(nTrianglesInFace)

    if (nTriangles == 0) return
    allocate(facePayloads(nTriangles))
    infoIdx = 0
    newEdgeIdx = edges % getObjectsNumber()
    do i = 1, nFaces
      if (nTrianglesInFace(i) == 0) cycle
      
      face = faces % getFaceBox(i)
      faceIdx = face % ptr % getIdx()
      faceVertices = face % ptr % getVertices()
      nVertices = size(faceVertices)
      
      ! Find the vertex in the face with the smallest index.
      minFaceVertexIdx = faceVertices(1) % ptr % getIdx()
      minFaceVertexIdxLoc = 1
      do j = 2, nVertices
        idx = faceVertices(j) % ptr % getIdx()
        if (idx < minFaceVertexIdx) then
          minFaceVertexIdx = idx
          minFaceVertexIdxLoc = j

        end if

      end do

      do j = 1, nTrianglesInFace(i)
        infoIdx = infoIdx + 1
        newFaceIdx = nFaces + infoIdx
        facePayloads(infoIdx) % idx = newFaceIdx
        facePayloads(infoIdx) % parentIdx = faceIdx
        facePayloads(infoIdx) % isBoundary = face % ptr % getIsBoundary()
        facePayloads(infoIdx) % boundaryConditions = face % ptr % getBoundaryConditions()
        facePayloads(infoIdx) % boundaryValues = face % ptr % getBoundaryValues()

        allocate(facePayloads(infoIdx) % vertices(3))
        facePayloads(infoIdx) % vertices = faceVertices([minFaceVertexIdxLoc, &
                                                         mod(minFaceVertexIdxLoc + j - 1, nVertices) + 1, &
                                                         mod(minFaceVertexIdxLoc + j, nVertices) + 1])

        ! Add this child to the face.
        call face % ptr % addChildIdx(newFaceIdx)

      end do
      ! Deactivate face.
      call face % ptr % deactivate()

    end do

    ! Create new faces.
    call self % buildTrianglesFromVertices(facePayloads, edges, faces)

  end subroutine decomposeFaces

  !!
  !!
  !!
  subroutine generateTetrahedra(self, elements, faces, vertices, payloads)
    class(centroidTriangulationMethod), intent(in)         :: self
    type(topologicalObjectShelf), intent(in)               :: elements, faces
    type(topologicalObjectShelf), intent(inout)            :: vertices
    type(buildElementPayload), dimension(:), intent(inout) :: payloads
    integer(shortInt)                                      :: i, infoIdx, j, k, nElements, newElementIdx, newVertexIdx
    type(elementBox)                                       :: element
    type(buildVertexPayload)                               :: vertexPayload
    type(orientatedFaceBox), dimension(:), allocatable     :: elementOrientatedFaces
    integer(shortInt), dimension(:), allocatable           :: childrenIdxs
    type(faceBox)                                          :: triangle

    ! Initialise nElements and newVertexIdx.
    nElements = elements % getObjectsNumber()
    newVertexIdx = vertices % getObjectsNumber()
    infoIdx = 0
    do i = 1, nElements
      element = elements % getElementBox(i)
      if (size(element % ptr % getVertices()) == 4) cycle

      ! Create a new vertex corresponding to the centroid of the current element.
      newVertexIdx = newVertexIdx + 1
      vertexPayload % idx = newVertexIdx
      vertexPayload % coordinates = element % ptr % getCentroid()
      vertexPayload % ratintCoordinates = element % ptr % getRationalCentroid()
      call vertices % initObject(vertexPayload)

      elementOrientatedFaces = element % ptr % getOrientatedFaces()
      do j = 1, size(elementOrientatedFaces)
        childrenIdxs = elementOrientatedFaces(j) % face % ptr % getChildrenIdxs()
        do k = 1, size(childrenIdxs)
          infoIdx = infoIdx + 1
          newElementIdx = nElements + infoIdx
          payloads(infoIdx) % idx = newElementIdx
          payloads(infoIdx) % localId = element % ptr % getLocalId()
          payloads(infoIdx) % parentIdx = element % ptr % getIdx()

          triangle = faces % getFaceBox(childrenIdxs(k))
          
          ! Allocate number of vertices for the current tetrahedron and populate them.
          allocate(payloads(infoIdx) % vertices(4))
          payloads(infoIdx) % vertices(1:3) = triangle % ptr % getVertices()
          payloads(infoIdx) % vertices(4) = vertices % getVertexBox(newVertexIdx)

        end do

      end do
      ! Deactivate current element.
      call element % ptr % deactivate()

    end do

  end subroutine generateTetrahedra

  !!
  !!
  !!
  subroutine triangulate(self, edges, elements, faces, vertices)
    class(centroidTriangulationMethod), intent(in)       :: self
    type(topologicalObjectShelf), intent(inout)          :: edges, elements, faces, vertices
    integer(shortInt)                                    :: i, j, nElements, nTetrahedra
    type(elementBox)                                     :: element
    type(orientatedFaceBox), dimension(:), allocatable   :: elementOrientatedFaces
    type(buildElementPayload), dimension(:), allocatable :: tetrahedraPayloads

    ! Count the number of tetrahedra to be genetated.
    nTetrahedra = 0
    nElements = elements % getObjectsNumber()
    do i = 1, nElements
      element = elements % getElementBox(i)
      if (size(element % ptr % getVertices()) == 4) cycle
      
      elementOrientatedFaces = element % ptr % getOrientatedFaces()
      do j = 1, size(elementOrientatedFaces)
        nTetrahedra = nTetrahedra + size(elementOrientatedFaces(j) % face % ptr % getVertices()) - 2

      end do

    end do

    ! Return early if no tetrahedra need to be generated. Allocate elementInfos otherwise.
    if (nTetrahedra == 0) return
    allocate(tetrahedraPayloads(nTetrahedra))

    ! Triangulate faces first.
    call self % decomposeFaces(edges, faces)

    ! Populate vertices within each elementInfo.
    call self % generateTetrahedra(elements, faces, vertices, tetrahedraPayloads)

    ! Call superclass to build tetrahedra.
    call self % buildTetrahedraFromVertices(tetrahedraPayloads, edges, elements, faces)

  end subroutine triangulate

end module centroidTriangulationMethod_class