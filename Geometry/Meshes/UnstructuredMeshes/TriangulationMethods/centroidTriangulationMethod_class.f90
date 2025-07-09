module centroidTriangulationMethod_class

  use edgeShelf_class,           only : edgeShelf
  use element_class,             only : buildElementInfo, elementBox
  use elementShelf_class,        only : elementShelf
  use face_class,                only : buildFaceInfo, faceBox, orientatedFaceBox
  use faceShelf_class,           only : faceShelf
  use numPrecision
  use triangulationMethod_inter, only : triangulationMethod
  use vertex_class,              only : vertexBox
  use vertexShelf_class,         only : vertexShelf

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
    class(centroidTriangulationMethod), intent(in) :: self
    type(edgeShelf), intent(inout)                 :: edges
    type(faceShelf), intent(inout)                 :: faces
    integer(shortInt)                              :: faceIdx, i, idx, infoIdx, j, minFaceVertexIdx, minFaceVertexIdxLoc, &
                                                      newEdgeIdx, nFaces, newFaceIdx, nTriangles, nVertices
    integer(shortInt), dimension(:), allocatable   :: nTrianglesInFace
    type(faceBox)                                  :: face
    type(buildFaceInfo), dimension(:), allocatable :: faceInfos
    type(vertexBox), dimension(:), allocatable     :: faceVertices

    ! Compute number of triangles to be created first.
    nFaces = faces % getObjectsNumber()
    allocate(nTrianglesInFace(nFaces))
    nTrianglesInFace = 0
    do i = 1, nFaces
      nVertices = size(faces % getFaceVertices(i))
      if (nVertices == 3) cycle
      nTrianglesInFace(i) = nVertices - 2

    end do

    ! Allocate number of new triangles to be generated and populate infos.
    nTriangles = sum(nTrianglesInFace)

    if (nTriangles == 0) return
    allocate(faceInfos(nTriangles))
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
        faceInfos(infoIdx) % idx = newFaceIdx
        faceInfos(infoIdx) % parentIdx = faceIdx
        faceInfos(infoIdx) % isBoundary = face % ptr % getIsBoundary()

        allocate(faceInfos(infoIdx) % vertices(3))
        faceInfos(infoIdx) % vertices = faceVertices([minFaceVertexIdxLoc, &
                                                      mod(minFaceVertexIdxLoc + j - 1, nVertices) + 1, &
                                                      mod(minFaceVertexIdxLoc + j, nVertices) + 1])

        ! Add this child to the face.
        call face % ptr % addChildIdx(newFaceIdx)

      end do
      ! Deactivate face.
      call face % ptr % deactivate()

    end do

    ! Create new faces.
    call self % buildTrianglesFromVertices(faceInfos, edges, faces)

  end subroutine decomposeFaces

  !!
  !!
  !!
  subroutine generateTetrahedra(self, elements, faces, vertices, infos)
    class(centroidTriangulationMethod), intent(in)      :: self
    type(elementShelf), intent(in)                      :: elements
    type(faceShelf), intent(in)                         :: faces
    type(vertexShelf), intent(inout)                    :: vertices
    type(buildElementInfo), dimension(:), intent(inout) :: infos
    integer(shortInt)                                   :: i, infoIdx, j, k, nElements, newElementIdx, newVertexIdx
    type(elementBox)                                    :: element
    type(orientatedFaceBox), dimension(:), allocatable  :: elementOrientatedFaces
    integer(shortInt), dimension(:), allocatable        :: childrenIdxs
    type(faceBox)                                       :: triangle

    ! Initialise nElements and newVertexIdx.
    nElements = elements % getObjectsNumber()
    newVertexIdx = vertices % getObjectsNumber()
    infoIdx = 0
    do i = 1, nElements
      element = elements % getElementBox(i)
      if (size(element % ptr % getVertices()) == 4) cycle

      ! Create a new vertex corresponding to the centroid of the current element.
      newVertexIdx = newVertexIdx + 1
      call vertices % initVertex(newVertexIdx, element % ptr % getCentroid())

      elementOrientatedFaces = element % ptr % getOrientatedFaces()
      do j = 1, size(elementOrientatedFaces)
        childrenIdxs = elementOrientatedFaces(j) % face % ptr % getChildrenIdxs()
        do k = 1, size(childrenIdxs)
          infoIdx = infoIdx + 1
          newElementIdx = nElements + infoIdx
          infos(infoIdx) % idx = newElementIdx
          infos(infoIdx) % localId = element % ptr % getLocalId()
          infos(infoIdx) % parentIdx = element % ptr % getIdx()

          triangle = faces % getFaceBox(childrenIdxs(k))
          
          ! Allocate number of vertices for the current tetrahedron and populate them.
          allocate(infos(infoIdx) % vertices(4))
          infos(infoIdx) % vertices(1:3) = triangle % ptr % getVertices()
          infos(infoIdx) % vertices(4) = vertices % getVertexBox(newVertexIdx)

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
    class(centroidTriangulationMethod), intent(in)     :: self
    type(edgeShelf), intent(inout)                     :: edges
    type(elementShelf), intent(inout)                  :: elements
    type(faceShelf), intent(inout)                     :: faces
    type(vertexShelf), intent(inout)                   :: vertices
    integer(shortInt)                                  :: i, j, nElements, nTetrahedra
    type(elementBox)                                   :: element
    type(orientatedFaceBox), dimension(:), allocatable :: elementOrientatedFaces
    type(buildElementInfo), dimension(:), allocatable  :: elementInfos

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
    allocate(elementInfos(nTetrahedra))

    ! Triangulate faces first.
    call self % decomposeFaces(edges, faces)

    ! Populate vertices within each elementInfo.
    call self % generateTetrahedra(elements, faces, vertices, elementInfos)

    ! Call superclass to build tetrahedra.
    call self % buildTetrahedraFromVertices(elementInfos, edges, elements, faces)

  end subroutine triangulate

end module centroidTriangulationMethod_class