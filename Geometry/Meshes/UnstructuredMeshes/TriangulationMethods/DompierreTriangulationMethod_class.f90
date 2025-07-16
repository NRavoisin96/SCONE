module DompierreTriangulationMethod_class

  use edge_class,                   only : edge, edgeBox
  use element_class,                only : buildElementPayload, elementBox
  use face_class,                   only : buildFacePayload, faceBox, orientatedFaceBox
  use genericProcedures,            only : areEqual, fatalError, numToChar, quickSort
  use numPrecision
  use topologicalObject_inter,      only : topologicalObjectBox
  use topologicalObjectShelf_class, only : topologicalObjectShelf
  use triangulationMethod_inter,    only : triangulationMethod
  use universalVariables,           only : NOT_PRESENT
  use vertex_class,                 only : vertexBox

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(triangulationMethod) :: DompierreTriangulationMethod
    private
  contains
    procedure          :: decomposeFaces
    procedure, private :: generateTetrahedraFromHexahedron
    procedure          :: triangulate
  end type DompierreTriangulationMethod

contains
  !!
  !!
  !!
  subroutine decomposeFaces(self, edges, faces)
    class(DompierreTriangulationMethod), intent(in)   :: self
    type(topologicalObjectShelf), intent(inout)       :: edges, faces
    integer(shortInt)                                 :: i, infoIdx, j, &
                                                         newEdgeIdx, newFaceIdx, nFaces, nVertices, nTriangles
    type(faceBox)                                     :: face
    type(vertexBox), dimension(:), allocatable        :: faceVertices
    type(buildFacePayload), dimension(:), allocatable :: trianglePayloads
    integer(shortInt), dimension(4)                   :: vertexIdxs
    real(defReal), dimension(3)                       :: diff_1, diff_2
    real(defReal)                                     :: dSquared_1, dSquared_2
    logical(defBool)                                  :: splitAlong13
    character(*), parameter                           :: here = 'decomposeFaces (DompierreTriangulationMethod_class.f90)'

    ! Do a first pass and count the number of triangles to be generated.
    nTriangles = 0
    nFaces = faces % getObjectsNumber()
    do i = 1, nFaces
      face = faces % getFaceBox(i)
      faceVertices = face % ptr % getVertices()
      nVertices = size(faceVertices)
      select case(nVertices)
        case(3)
          ! Do nothing.

        case(4)
          ! Add two new triangles.
          nTriangles = nTriangles + 2

        case default
          ! Call fatalError.
          call fatalError(here, 'Invalid number of vertices: '//numToChar(nVertices)//' for face: '&
                          //numToChar(face % ptr % getIdx())//'.')

      end select

    end do

    ! Return if no new triangles need to be generated.
    if (nTriangles == 0) return
    allocate(trianglePayloads(nTriangles))

    infoIdx = 0
    newEdgeIdx = edges % getObjectsNumber()
    do i = 1, nFaces
      face = faces % getFaceBox(i)
      faceVertices = face % ptr % getVertices()
      nVertices = size(faceVertices)

      if (nVertices == 4) then
        ! Decompose the quadrilateral into two triangles. Decompose along the shortest diagonal. If diagonals have
        ! the same length, use a tie-breaker based on the smallest index sum.
        do j = 1, 4
          vertexIdxs(j) = faceVertices(j) % ptr % getIdx()

        end do
        diff_1 = faceVertices(3) % ptr % getCoordinates() - faceVertices(1) % ptr % getCoordinates()
        diff_2 = faceVertices(4) % ptr % getCoordinates() - faceVertices(2) % ptr % getCoordinates()
        dSquared_1 = dot_product(diff_1, diff_1)
        dSquared_2 = dot_product(diff_2, diff_2)
        splitAlong13 = dSquared_1 < dSquared_2 .or. (areEqual(dSquared_1, dSquared_2) .and. &
                       sum(vertexIdxs([1, 3])) <= sum(vertexIdxs([2, 4])))

        do j = 1, 2
          infoIdx = infoIdx + 1
          newFaceIdx = nFaces + infoIdx
          trianglePayloads(infoIdx) % idx = newFaceIdx
          trianglePayloads(infoIdx) % parentIdx = face % ptr % getIdx()
          trianglePayloads(infoIdx) % isBoundary = face % ptr % getIsBoundary()

          allocate(trianglePayloads(infoIdx) % vertices(3))
          if (splitAlong13) then
            trianglePayloads(infoIdx) % vertices = faceVertices([1, j + 1, j + 2])

          else
            trianglePayloads(infoIdx) % vertices = faceVertices([2, j + 2, merge(1, j + 3, j == 2)])

          end if

          ! Add this child to the face.
          call face % ptr % addChildIdx(newFaceIdx)

        end do
        ! Deactivate face.
        call face % ptr % deactivate()

      end if

    end do
    
    ! Create new faces.
    call self % buildTrianglesFromVertices(trianglePayloads, edges, faces)

  end subroutine decomposeFaces

  !!
  !!
  !!
  subroutine generateTetrahedraFromHexahedron(self, hexahedron, vertices, edges, infoIdx, newElementIdx, tetrahedraPayloads)
    class(DompierreTriangulationMethod), intent(in)        :: self
    type(elementBox), intent(in)                           :: hexahedron
    type(topologicalObjectShelf), intent(in)               :: vertices
    type(topologicalObjectShelf), intent(inout)            :: edges
    integer(shortInt), intent(inout)                       :: infoIdx, newElementIdx
    type(buildElementPayload), dimension(:), intent(inout) :: tetrahedraPayloads
    type(vertexBox), dimension(8)                          :: hexahedronVertices, sortedVertices
    integer(shortInt)                                      :: i, idx, j, nEdges, parity
    integer(shortInt), dimension(3)                        :: edgeVertexIdxs
    integer(shortInt), dimension(8)                        :: vertexIdxs
    type(topologicalObjectBox), dimension(:), allocatable  :: vertexEdges
    integer(shortInt), dimension(:), allocatable           :: edgeIdxs
    type(edgeBox), dimension(12)                           :: hexahedronEdges
    type(vertexBox), dimension(2)                          :: edgeVertices
    type(orientatedFaceBox), dimension(6)                  :: hexahedronOrientatedFaces
    type(vertexBox), dimension(4)                          :: orientatedFaceVertices
    integer(shortInt), dimension(4)                        :: orientatedFaceVertexIdxs
    logical(defBool)                                       :: isEven
    character(*), parameter :: here = 'generateTetrahedraFromHexahedron (DompierreTriangulationMethod_class.f90)'

    hexahedronVertices = hexahedron % ptr % getVertices()
    do i = 1, 8
      vertexIdxs(i) = hexahedronVertices(i) % ptr % getIdx() 

    end do

    ! First sorted vertex is the vertex with smallest global index.
    sortedVertices(1) = hexahedronVertices(minloc(vertexIdxs, 1))
    vertexEdges = sortedVertices(1) % ptr % getSharingEdges()
    nEdges = size(vertexEdges)
    allocate(edgeIdxs(nEdges))
    do i = 1, nEdges
      ! Downcast to correct type.
      select type(ptr => vertexEdges(i) % ptr)
        type is(edge)
          edgeIdxs(i) = ptr % getIdx()

        class default
          call fatalError(here, 'Edge with index: '//numToChar(ptr % getIdx())//' is not an edge.')

      end select

    end do

    ! Loop through all the edges containing the first vertex and find the second, fourth, and fifth vertices.
    hexahedronEdges = hexahedron % ptr % getEdges()
    idx = 0
    do i = 1, 12
      if (any(edgeIdxs == hexahedronEdges(i) % ptr % getIdx())) then
        idx = idx + 1
        edgeVertices = hexahedronEdges(i) % ptr % getVertices()
        if (associated(sortedVertices(1) % ptr, edgeVertices(1) % ptr)) then
          edgeVertexIdxs(idx) = edgeVertices(2) % ptr % getIdx()

        else
          edgeVertexIdxs(idx) = edgeVertices(1) % ptr % getIdx()

        end if

      end if

    end do

    ! Sort edgeVertexIdxs and assign the second, fourth, and fifth vertices.
    call quickSort(edgeVertexIdxs)
    sortedVertices([2, 4, 5]) = vertices % getVertexBox(edgeVertexIdxs)

    ! Now find the third vertex. It is the vertex which shares the face with vertices 1, 2, and 4.
    hexahedronOrientatedFaces = hexahedron % ptr % getOrientatedFaces()
    do i = 1, 6
      orientatedFaceVertices = hexahedronOrientatedFaces(i) % face % ptr % getVertices()
      do j = 1, 4
        orientatedFaceVertexIdxs(j) = orientatedFaceVertices(j) % ptr % getIdx()

      end do
      if (any(orientatedFaceVertexIdxs == sortedVertices(1) % ptr % getIdx()) .and. &
          any(orientatedFaceVertexIdxs == sortedVertices(2) % ptr % getIdx()) .and. &
          any(orientatedFaceVertexIdxs == sortedVertices(4) % ptr % getIdx())) then

        ! We have found our face. Now assign vertex 3.
        do j = 1, 4
          if (associated(orientatedFaceVertices(j) % ptr, sortedVertices(1) % ptr) .or. &
              associated(orientatedFaceVertices(j) % ptr, sortedVertices(2) % ptr) .or. &
              associated(orientatedFaceVertices(j) % ptr, sortedVertices(4) % ptr)) cycle
          sortedVertices(3) = orientatedFaceVertices(j)

        end do

      end if

    end do

    ! Now assign vertices 6, 7, and 8. Vertex 6 is connected to vertex 2.
    do i = 2, 4
      if (allocated(edgeIdxs)) deallocate(edgeIdxs)
      vertexEdges = sortedVertices(i) % ptr % getSharingEdges()
      nEdges = size(vertexEdges)
      allocate(edgeIdxs(nEdges))
      do j = 1, nEdges
        ! Downcast to correct type.
        select type(ptr => vertexEdges(j) % ptr)
          type is(edge)
            edgeIdxs(j) = ptr % getIdx()

          class default
            call fatalError(here, 'Edge with index: '//numToChar(ptr % getIdx())//' is not an edge.')

        end select

      end do

      do j = 1, 12
        if (any(edgeIdxs == hexahedronEdges(j) % ptr % getIdx())) then
          edgeVertices = hexahedronEdges(j) % ptr % getVertices()
          if (associated(sortedVertices(i) % ptr, edgeVertices(1) % ptr)) then
            if (associated(sortedVertices(i - 1) % ptr, edgeVertices(2) % ptr) .or. &
                associated(sortedVertices(merge(1, i + 1, i == 4)) % ptr, edgeVertices(2) % ptr)) cycle

              ! We have found our vertex.
              sortedVertices(i + 4) = edgeVertices(2)

          else
            if (associated(sortedVertices(i - 1) % ptr, edgeVertices(1) % ptr) .or. &
                associated(sortedVertices(merge(1, i + 1, i == 4)) % ptr, edgeVertices(1) % ptr)) cycle

              ! We have found our vertex.
              sortedVertices(i + 4) = edgeVertices(1)

          end if

        end if

      end do

    end do

    ! Now that we have sorted the vertices in the correct order, simply compute the parity of the
    ! hexahedron.
    parity = 0
    if (edges % getObjectIdxOrDefault(sortedVertices([1, 3]), NOT_PRESENT) == NOT_PRESENT) parity = parity + 1
    if (edges % getObjectIdxOrDefault(sortedVertices([1, 6]), NOT_PRESENT) /= NOT_PRESENT) parity = parity + 1
    if (edges % getObjectIdxOrDefault(sortedVertices([1, 8]), NOT_PRESENT) /= NOT_PRESENT) parity = parity + 1
    isEven = mod(parity, 2) == 0
    do i = 1, 6
      infoIdx = infoIdx + 1
      newElementIdx = newElementIdx + 1
      tetrahedraPayloads(infoIdx) % idx = newElementIdx
      tetrahedraPayloads(infoIdx) % parentIdx = hexahedron % ptr % getIdx()
      tetrahedraPayloads(infoIdx) % localId = hexahedron % ptr % getLocalId()
      allocate(tetrahedraPayloads(infoIdx) % vertices(4))
      select case(i)
        case(1)
          if (isEven) then
            tetrahedraPayloads(infoIdx) % vertices = sortedVertices([1, 2, 5, 8])

          else
            tetrahedraPayloads(infoIdx) % vertices = sortedVertices([1, 2, 3, 7])

          end if

        case(2)
          if (isEven) then
            tetrahedraPayloads(infoIdx) % vertices = sortedVertices([1, 2, 4, 8])

          else
            tetrahedraPayloads(infoIdx) % vertices = sortedVertices([1, 3, 4, 7])

          end if

        case(3)
          if (isEven) then
            tetrahedraPayloads(infoIdx) % vertices = sortedVertices([2, 3, 4, 8])

          else
            tetrahedraPayloads(infoIdx) % vertices = sortedVertices([1, 4, 8, 7])

          end if

        case(4)
          if (isEven) then
            tetrahedraPayloads(infoIdx) % vertices = sortedVertices([2, 5, 6, 8])

          else
            tetrahedraPayloads(infoIdx) % vertices = sortedVertices([1, 8, 5, 7])

          end if

        case(5)
          if (isEven) then
            tetrahedraPayloads(infoIdx) % vertices = sortedVertices([2, 3, 6, 8])

          else
            tetrahedraPayloads(infoIdx) % vertices = sortedVertices([1, 5, 6, 7])

          end if

        case(6)
          if (isEven) then
            tetrahedraPayloads(infoIdx) % vertices = sortedVertices([3, 6, 7, 8])

          else
            tetrahedraPayloads(infoIdx) % vertices = sortedVertices([1, 6, 2, 7])

          end if

      end select

    end do

  end subroutine generateTetrahedraFromHexahedron

  !!
  !!
  !!
  subroutine triangulate(self, edges, elements, faces, vertices)
    class(DompierreTriangulationMethod), intent(in)      :: self
    type(topologicalObjectShelf), intent(inout)          :: edges, elements, faces, vertices
    integer(shortInt)                                    :: i, infoIdx, nElements, newElementIdx, nFaces, nTetrahedra, nVertices
    type(elementBox)                                     :: element
    type(orientatedFaceBox), dimension(:), allocatable   :: elementOrientatedFaces
    type(vertexBox), dimension(:), allocatable           :: elementVertices
    type(buildElementPayload), dimension(:), allocatable :: tetrahedraPayloads
    character(*), parameter                              :: here = 'triangulate (DompierreTriangulationMethod_class.f90)'

    ! Do a first pass to count the number of tetrahedra to be generated.
    nElements = elements % getObjectsNumber()
    nTetrahedra = 0
    do i = 1, nElements
      element = elements % getElementBox(i)
      elementOrientatedFaces = element % ptr % getOrientatedFaces()
      nFaces = size(elementOrientatedFaces)
      elementVertices = element % ptr % getVertices()
      nVertices = size(elementVertices)

      select case(nVertices)
        case(4)
          ! Do nothing.

        case(5)
          if (nFaces /= 5) call fatalError(here, 'Invalid element type.')
          nTetrahedra = nTetrahedra + 2

        case(6)
          if (nFaces /= 5) call fatalError(here, 'Invalid element type.')
          nTetrahedra = nTetrahedra + 3

        case(8)
          if (nFaces /= 6) call fatalError(here, 'Invalid element type.')
          nTetrahedra = nTetrahedra + 6

        case default
          call fatalError(here, 'Invalid element type.')

      end select

    end do

    ! If no new tetrahedra need to be generated return early.
    if (nTetrahedra == 0) return
    allocate(tetrahedraPayloads(nTetrahedra))
    
    ! Decompose faces first.
    call self % decomposeFaces(edges, faces)

    ! Loop through all elements and generate tetrahedra.
    infoIdx = 0
    newElementIdx = elements % getObjectsNumber()
    do i = 1, nElements
      element = elements % getElementBox(i)
      elementVertices = element % ptr % getVertices()

      select case(size(elementVertices))
        case(8)
          call self % generateTetrahedraFromHexahedron(element, vertices, edges, infoIdx, newElementIdx, tetrahedraPayloads)
          

      end select
      ! Deactivate element.
      call element % ptr % deactivate()

    end do
    
    ! Call superclass to build tetrahedra.
    call self % buildTetrahedraFromVertices(tetrahedraPayloads, edges, elements, faces)

  end subroutine triangulate

end module DompierreTriangulationMethod_class