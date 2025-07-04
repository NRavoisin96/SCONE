module centroidTriangulationMethod_class

  use edge_class,                   only : buildEdgeInfo, edgeBox
  use edgeShelf_class,              only : edgeShelf
  use element_class,                only : buildElementInfo, elementBox
  use elementShelf_class,           only : elementShelf
  use face_class,                   only : buildFaceInfo, faceBox, orientatedFaceBox
  use faceShelf_class,              only : faceShelf
  use fanTriangulationMethod_inter, only : fanTriangulationMethod
  use genericProcedures,            only : fatalError, numToChar
  use numPrecision
  use universalVariables,           only : NOT_PRESENT
  use vertex_class,                 only : vertexBox
  use vertexShelf_class,            only : vertexShelf

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(fanTriangulationMethod) :: centroidTriangulationMethod
    private
  contains
    procedure :: buildTetrahedraVertices
    procedure :: triangulate
  end type centroidTriangulationMethod

contains
  !!
  !!
  !!
  subroutine buildTetrahedraVertices(self, elements, faces, vertices, infos)
    class(centroidTriangulationMethod), intent(in)      :: self
    type(elementShelf), intent(in)                      :: elements
    type(faceShelf), intent(in)                         :: faces
    type(vertexShelf), intent(inout)                    :: vertices
    type(buildElementInfo), dimension(:), intent(inout) :: infos
    integer(shortInt)                                   :: i, j, k, l, nElements, newElementIdx, newVertexIdx
    type(elementBox)                                    :: element
    type(orientatedFaceBox), dimension(:), allocatable  :: elementOrientatedFaces
    integer(shortInt), dimension(:), allocatable        :: childrenIdxs
    type(faceBox)                                       :: triangle
    real(defReal), dimension(3)                         :: outwardNormal

    ! Initialise nElements and newVertexIdx.
    nElements = elements % getObjectsNumber()
    newElementIdx = nElements
    newVertexIdx = vertices % getObjectsNumber()
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
          newElementIdx = newElementIdx + 1
          infos(newElementIdx) % idx = newElementIdx
          infos(newElementIdx) % parentIdx = element % ptr % getIdx()

          triangle = faces % getFaceBox(childrenIdxs(k))
          
          ! Allocate number of vertices for the current tetrahedron and populate them.
          allocate(infos(newElementIdx) % vertices(4))
          infos(newElementIdx) % vertices(1:3) = triangle % ptr % getVertices()
          infos(newElementIdx) % vertices(4) = vertices % getVertexBox(newVertexIdx)

          ! Update connectivity for vertices.
          do l = 1, 4
            call infos(newElementIdx) % vertices(l) % ptr % addElementIdx(newElementIdx)

          end do

          ! Allocate number of faces for the current tetrahedron and create the first face.
          allocate(infos(newElementIdx) % orientatedFaces(4))
          infos(newElementIdx) % orientatedFaces(1) % face = triangle
          if (elementOrientatedFaces(j) % isOwner) infos(newElementIdx) % orientatedFaces(1) % isOwner = .true.
          outwardNormal = triangle % ptr % getNormal()
          if (.not. infos(newElementIdx) % orientatedFaces(1) % isOwner) outwardNormal = -outwardNormal
          infos(newElementIdx) % orientatedFaces(1) % outwardNormal = outwardNormal

          ! Update connectivity for this face.
          call infos(newElementIdx) % orientatedFaces(1) % face % ptr % addElementIdx(newElementIdx)

          ! Allocate number of edges for the current tetrahedron and populate the first three.
          allocate(infos(newElementIdx) % edges(6))
          infos(newElementIdx) % edges(1:3) = triangle % ptr % getEdges()

          ! Update connectivity for the three edges.
          do l = 1, 3
            call infos(newElementIdx) % edges(l) % ptr % addElementIdx(newElementIdx)

          end do

        end do

      end do
      ! Deactivate current element.
      call element % ptr % deactivate()

    end do

  end subroutine buildTetrahedraVertices

  !!
  !!
  !!
  subroutine triangulate(self, edges, elements, faces, vertices)
    class(centroidTriangulationMethod), intent(in)     :: self
    type(edgeShelf), intent(inout)                     :: edges
    type(elementShelf), intent(inout)                  :: elements
    type(faceShelf), intent(inout)                     :: faces
    type(vertexShelf), intent(inout)                   :: vertices
    integer(shortInt)                                  :: edgeIdx, i, j, k, l, nElements, newEdgeIdx, newFaceIdx, &
                                                          nTetrahedra, triangleIdx
    type(elementBox)                                   :: element
    type(orientatedFaceBox), dimension(:), allocatable :: elementOrientatedFaces
    type(faceBox)                                      :: triangle
    type(buildElementInfo), dimension(:), allocatable  :: elementInfos
    type(buildFaceInfo)                                :: faceInfo
    type(buildEdgeInfo)                                :: edgeInfo
    real(defReal), dimension(3)                        :: outwardNormal
    character(*), parameter                            :: here = 'triangulate (centroidTriangulationMethod_class.f90)'

    ! Triangulate faces first.
    call self % decomposeFaces(edges, faces)

    ! Count the number of tetrahedra to be genetated.
    nTetrahedra = 0
    nElements = elements % getObjectsNumber()
    do i = 1, nElements
      element = elements % getElementBox(i)
      if (size(element % ptr % getVertices()) == 4) cycle
      
      elementOrientatedFaces = element % ptr % getOrientatedFaces()
      do j = 1, size(elementOrientatedFaces)
        nTetrahedra = nTetrahedra + size(elementOrientatedFaces(j) % face % ptr % getChildrenIdxs())

      end do

    end do

    ! Return early if no tetrahedra need to be generated. Allocate elementInfos otherwise.
    if (nTetrahedra == 0) return
    allocate(elementInfos(nTetrahedra))

    ! Populate vertices within each elementInfo.
    call self % buildTetrahedraVertices(elements, faces, vertices, elementInfos)

    ! Initialise faceInfo.
    allocate(faceInfo % vertices(3))
    allocate(faceInfo % edges(3))
    faceInfo % testNormal = .true.

    ! Generate new tetrahedra.
    newEdgeIdx = edges % getObjectsNumber()
    newFaceIdx = faces % getObjectsNumber()
    do i = 1, nTetrahedra
      ! Compute faceInfo % testCentroid.
      do j = 1, 4
        faceInfo % testCentroid = faceInfo % testCentroid + elementInfos(i) % vertices(j) % ptr % getCoordinates()

      end do
      faceInfo % testCentroid = FOURTH * faceInfo % testCentroid
      
      ! Create three internal triangles and their edges.
      faceInfo % vertices(3) = elementInfos(i) % vertices(4)
      do j = 1, 3
        ! Assign the three vertices for the current triangle.
        faceInfo % vertices(1) = elementInfos(i) % vertices(j)
        faceInfo % vertices(2) = elementInfos(i) % vertices(merge(1, j + 1, j == 3))

        ! Check if the current triangle already exists.
        triangleIdx = faces % getFaceIdxOrDefault(faceInfo % vertices, NOT_PRESENT)

        ! If triangle is not already present, we need to create it.
        if (triangleIdx == NOT_PRESENT) then
          newFaceIdx = newFaceIdx + 1
          faceInfo % idx = newFaceIdx

          ! Create edges.
          do k = 1, 3
            edgeIdx = edges % getEdgeIdxOrDefault([faceInfo % vertices(merge(2, 1, 2 < k)), &
                                                   faceInfo % vertices(merge(3, 2, 1 < k))], NOT_PRESENT)

            ! If edge is not already present we need to create it.
            if (edgeIdx == NOT_PRESENT) then
              newEdgeIdx = newEdgeIdx + 1
              edgeInfo % idx = newEdgeIdx
              edgeInfo % vertices(1) = faceInfo % vertices(merge(2, 1, 2 < k))
              edgeInfo % vertices(2) = faceInfo % vertices(merge(3, 2, 1 < k))
              call edges % initEdge(edgeInfo)

              ! Update connectivity.
              do l = 1, 2
                call edgeInfo % vertices(l) % ptr % addEdgeIdx(newEdgeIdx)

              end do

              edgeIdx = newEdgeIdx

            end if

            ! Now retrieve the edge from the shelf.
            faceInfo % edges(k) = edges % getEdgeBox(edgeIdx)

            ! Update connectivity.
            call faceInfo % vertices(k) % ptr % addFaceIdx(newFaceIdx)
            call faceInfo % edges(k) % ptr % addFaceIdx(newFaceIdx)
  
          end do

          ! Create new triangle.
          call faces % addFace(faceInfo)

          ! Since the triangle was created, the current tetrahedron is its owner.
          elementInfos(i) % orientatedFaces(j + 1) % isOwner = .true.
          triangleIdx = newFaceIdx

        end if

        ! Now retrieve the triangle and construct orientated face box for the tetrahedron.
        triangle = faces % getFaceBox(triangleIdx)
        elementInfos(i) % orientatedFaces(j + 1) % face = triangle
        outwardNormal = triangle % ptr % getNormal()
        if (.not. elementInfos(i) % orientatedFaces(j + 1) % isOwner) outwardNormal = -outwardNormal
        elementInfos(i) % orientatedFaces(j + 1) % outwardNormal = outwardNormal

        ! Update connectivity for this face.
        call elementInfos(i) % orientatedFaces(j + 1) % face % ptr % addElementIdx(nElements + i)

        ! Now assign the remaining edges to the tetrahedron.
        do k = 1, 3
          edgeIdx = edges % getEdgeIdxOrDefault([elementInfos(i) % vertices(k), elementInfos(i) % vertices(4)], NOT_PRESENT)

          ! If edge was not found call fatalError.
          if (edgeIdx == NOT_PRESENT) call fatalError(here, 'Unable to retrieve tetrahedron edge after initialisation.')
          elementInfos(i) % edges(k + 3) = edges % getEdgeBox(edgeIdx)

          ! Update connectivity.
          call elementInfos(i) % edges(k + 3) % ptr % addElementIdx(nElements + i)
        
        end do

      end do

    end do
    ! Add new tetrahedra to the shelf.
    call elements % addElement(elementInfos)

  end subroutine triangulate

end module centroidTriangulationMethod_class