module fanTriangulationMethod_inter

  use edge_class,                only : buildEdgeInfo, edgeBox
  use edgeShelf_class,           only : edgeShelf
  use face_class,                only : buildFaceInfo, faceBox
  use faceShelf_class,           only : faceShelf
  use numPrecision
  use triangulationMethod_inter, only : triangulationMethod
  use universalVariables,        only : NOT_PRESENT
  use vertex_class,              only : vertexBox

  implicit none
  private

  !!
  !!
  !!
  type, public, abstract, extends(triangulationMethod) :: fanTriangulationMethod
    private
  contains
    procedure :: decomposeFaces
  end type fanTriangulationMethod

contains
  !!
  !!
  !!
  subroutine decomposeFaces(self, edges, faces)
    class(fanTriangulationMethod), intent(in)      :: self
    type(edgeShelf), intent(inout)                 :: edges
    type(faceShelf), intent(inout)                 :: faces
    integer(shortInt)                              :: edgeIdx, faceIdx, i, j, k, l, newEdgeIdx, newFaceIdx, &
                                                      nFaces, nTriangles, nVertices
    integer(shortInt), dimension(:), allocatable   :: nTrianglesInFace
    type(faceBox)                                  :: face
    type(buildFaceInfo), dimension(:), allocatable :: faceInfos
    type(vertexBox), dimension(:), allocatable     :: faceVertices
    type(buildEdgeInfo)                            :: edgeInfo

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
    newFaceIdx = 0
    newEdgeIdx = edges % getObjectsNumber()
    do i = 1, nFaces
      face = faces % getFaceBox(i)
      faceIdx = face % ptr % getIdx()
      faceVertices = face % ptr % getVertices()
      do j = 1, nTrianglesInFace(i)
        newFaceIdx = newFaceIdx + 1
        faceInfos(newFaceIdx) % idx = nFaces + newFaceIdx
        faceInfos(newFaceIdx) % parentIdx = faceIdx
        faceInfos(newFaceIdx) % isBoundary = face % ptr % getIsBoundary()

        allocate(faceInfos(newFaceIdx) % vertices(3))
        faceInfos(newFaceIdx) % vertices(1) = faceVertices(1)
        faceInfos(newFaceIdx) % vertices(2) = faceVertices(j + 1)
        faceInfos(newFaceIdx) % vertices(3) = faceVertices(j + 2)

        allocate(faceInfos(newFaceIdx) % edges(3))
        do k = 1, 3
          edgeIdx = edges % getEdgeIdxOrDefault([faceInfos(newFaceIdx) % vertices(k), &
                                                 faceInfos(newFaceIdx) % vertices(merge(1, k + 1, k == 3))], NOT_PRESENT)

          ! If edge is not already present in shelf we need to create it.
          if (edgeIdx == NOT_PRESENT) then
            newEdgeIdx = newEdgeIdx + 1
            edgeInfo % idx = newEdgeIdx
            edgeInfo % vertices(1) = faceInfos(newFaceIdx) % vertices(k)
            edgeInfo % vertices(2) = faceInfos(newFaceIdx) % vertices(merge(1, k + 1, k == 3))
            call edges % initEdge(edgeInfo)

            ! Update connectivity.
            do l = 1, 2
              call edgeInfo % vertices(l) % ptr % addEdgeIdx(newEdgeIdx)

            end do

            edgeIdx = newEdgeIdx

          end if

          ! Now add the edge to the new triangle.
          faceInfos(newFaceIdx) % edges(k) = edges % getEdgeBox(edgeIdx)

          ! Update connectivity.
          call faceInfos(newFaceIdx) % vertices(k) % ptr % addFaceIdx(newFaceIdx)
          call faceInfos(newFaceIdx) % edges(k) % ptr % addFaceIdx(newFaceIdx)

        end do

        ! Add this child to the face.
        call face % ptr % addChildIdx(newFaceIdx)

      end do

      ! Deactivate face.
      call face % ptr % deactivate()

    end do

    ! Create new faces.
    call faces % addFace(faceInfos)

  end subroutine decomposeFaces

end module fanTriangulationMethod_inter