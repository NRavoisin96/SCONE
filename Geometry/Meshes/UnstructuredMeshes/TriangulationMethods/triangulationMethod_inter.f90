module triangulationMethod_inter

  use edge_class,         only : buildEdgeInfo
  use edgeShelf_class,    only : edgeShelf
  use element_class,      only : buildElementInfo, elementBox
  use elementShelf_class, only : elementShelf
  use face_class,         only : buildFaceInfo, faceBox, orientatedFaceBox
  use faceShelf_class,    only : faceShelf
  use genericProcedures,  only : fatalError
  use numPrecision
  use universalVariables, only : FOURTH, NOT_PRESENT
  use vertexShelf_class,  only : vertexShelf

  implicit none
  private

  !!
  !!
  !!
  type, public, abstract :: triangulationMethod
    private
  contains
    procedure                        :: buildTetrahedraFromVertices
    procedure                        :: buildTrianglesFromVertices
    procedure(triangulate), deferred :: triangulate
  end type triangulationMethod

  !!
  !!
  !!
  abstract interface
    !!
    !!
    !!
    subroutine triangulate(self, edges, elements, faces, vertices)
      import                                 :: edgeShelf, elementShelf, faceShelf, triangulationMethod, vertexShelf
      class(triangulationMethod), intent(in) :: self
      type(edgeShelf), intent(inout)         :: edges
      type(elementShelf), intent(inout)      :: elements
      type(faceShelf), intent(inout)         :: faces
      type(vertexShelf), intent(inout)       :: vertices
    end subroutine triangulate

  end interface

contains
  !!
  !!
  !!
  subroutine buildTetrahedraFromVertices(self, tetrahedraInfos, edges, elements, faces)
    class(triangulationMethod), intent(in)              :: self
    type(buildElementInfo), dimension(:), intent(inout) :: tetrahedraInfos
    type(edgeShelf), intent(inout)                      :: edges
    type(elementShelf), intent(inout)                   :: elements
    type(faceShelf), intent(inout)                      :: faces
    type(buildFaceInfo)                                 :: faceInfo
    integer(shortInt)                                   :: edgeIdx, i, j, k, l, newEdgeIdx, newFaceIdx, triangleIdx
    type(buildEdgeInfo)                                 :: edgeInfo
    type(elementBox)                                    :: element
    type(orientatedFaceBox), dimension(:), allocatable  :: elementOrientatedFaces
    integer(shortInt), dimension(:), allocatable        :: childrenIdxs
    type(faceBox)                                       :: triangle
    real(defReal), dimension(3)                         :: outwardNormal
    character(*), parameter :: here = 'buildTetrahedraFromVertices (triangulationMethod_inter.f90)'

    ! Initialise variables.
    faceInfo % testNormal = .true.
    allocate(faceInfo % edges(3))
    allocate(faceInfo % vertices(3))
    newEdgeIdx = edges % getObjectsNumber()
    newFaceIdx = faces % getObjectsNumber()
    
    ! Loop through all tetrahedra to be generated.
    do i = 1, size(tetrahedraInfos)
      ! Compute centroid of the current tetrahedron.
      faceInfo % testCentroid = ZERO
      do j = 1, 4
        faceInfo % testCentroid = faceInfo % testCentroid + tetrahedraInfos(i) % vertices(j) % ptr % getCoordinates()

      end do
      faceInfo % testCentroid = FOURTH * faceInfo % testCentroid

      ! Build topology. Triangles first.
      allocate(tetrahedraInfos(i) % orientatedFaces(4))
      do j = 1, 4
        select case(j)
          case(1)
            faceInfo % vertices = tetrahedraInfos(i) % vertices([1, 2, 3])

          case(2)
            faceInfo % vertices = tetrahedraInfos(i) % vertices([1, 2, 4])

          case(3)
            faceInfo % vertices = tetrahedraInfos(i) % vertices([1, 3, 4])

          case(4)
            faceInfo % vertices = tetrahedraInfos(i) % vertices([2, 3, 4])

        end select
        triangleIdx = faces % getFaceIdxOrDefault(faceInfo % vertices, NOT_PRESENT)

        ! If triangle is not already present, we need to create it.
        if (triangleIdx == NOT_PRESENT) then
          newFaceIdx = newFaceIdx + 1
          faceInfo % idx = newFaceIdx

          do k = 1, 3
            select case(k)
              case(1)
                edgeInfo % vertices = faceInfo % vertices([1, 2])

              case(2)
                edgeInfo % vertices = faceInfo % vertices([1, 3])

              case(3)
                edgeInfo % vertices = faceInfo % vertices([2, 3])

            end select
            edgeIdx = edges % getEdgeIdxOrDefault(edgeInfo % vertices, NOT_PRESENT)
            if (edgeIdx == NOT_PRESENT) then
              newEdgeIdx = newEdgeIdx + 1
              edgeInfo % idx = newEdgeIdx
              call edges % initEdge(edgeInfo)

              ! Update connectivity.
              do l = 1, 2
                call edgeInfo % vertices(l) % ptr % addEdgeIdx(newEdgeIdx)

              end do
              edgeIdx = newEdgeIdx

            end if
            faceInfo % edges(k) = edges % getEdgeBox(edgeIdx)
            ! Update connectivity.
            call faceInfo % edges(k) % ptr % addFaceIdx(newFaceIdx)
            call faceInfo % vertices(k) % ptr % addFaceIdx(newFaceIdx)

          end do

          ! Since the triangle was created, the current tetrahedron is its owner.
          call faces % addFace(faceInfo)
          tetrahedraInfos(i) % orientatedFaces(j) % isOwner = .true.
          triangleIdx = newFaceIdx

        else
          ! If triangle already exists, check if the parent element owns it.
          element = elements % getElementBox(tetrahedraInfos(i) % parentIdx)
          elementOrientatedFaces = element % ptr % getOrientatedFaces()
          do k = 1, size(elementOrientatedFaces)
            childrenIdxs = elementOrientatedFaces(k) % face % ptr % getChildrenIdxs()
            if (any(childrenIdxs == triangleIdx)) then
              if (elementOrientatedFaces(k) % isOwner) tetrahedraInfos(i) % orientatedFaces(j) % isOwner = .true.

            end if

          end do

        end if

        ! Now retrieve the triangle and construct orientated face box for the tetrahedron.
        triangle = faces % getFaceBox(triangleIdx)
        tetrahedraInfos(i) % orientatedFaces(j) % face = triangle
        outwardNormal = triangle % ptr % getNormal()
        if (.not. tetrahedraInfos(i) % orientatedFaces(j) % isOwner) outwardNormal = -outwardNormal
        tetrahedraInfos(i) % orientatedFaces(j) % outwardNormal = outwardNormal

      end do

      ! Now build edges.
      allocate(tetrahedraInfos(i) % edges(6))
      do j = 1, 6
        select case(j)
          case(1)
            edgeIdx = edges % getEdgeIdxOrDefault([tetrahedraInfos(i) % vertices(1), &
                                                   tetrahedraInfos(i) % vertices(2)], NOT_PRESENT)

          case(2)
            edgeIdx = edges % getEdgeIdxOrDefault([tetrahedraInfos(i) % vertices(1), &
                                                   tetrahedraInfos(i) % vertices(3)], NOT_PRESENT)

          case(3)
            edgeIdx = edges % getEdgeIdxOrDefault([tetrahedraInfos(i) % vertices(1), &
                                                   tetrahedraInfos(i) % vertices(4)], NOT_PRESENT)
                                    
          case(4)
            edgeIdx = edges % getEdgeIdxOrDefault([tetrahedraInfos(i) % vertices(2), &
                                                   tetrahedraInfos(i) % vertices(3)], NOT_PRESENT)

          case(5)
            edgeIdx = edges % getEdgeIdxOrDefault([tetrahedraInfos(i) % vertices(2), &
                                                   tetrahedraInfos(i) % vertices(4)], NOT_PRESENT)

          case(6)
            edgeIdx = edges % getEdgeIdxOrDefault([tetrahedraInfos(i) % vertices(3), &
                                                   tetrahedraInfos(i) % vertices(4)], NOT_PRESENT)

        end select
        if (edgeIdx == NOT_PRESENT) call fatalError(here, 'Unable to find edge during tetrahedron initialisation.')
        tetrahedraInfos(i) % edges(j) = edges % getEdgeBox(edgeIdx)

      end do

    end do

    ! Ship payload to elementShelf.
    call elements % addElement(tetrahedraInfos)

  end subroutine buildTetrahedraFromVertices

  !!
  !!
  !!
  subroutine buildTrianglesFromVertices(self, triangleInfos, edges, faces)
    class(triangulationMethod), intent(in)           :: self
    type(buildFaceInfo), dimension(:), intent(inout) :: triangleInfos
    type(edgeShelf), intent(inout)                   :: edges
    type(faceShelf), intent(inout)                   :: faces
    integer(shortInt)                                :: edgeIdx, i, j, k, newEdgeIdx
    type(buildEdgeInfo)                              :: edgeInfo

    ! Initialise variables.
    newEdgeIdx = edges % getObjectsNumber()
    
    ! Generate all triangles from vertices.
    do i = 1, size(triangleInfos)
      allocate(triangleInfos(i) % edges(3))
      do j = 1, 3
        select case(j)
          case(1)
            edgeInfo % vertices = triangleInfos(i) % vertices([1, 2])

          case(2)
            edgeInfo % vertices = triangleInfos(i) % vertices([1, 3])

          case(3)
            edgeInfo % vertices = triangleInfos(i) % vertices([2, 3])

        end select
        edgeIdx = edges % getEdgeIdxOrDefault(edgeInfo % vertices, NOT_PRESENT)

        ! If edge is not already present, we need to create it.
        if (edgeIdx == NOT_PRESENT) then
          newEdgeIdx = newEdgeIdx + 1
          edgeInfo % idx = newEdgeIdx
          call edges % initEdge(edgeInfo)

          ! Update connectivity.
          do k = 1, 2
            call edgeInfo % vertices(k) % ptr % addEdgeIdx(newEdgeIdx)

          end do
          edgeIdx = newEdgeIdx

        end if
        triangleInfos(i) % edges(j) = edges % getEdgeBox(edgeIdx)
        ! Update connectivity.
        call triangleInfos(i) % vertices(j) % ptr % addFaceIdx(triangleInfos(i) % idx)
        call triangleInfos(i) % edges(j) % ptr % addFaceIdx(triangleInfos(i) % idx)

      end do

    end do

    ! Ship payload to faceShelf.
    call faces % addFace(triangleInfos)

  end subroutine buildTrianglesFromVertices

end module triangulationMethod_inter