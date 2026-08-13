module triangulationMethod_inter

  use element_class,                 only : buildElementPayload, elementBox
  use extentTopologicalObject_inter, only : buildExtentTopologicalObjectPayload
  use face_class,                    only : buildFacePayload, faceBox, orientatedFaceBox
  use genericProcedures,             only : fatalError
  use numPrecision
  use ratint
  use topologicalObjectShelf_class,  only : topologicalObjectShelf
  use universalVariables,            only : FOURTH, NOT_PRESENT

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
      import                                      :: topologicalObjectShelf, triangulationMethod
      class(triangulationMethod), intent(in)      :: self
      type(topologicalObjectShelf), intent(inout) :: edges, elements, faces, vertices
    end subroutine triangulate

  end interface

contains
  !!
  !!
  !!
  subroutine buildTetrahedraFromVertices(self, tetrahedraPayloads, edges, elements, faces)
    class(triangulationMethod), intent(in)                 :: self
    type(buildElementPayload), dimension(:), intent(inout) :: tetrahedraPayloads
    type(topologicalObjectShelf), intent(inout)            :: edges, elements, faces
    integer(shortInt)                                      :: edgeIdx, i, j, k, newEdgeIdx, newFaceIdx, triangleIdx
    integer(shortInt), dimension(:), allocatable           :: childrenIdxs
    real(defReal), dimension(3)                            :: outwardNormal
    type(buildExtentTopologicalObjectPayload)              :: edgePayload
    type(buildFacePayload)                                 :: facePayload
    type(elementBox)                                       :: element
    type(faceBox)                                          :: triangle
    type(orientatedFaceBox), dimension(:), allocatable     :: elementOrientatedFaces
    type(ratint_t), dimension(3)                           :: rationalOutwardNormal
    character(*), parameter :: here = 'buildTetrahedraFromVertices (triangulationMethod_inter.f90)'

    ! Initialise variables.
    facePayload % testNormal = .true.
    allocate(facePayload % edges(3))
    allocate(facePayload % vertices(3))
    newEdgeIdx = edges % getObjectsNumber()
    newFaceIdx = faces % getObjectsNumber()
    
    ! Loop through all tetrahedra to be generated.
    do i = 1, size(tetrahedraPayloads)
      ! Add this tetrahedron to its parent.
      element = elements % getElementBox(tetrahedraPayloads(i) % parentIdx)
      call element % ptr % addChildIdx(tetrahedraPayloads(i) % idx)

      ! Compute centroid of the current tetrahedron.
      facePayload % testCentroid = ZERO
      do j = 1, 4
        facePayload % testCentroid = facePayload % testCentroid + tetrahedraPayloads(i) % vertices(j) % ptr % getCoordinates()

      end do
      facePayload % testCentroid = FOURTH * facePayload % testCentroid

      ! Build topology. Triangles first.
      allocate(tetrahedraPayloads(i) % orientatedFaces(4))
      do j = 1, 4
        select case(j)
          case(1)
            facePayload % vertices = tetrahedraPayloads(i) % vertices([1, 2, 3])

          case(2)
            facePayload % vertices = tetrahedraPayloads(i) % vertices([1, 2, 4])

          case(3)
            facePayload % vertices = tetrahedraPayloads(i) % vertices([1, 3, 4])

          case(4)
            facePayload % vertices = tetrahedraPayloads(i) % vertices([2, 3, 4])

        end select
        triangleIdx = faces % getObjectIdxOrDefault(facePayload % vertices, NOT_PRESENT)

        ! If triangle is not already present, we need to create it.
        if (triangleIdx == NOT_PRESENT) then
          newFaceIdx = newFaceIdx + 1
          facePayload % idx = newFaceIdx

          do k = 1, 3
            select case(k)
              case(1)
                edgePayload % vertices = facePayload % vertices([1, 2])

              case(2)
                edgePayload % vertices = facePayload % vertices([1, 3])

              case(3)
                edgePayload % vertices = facePayload % vertices([2, 3])

            end select
            edgeIdx = edges % getObjectIdxOrDefault(edgePayload % vertices, NOT_PRESENT)
            if (edgeIdx == NOT_PRESENT) then
              newEdgeIdx = newEdgeIdx + 1
              edgePayload % idx = newEdgeIdx
              call edges % initObject(edgePayload)
              edgeIdx = newEdgeIdx

            end if
            facePayload % edges(k) = edges % getEdgeBox(edgeIdx)

          end do

          ! Since the triangle was created, the current tetrahedron is its owner.
          call faces % initObject(facePayload)
          tetrahedraPayloads(i) % orientatedFaces(j) % isOwner = .true.
          triangleIdx = newFaceIdx

        else
          ! If triangle already exists, check if the parent element owns it.
          elementOrientatedFaces = element % ptr % getOrientatedFaces()
          do k = 1, size(elementOrientatedFaces)
            childrenIdxs = elementOrientatedFaces(k) % face % ptr % getChildrenIdxs()
            if (any(childrenIdxs == triangleIdx)) then
              if (elementOrientatedFaces(k) % isOwner) tetrahedraPayloads(i) % orientatedFaces(j) % isOwner = .true.

            end if

          end do

        end if

        ! Now retrieve the triangle and construct orientated face box for the tetrahedron.
        triangle = faces % getFaceBox(triangleIdx)
        tetrahedraPayloads(i) % orientatedFaces(j) % face = triangle
        outwardNormal = triangle % ptr % getNormal()
        rationalOutwardNormal = triangle % ptr % getRatintNormal()
        if(.not. tetrahedraPayloads(i) % orientatedFaces(j) % isOwner) then
          outwardNormal = -outwardNormal
          call swapSign(rationalOutwardNormal)

        end if
        tetrahedraPayloads(i) % orientatedFaces(j) % outwardNormal = outwardNormal
        tetrahedraPayloads(i) % orientatedFaces(j) % ratintOutwardNormal = rationalOutwardNormal

      end do

      ! Now build edges.
      allocate(tetrahedraPayloads(i) % edges(6))
      do j = 1, 6
        select case(j)
          case(1)
            edgeIdx = edges % getObjectIdxOrDefault([tetrahedraPayloads(i) % vertices(1), &
                                                     tetrahedraPayloads(i) % vertices(2)], NOT_PRESENT)

          case(2)
            edgeIdx = edges % getObjectIdxOrDefault([tetrahedraPayloads(i) % vertices(1), &
                                                     tetrahedraPayloads(i) % vertices(3)], NOT_PRESENT)

          case(3)
            edgeIdx = edges % getObjectIdxOrDefault([tetrahedraPayloads(i) % vertices(1), &
                                                     tetrahedraPayloads(i) % vertices(4)], NOT_PRESENT)
                                    
          case(4)
            edgeIdx = edges % getObjectIdxOrDefault([tetrahedraPayloads(i) % vertices(2), &
                                                     tetrahedraPayloads(i) % vertices(3)], NOT_PRESENT)

          case(5)
            edgeIdx = edges % getObjectIdxOrDefault([tetrahedraPayloads(i) % vertices(2), &
                                                     tetrahedraPayloads(i) % vertices(4)], NOT_PRESENT)

          case(6)
            edgeIdx = edges % getObjectIdxOrDefault([tetrahedraPayloads(i) % vertices(3), &
                                                     tetrahedraPayloads(i) % vertices(4)], NOT_PRESENT)

        end select
        if (edgeIdx == NOT_PRESENT) call fatalError(here, 'Unable to find edge during tetrahedron initialisation.')
        tetrahedraPayloads(i) % edges(j) = edges % getEdgeBox(edgeIdx)

      end do

    end do

    ! Ship payload to elementShelf.
    call elements % initObject(tetrahedraPayloads)

  end subroutine buildTetrahedraFromVertices

  !!
  !!
  !!
  subroutine buildTrianglesFromVertices(self, trianglePayloads, edges, faces)
    class(triangulationMethod), intent(in)              :: self
    type(buildFacePayload), dimension(:), intent(inout) :: trianglePayloads
    type(topologicalObjectShelf), intent(inout)         :: edges, faces
    integer(shortInt)                                   :: edgeIdx, i, j, newEdgeIdx
    type(buildExtentTopologicalObjectPayload)           :: edgePayload

    ! Initialise variables.
    newEdgeIdx = edges % getObjectsNumber()
    
    ! Generate all triangles from vertices.
    do i = 1, size(trianglePayloads)
      allocate(trianglePayloads(i) % edges(3))
      do j = 1, 3
        select case(j)
          case(1)
            edgePayload % vertices = trianglePayloads(i) % vertices([1, 2])

          case(2)
            edgePayload % vertices = trianglePayloads(i) % vertices([1, 3])

          case(3)
            edgePayload % vertices = trianglePayloads(i) % vertices([2, 3])

        end select
        edgeIdx = edges % getObjectIdxOrDefault(edgePayload % vertices, NOT_PRESENT)

        ! If edge is not already present, we need to create it.
        if (edgeIdx == NOT_PRESENT) then
          newEdgeIdx = newEdgeIdx + 1
          edgePayload % idx = newEdgeIdx
          call edges % initObject(edgePayload)
          edgeIdx = newEdgeIdx

        end if
        trianglePayloads(i) % edges(j) = edges % getEdgeBox(edgeIdx)

      end do

    end do

    ! Ship payload to faceShelf.
    call faces % initObject(trianglePayloads)

  end subroutine buildTrianglesFromVertices

end module triangulationMethod_inter