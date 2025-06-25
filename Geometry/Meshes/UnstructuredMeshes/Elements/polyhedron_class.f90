module polyhedron_class
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use edgeShelf_class,              only : edgeShelf
  use element_inter,                only : element, elementBox, kill_super => kill
  use face_inter,                   only : faceBox
  use faceShelf_class,              only : faceShelf
  use genericProcedures,            only : append, areEqual, computePyramidCentre, computePyramidVolume, &
                                           computeTetrahedronCentre, computeTetrahedronVolume, findCommon, &
                                           fatalError, numToChar
  use numPrecision
  use tetrahedron_class,            only : tetrahedron
  use triangle_class,               only : triangle
  use universalVariables,           only : SURF_TOL, INF, ZERO
  use vertexShelf_class,            only : vertexShelf
  
  implicit none
  private
  
  !!
  !! Element (cell) of an OpenFOAM mesh. Consists of a list of vertices and faces indices, as well
  !! as a list of tetrahedra indices into which the element is decomposed.
  !!
  !! Private members:
  !!   idx      -> Index of the element.
  !!   vertices -> Array of vertices indices making the element up.
  !!   faces    -> Array of faces indices making the element up.
  !!   Volume   -> Volume of the element.
  !!   Centroid -> Vector pointing to the centroid of the element.
  !!
  type, public, extends(element) :: polyhedron
    private
  contains
    ! Build procedures.
    procedure                    :: computeComponents
    procedure                    :: split
    ! Runtime procedures.
    procedure                    :: kill
  end type polyhedron

contains
  
  !! Subroutine 'computeVolumeAndCentroid'
  !!
  !! Basic description:
  !!   Computes the volume and volume-weighted centroid of the element.
  !!
  !! Detailed description:
  !!   First estimates the centroid of the element by taking an area-weighted average of the faces' 
  !!   centroids. Using this estimate, the element is decomposed into a number of pyramids whose 
  !!   apices are the estimated centroid of the element. Looping through all the faces, the volume 
  !!   and centroid of each pyramid is computed (see pyramid_class for more details); the centroid 
  !!   of the element is then obtained by taking a volume-weighted average of the pyramids' 
  !!   centroids.
  !!
  !! Arguments:
  !!   faces [in]    -> An array of face structures making the element up.
  !!   vertices [in] -> An array of vertex structures making the element up.
  !!
  !! Error:
  !!   fatalError if the volume is element is negative or infinite.
  !!
  pure subroutine computeComponents(self, faceIdxs, vertexIdxs, faces, vertices, centroid, volume)
    class(polyhedron), intent(inout)             :: self
    integer(shortInt), dimension(:), intent(in)  :: faceIdxs, vertexIdxs
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    real(defReal), dimension(3), intent(out)     :: centroid
    real(defReal), intent(out)                   :: volume
    integer(shortInt)                            :: i, nFaces, nVertices, absFaceIdx
    real(defReal)                                :: area, sumAreas, sumVolumes
    real(defReal), dimension(3)                  :: sumVolumesCentroid, C
    real(defReal), dimension(:, :), allocatable  :: array
    character(100), parameter                    :: Here = 'computeVolumeAndCentroid (element_class.f90)'

    ! Retrieve the number of faces and vertices in the element and initialise variables.
    nFaces = size(faceIdxs)
    nVertices = size(vertexIdxs)
    C = ZERO
    sumAreas = ZERO
    
    ! Loop through all faces and compute the element's approximate centroid
    ! by performing an area-weighted average of the different faces' centroids.
    do i = 1, nFaces
      ! Retrieve the area of the current face.
      absFaceIdx = abs(faceIdxs(i))
      area = faces % getFaceArea(absFaceIdx)
      ! Update the area-weighted centroid and the sum of faces' areas.
      C = C + faces % getFaceCentroid(absFaceIdx) * area
      sumAreas = sumAreas + area

    end do
    
    ! Using the approximate centroid, compute the actual centroid by performing a volume-weighted
    ! average of the different pyramids' centroids.
    sumVolumes = ZERO
    sumVolumesCentroid = ZERO
    allocate(array(2, 3))
    C = C / sumAreas
    
    ! Loop through all faces (pyramids).
    do i = 1, nFaces
      ! Retrieve the volume of the current pyramid and update the volume-weighted centroid and the sum of volumes.
      absFaceIdx = abs(faceIdxs(i))
      array(1, :) = faces % getFaceNormal(absFaceIdx) * faces % getFaceArea(absFaceIdx)
      array(2, :) = C - faces % getFaceCentroid(absFaceIdx)
      volume = computePyramidVolume(array)

      array(1, :) = 3.0_defReal * faces % getFaceCentroid(absFaceIdx)
      array(2, :) = C
      sumVolumesCentroid = sumVolumesCentroid + computePyramidCentre(array) * volume
      sumVolumes = sumVolumes + volume

    end do
    ! The volume of the element is simply the sum of volumes, while the centroid is the average of
    ! the volume-weighted sum.
    volume = sumVolumes
    centroid = sumVolumesCentroid / sumVolumes

  end subroutine computeComponents
  
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(polyhedron), intent(inout) :: self
    
    ! Element.
    call kill_super(self)

  end subroutine kill
  
  !! Subroutine 'split'
  !!
  !! Basic description:
  !!   Splits the element into a number of pyramids corresponding to the number of faces in the 
  !!   element.
  !!
  !! Detailed description:
  !!   The element is split into pyramids by subdividing it from its centroid: this then forms the
  !!   apex of each generated pyramid. For each face in the element, triangles are created by
  !!   joining each edge in the face to the common apex. During the creation of triangles a check is
  !!   made to ensure that their normal vectors point in the correct direction (from the pyramid of
  !!   lowest index to that of greatest index, just as OpenFOAM does for elements).
  !!
  !! Arguments:
  !!   faces [in]              -> A faceShelf.
  !!   edges [inout]           -> An edgeShelf.
  !!   vertices [inout]        -> A vertexShelf.
  !!   triangles [inout]       -> A triangleShelf.
  !!   pyramids [inout]        -> A pyramidShelf.
  !!   lastEdgeIdx [inout]     -> Index of the first free item in the edgeShelf.
  !!   lastTriangleIdx [inout] -> Index of the first free item in the triangleShelf.
  !!   lastPyramidIdx [inout]  -> Index of the first free item in the pyramidShelf.
  !!   lastVertexIdx [in]      -> Index of the first free item in the vertexShelf.
  !!
  subroutine split(self, faces, lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx, &
                        newEdges, newFaces, newVertices, tetrahedra, triangles)
    class(polyhedron), intent(inout)              :: self
    type(faceShelf), intent(inout)                :: faces, newFaces
    integer(shortInt), intent(inout)              :: lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx
    type(edgeShelf), intent(inout)                :: newEdges
    type(vertexShelf), intent(inout)              :: newVertices
    type(elementBox), dimension(:), intent(inout) :: tetrahedra
    type(faceBox), dimension(:), intent(inout)    :: triangles
    integer(shortInt)                             :: i, j, k, l, faceIdx, absFaceIdx, faceTriangleIdx, triangleVertexIdx, &
                                                     nFaces, nVertices, edgeIdx, commonTriangleIdx
    integer(shortInt), dimension(:), allocatable  :: faceIdxs, faceTriangleIdxs, vertexIdxs
    integer(shortInt), dimension(2)               :: edgeVertexIdxs
    integer(shortInt), dimension(3)               :: triangleVertexIdxs, testVertexIdxs
    integer(shortInt), dimension(4)               :: triangleIdxs
    integer(shortInt), dimension(6)               :: edgeIdxs
    real(defReal)                                 :: volume
    real(defReal), dimension(3)                   :: centroid
    real(defReal), dimension(4, 3)                :: array
    real(defReal), dimension(:, :), allocatable   :: vertexCoords
    type(axisAlignedBoundingBox)                  :: boundingBox
    
    ! Initialise a new vertex corresponding to the centroid of the polyhedron.
    lastNewVertexIdx = lastNewVertexIdx + 1
    call newVertices % initVertex(lastNewVertexIdx, self % getCentroid())
    array(4, :) = self % getCentroid()

    ! Retrieve the indices of the vertices in the polyhedron and compute the number of vertices.
    vertexIdxs = self % getVertexIdxs()
    nVertices = size(vertexIdxs)

    ! Loop through all the vertices in the polyhedron and create new edges joining each vertex to the centroid.
    do i = 1, nVertices
      lastNewEdgeIdx = lastNewEdgeIdx + 1
      edgeVertexIdxs = [vertexIdxs(i), lastNewVertexIdx]
      call newEdges % initEdge(lastNewEdgeIdx, [vertexIdxs(i), lastNewVertexIdx])

      do j = 1, 2
        call newVertices % addEdgeIdxToVertex(edgeVertexIdxs(j), lastNewEdgeIdx)

      end do

    end do

    ! Retrieve the indices of the faces in the polyhedron and compute the number of faces.
    faceIdxs = self % getFaceIdxs()
    nFaces = size(faceIdxs)

    ! Loop through all the faces in the polyhedron.
    do i = 1, nFaces
      ! Retrieve the current face index and create its absolute value.
      faceIdx = faceIdxs(i)
      absFaceIdx = abs(faceIdx)

      ! Retrieve the indices of the triangles created from the current face.
      faceTriangleIdxs = faces % getFaceTriangleIdxs(absFaceIdx)
      do j = 1, size(faceTriangleIdxs)
        ! Increment lastNewElementIdx and retrieve the indices of the vertices in the triangle.
        lastNewElementIdx = lastNewElementIdx + 1
        faceTriangleIdx = faceTriangleIdxs(j)
        call newFaces % addElementIdxToFace(faceTriangleIdx, lastNewElementIdx)
        
        triangleVertexIdxs = triangles(faceTriangleIdx) % item % getVertexIdxs()
        edgeIdxs(1:3) = triangles(faceTriangleIdx) % item % getEdgeIdxs()

        ! Create the indices of the vertices in the new tetrahedron.
        vertexIdxs = [triangleVertexIdxs, lastNewVertexIdx]
        call newVertices % addElementIdxToVertex(lastNewVertexIdx, lastNewElementIdx)
        triangleIdxs(1) = sign(faceTriangleIdx, faceIdx)

        ! Create remaining array entries and update mesh connectivity.
        do k = 1, 3
          triangleVertexIdx = triangleVertexIdxs(k)
          array(k, :) = newVertices % getVertexCoordinates(triangleVertexIdx)
          call newVertices % addElementIdxToVertex(triangleVertexIdx, lastNewElementIdx)

        end do

        ! Compute tetrahedron centroid and volume.
        centroid = computeTetrahedronCentre(array)
        volume = computeTetrahedronVolume(array)

        ! Loop through all the remaining faces in the new tetrahedron.
        if (allocated(vertexCoords)) deallocate(vertexCoords)
        allocate(vertexCoords(3, 3))
        do k = 1, 3
          edgeIdxs(k + 3) = newVertices % findCommonEdgeIdx(triangleVertexIdxs(k), lastNewVertexIdx)
          call newEdges % addElementIdxToEdge(edgeIdxs(k), lastNewElementIdx)
          call newEdges % addElementIdxToEdge(edgeIdxs(k + 3), lastNewElementIdx)
          
          ! Check if a triangle containing the three vertices already exists.
          testVertexIdxs = [triangleVertexIdxs(k), triangleVertexIdxs(mod(k, 3) + 1), lastNewVertexIdx]
          commonTriangleIdx = newVertices % findCommonFaceIdx(testVertexIdxs)
          if (commonTriangleIdx > 0) then
            call newFaces % addElementIdxToFace(commonTriangleIdx, lastNewElementIdx)
            triangleIdxs(k + 1) = -commonTriangleIdx
            cycle

          end if

          ! Create a new internal triangle.
          lastNewFaceIdx = lastNewFaceIdx + 1
          allocate(triangle :: triangles(lastNewFaceIdx) % item)
          triangleIdxs(k + 1) = lastNewFaceIdx

          ! Update mesh connectivity.
          do l = 1, 3
            edgeIdx = newVertices % findCommonEdgeIdx(testVertexIdxs(l), testVertexIdxs(mod(l, 3) + 1))
            call newEdges % addFaceIdxToEdge(edgeIdx, lastNewFaceIdx)
            call triangles(lastNewFaceIdx) % item % addEdgeIdx(edgeIdx)
            call newVertices % addFaceIdxToVertex(testVertexIdxs(l), lastNewFaceIdx)

            ! Retrieve the current vertex coordinates and update the triangle's bounding box.
            vertexCoords(:, l) = newVertices % getVertexCoordinates(testVertexIdxs(l))

          end do

          ! Build the new triangle and add it into the new faceShelf.
          call boundingBox % computeBounds(vertexCoords)
          call triangles(lastNewFaceIdx) % item % build(lastNewFaceIdx, 0, .false., testVertexIdxs, newVertices, 'Triangle', &
          boundingBox, centroid)
          call triangles(lastNewFaceIdx) % item % addElementIdx(lastNewElementIdx)
          call newFaces % addFace(lastNewFaceIdx, triangles(lastNewFaceIdx))

        end do

        ! Initialise new tetrahedron in the shelf.
        allocate(tetrahedron :: tetrahedra(lastNewElementIdx) % item)
        if (allocated(vertexCoords)) deallocate(vertexCoords)
        allocate(vertexCoords(3, 4))
        do k = 1, 4
          vertexCoords(:, k) = newVertices % getVertexCoordinates(vertexIdxs(k))

        end do
        call boundingBox % computeBounds(vertexCoords)
        call tetrahedra(lastNewElementIdx) % item % init(lastNewElementIdx, self % getIdx(), triangleIdxs, vertexIdxs, &
                                                         centroid, volume, .true., 'Tetrahedron', boundingBox, edgeIdxs)

      end do

    end do
  
  end subroutine split

end module polyhedron_class