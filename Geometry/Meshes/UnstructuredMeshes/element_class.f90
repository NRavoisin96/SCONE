module element_class
  
  use edgeShelf_class,     only : edgeShelf
  use face_class,          only : face
  use faceShelf_class,     only : faceShelf
  use genericProcedures,   only : append, areEqual, computePyramidCentre, computePyramidVolume, &
                                  computeTetrahedronCentre, computeTetrahedronVolume, findCommon, &
                                  fatalError, numToChar
  use numPrecision
  use pyramid_class,       only : pyramid
  use pyramidShelf_class,  only : pyramidShelf
  use tetrahedron_class,   only : tetrahedron
  use triangle_class,      only : triangle
  use triangleShelf_class, only : triangleShelf
  use universalVariables,  only : SURF_TOL, INF, ZERO
  use vertex_class,        only : vertex
  use vertexShelf_class,   only : vertexShelf
  
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
  type, public :: element
    private
    integer(shortInt)                            :: idx = 0
    integer(shortInt), dimension(:), allocatable :: edgeIdxs, faceIdxs, vertexIdxs, tetrahedronIdxs
    real(defReal)                                :: volume = ZERO
    real(defReal), dimension(3)                  :: centroid = ZERO
  contains
    ! Build procedures.
    procedure                                    :: addEdgeIdx
    procedure                                    :: addFaceToElement
    procedure                                    :: addVertexToElement
    procedure                                    :: computeVolumeAndCentroid
    procedure                                    :: isConvex
    procedure                                    :: setIdx
    procedure                                    :: split
    ! Runtime procedures.
    procedure                                    :: computeIntersectedFace
    procedure                                    :: computePotentialFaces
    procedure                                    :: getCentroid
    procedure                                    :: getEdgeIdxs
    procedure                                    :: getFaces
    procedure                                    :: getIdx
    procedure                                    :: getVertices
    procedure                                    :: getVolume
    procedure                                    :: kill
    procedure                                    :: testForInclusion
  end type element

contains

  !! Subroutine 'addEdgeIdx'
  !!
  !! Basic description:
  !!   Adds the index of an edge sharing the element.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the edge.
  !!
  elemental subroutine addEdgeIdx(self, idx)
    class(element), intent(inout)  :: self
    integer(shortInt), intent(in)  :: idx

    call append(self % edgeIdxs, idx, .true.)

  end subroutine addEdgeIdx
  
  !! Subroutine 'addFaceToElement'
  !!
  !! Basic description:
  !!   Adds the index of a face belonging to the element.
  !!
  !! Arguments:
  !!   faceIdx [in] -> Index of the face.
  !!
  elemental subroutine addFaceToElement(self, faceIdx)
    class(element), intent(inout) :: self
    integer(shortInt), intent(in) :: faceIdx
    
    call append(self % faceIdxs, faceIdx)

  end subroutine addFaceToElement
  
  !! Subroutine 'addVertexToElement'
  !!
  !! Basic description:
  !!   Adds the index of a vertex belonging to the element. Only adds it if the index is not already
  !!   present.
  !!
  !! Arguments:
  !!   vertexIdx [in] -> Index of the vertex.
  !!
  elemental subroutine addVertexToElement(self, vertexIdx)
    class(element), intent(inout) :: self
    integer(shortInt), intent(in) :: vertexIdx

    call append(self % vertexIdxs, vertexIdx, .true.)

  end subroutine addVertexToElement

  !! Subroutine 'computeIntersectedFace'
  !!
  !! Basic description:
  !!   Computes the element's face which is intersected by a line segment.
  !!
  !! Detailed description:
  !!   See Macpherson, et al. (2009). DOI: 10.1002/cnm.1128.
  !!
  !! Arguments:
  !!   startPos [in]            -> Beginning of the line segment.
  !!   endPos [in]              -> End of the line segment.
  !!   potentialFaceIdxs [in]   -> Array of potential faces intersected by the line segment.
  !!   intersectedFaceIdx [out] -> Index of the face intersected by the line segment.
  !!   lambda [out]             -> Fraction of the line segment to be traversed before reaching the
  !!                               intersection point.
  !!   faces [in]               -> A faceShelf.
  !!
  !! TODO: Review this. Especially for cases where a particle may end up on an edge or at a vertex.
  pure subroutine computeIntersectedFace(self, startPos, endPos, potentialFaceIdxs, &
                                         intersectedFaceIdx, lambda, faces)
    class(element), intent(in)                  :: self
    real(defReal), dimension(3), intent(in)     :: startPos, endPos
    integer(shortInt), dimension(:), intent(in) :: potentialFaceIdxs
    integer(shortInt), intent(out)              :: intersectedFaceIdx
    real(defReal), intent(out)                  :: lambda
    type(faceShelf), intent(in)                 :: faces
    integer(shortInt)                           :: i, faceIdx, absFaceIdx
    real(defReal), dimension(3)                 :: centroid, normal
    real(defReal)                               :: faceLambda
    
    ! Initialise lambda = INF.
    lambda = INF
    
    ! Loop over all potentially intersected triangles.
    do i = 1, size(potentialFaceIdxs)
      ! Retrieve the current face's signed normal vector and flip it if necessary.
      faceIdx = potentialFaceIdxs(i)
      absFaceIdx = abs(faceIdx)
      normal = faces % getFaceNormal(faceIdx)
      
      ! Retrieve the current triangle's centre and compute triangleLambda.
      centroid = faces % getFaceCentroid(absFaceIdx)
      faceLambda = dot_product(centroid - startPos, normal) / dot_product(endPos - startPos, normal)
      
      ! If triangleLambda < lambda, update lambda and intersectedTriangleIdx.
      if (faceLambda < lambda) then
        lambda = faceLambda
        intersectedFaceIdx = absFaceIdx

      end if

    end do

  end subroutine computeIntersectedFace

  !! Function 'computePotentialFaces'
  !!
  !! Basic description:
  !!   Computes a list of indices of the potentially intersected faces using the element's 
  !!   centroid and the end of a line segment.
  !!
  !! Detailed description:
  !!   See Macpherson, et al. (2009). DOI: 10.1002/cnm.1128.
  !!
  !! Arguments:
  !!   endPos [in]       -> End of the line segment.
  !!   faces [in]        -> A faceShelf.
  !!
  !! Result:
  !!   potentialFaceIdxs -> Array of indices of the potential faces intersected by the line segment.
  !!
  pure function computePotentialFaces(self, endPos, faces) result(potentialFaceIdxs)
    class(element), intent(in)                   :: self
    real(defReal), dimension(3), intent(in)      :: endPos
    type(faceShelf), intent(in)                  :: faces
    integer(shortInt), dimension(:), allocatable :: potentialFaceIdxs
    integer(shortInt)                            :: i, faceIdx, absFaceIdx
    real(defReal)                                :: lambda
    real(defReal), dimension(3)                  :: centroid, faceCentroid, normal
    
    ! Retrieve element's centroid then loop over all faces in the tetrahedron.
    allocate(potentialFaceIdxs(0))
    centroid = self % centroid
    do i = 1, size(self % faceIdxs)
      ! Retrieve current face index and create its absolute index.
      faceIdx = self % faceIdxs(i)
      absFaceIdx = abs(faceIdx)
      
      ! Retrieve the signed normal vector of the current face.
      normal = faces % getFaceNormal(faceIdx)
      
      ! Retrieve the centre of the current face and compute lambda.
      faceCentroid = faces % getFaceCentroid(absFaceIdx)
      lambda = dot_product(faceCentroid - centroid, normal) / dot_product(endPos - centroid, normal)
      
      ! If ZERO <= lambda <= ONE, append the current face to the list of potentially intersected faces.
      if (ZERO <= lambda .and. lambda <= ONE) call append(potentialFaceIdxs, faceIdx)

    end do

  end function computePotentialFaces
  
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
  subroutine computeVolumeAndCentroid(self, vertices, faces)
    class(element), intent(inout)               :: self
    type(vertexShelf), intent(in)               :: vertices
    type(faceShelf), intent(in)                 :: faces
    integer(shortInt)                           :: i, nFaces, nVertices, absFaceIdx
    real(defReal)                               :: area, sumAreas, volume, sumVolumes
    real(defReal), dimension(3)                 :: sumVolumesCentroid, C
    real(defReal), dimension(:, :), allocatable :: array
    character(100), parameter                   :: Here = 'computeVolumeAndCentroid (element_class.f90)'
    
    ! Retrieve the number of vertices in the element.
    nVertices = size(self % vertexIdxs)
    
    ! If the element is a tetrahedron perform a direct computation to minimise round-off errors.
    if (nVertices == 4) then
      allocate(array(4, 3))
      do i = 1, 4
        array(i, :) = vertices % getVertexCoordinates(self % vertexIdxs(i))

      end do
      self % volume = computeTetrahedronVolume(array)
      self % centroid = computeTetrahedronCentre(array)

    else
      ! Retrieve the number of faces in the element and initialise variables.
      nFaces = size(self % faceIdxs)
      C = ZERO
      sumAreas = ZERO
      
      ! Loop through all faces and compute the element's approximate centroid
      ! by performing an area-weighted average of the different faces' centroids.
      do i = 1, nFaces
        ! Retrieve the area of the current face.
        absFaceIdx = abs(self % faceIdxs(i))
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
        absFaceIdx = abs(self % faceIdxs(i))
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
      self % volume = sumVolumes
      self % centroid = sumVolumesCentroid / sumVolumes

    end if
    
    ! Check that the volume is valid.
    if (self % volume <= ZERO) call fatalError(Here, 'Negative volume for element with index: '//numToChar(self % idx)//'.')
    if (self % volume >= INF) call fatalError(Here, 'Infinite volume for element with index: '//numToChar(self % idx)//'.')

  end subroutine computeVolumeAndCentroid
  
  !! Function 'isConvex'
  !!
  !! Basic description:
  !!   Checks whether the element is convex.
  !!
  !! Detailed description:
  !!    Convexity is checked by taking each vertex in the a given face and creating a vector 
  !!    connecting said vertex to each vertex in the element not in the current face. If the element 
  !!    is convex then all the vertices not in the current face must lie on the same side of the 
  !!    face, hence the dot product between the current face's normal vector and the test vector 
  !!    must be negative. If at any point the dot product is found to be positive the check is 
  !!    aborted.
  !!
  !! Arguments:
  !!   vertices [in] -> A vertexShelf.
  !!   faces [in]    -> A faceShelf.
  !!
  !! Result:
  !!   isIt          -> .true. if the element is convex.
  !!
  elemental function isConvex(self, vertices, faces) result(isIt)
    class(element), intent(in)                   :: self
    type(vertexShelf), intent(in)                :: vertices
    type(faceShelf), intent(in)                  :: faces
    integer(shortInt)                            :: i, j, k, faceIdx, absFaceIdx, vertexIdx
    integer(shortInt), dimension(:), allocatable :: faceVertexIdxs
    logical(defBool)                             :: isIt
    real(defReal), dimension(3)                  :: normal, faceVertexCoords
    
    ! Initialise isIt = .false.
    isIt = .false.
    
    ! Now loop through all the faces in the element.
    do i = 1, size(self % faceIdxs)
      ! Retrieve the current face's vertices and signed normal vector.
      faceIdx = self % faceIdxs(i)
      absFaceIdx = abs(faceIdx)
      faceVertexIdxs = faces % getFaceVertexIdxs(absFaceIdx)
      normal = faces % getFaceNormal(faceIdx)
      
      ! Loop through all the vertices in the current face.
      do j = 1, size(faceVertexIdxs)
        ! Retrieve the coordinates of the current face vertex and loop through all the vertices in the element.
        faceVertexCoords = vertices % getVertexCoordinates(faceVertexIdxs(j))
        do k = 1, size(self % vertexIdxs)
          ! Cycle to the next vertex if the current vertex index corresponds to the index of a vertex in the current face.
          vertexIdx = self % vertexIdxs(k)
          if (any(faceVertexIdxs == vertexIdx)) cycle
          
          ! Assemble the test vector and check if normal .dot. testVector > ZERO. If yes, the element
          ! is concave and we can return early.
          if (dot_product(normal, vertices % getVertexCoordinates(vertexIdx) - faceVertexCoords) > ZERO) return

        end do

      end do

    end do
    
    ! If reached this point the element is convex. Update isIt = .true.
    isIt = .true.

  end function isConvex
  
  !! Function 'getCentroid'
  !!
  !! Basic description:
  !!   Returns the centroid of the element.
  !!
  !! Result:
  !!   centroid -> A vector pointing to the centroid of the element.
  !!
  pure function getCentroid(self) result(centroid)
    class(element), intent(in)  :: self
    real(defReal), dimension(3) :: centroid
    
    centroid = self % centroid

  end function getCentroid

  !! Function 'getEdgeIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the edges in the element.
  !!
  !! Result:
  !!   edgeIdxs -> Indices of the edges in the element.
  !!
  pure function getEdgeIdxs(self) result(edgeIdxs)
    class(element), intent(in)                          :: self
    integer(shortInt), dimension(size(self % edgeIdxs)) :: edgeIdxs

    edgeIdxs = self % edgeIdxs

  end function getEdgeIdxs
  
  !! Function 'getFaces'
  !!
  !! Basic description:
  !!   Returns the indices of the faces in the element.
  !!
  !! Result:
  !!   faceIdxs -> An array listing the indices of the faces in the element.
  !!
  pure function getFaces(self) result(faceIdxs)
    class(element), intent(in)                          :: self
    integer(shortInt), dimension(size(self % faceIdxs)) :: faceIdxs
    
    faceIdxs = self % faceIdxs

  end function getFaces
  
  !! Function 'getIdx'
  !!
  !! Basic description:
  !!   Returns the index of the element.
  !!
  !! Result:
  !!   idx -> Index of the element.
  !!
  elemental function getIdx(self) result(idx)
    class(element), intent(in) :: self
    integer(shortInt)          :: idx
    
    idx = self % idx

  end function getIdx
  
  !! Function 'getVertices'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in the element.
  !!
  !! Result:
  !!   vertexIdxs -> An array listing indices of the vertices in the element.
  !!
  pure function getVertices(self) result(vertexIdxs)
    class(element), intent(in)                          :: self
    integer(shortInt), dimension(size(self % vertexIdxs)) :: vertexIdxs
    
    vertexIdxs = self % vertexIdxs

  end function getVertices
  
  !! Function 'getVolume'
  !!
  !! Basic description:
  !!   Returns the volume of the element.
  !!
  !! Result:
  !!   volume -> Volume of the element.
  !!
  elemental function getVolume(self) result(volume)
    class(element), intent(in) :: self
    real(defReal)              :: volume
    
    volume = self % volume

  end function getVolume
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(element), intent(inout) :: self
    
    self % idx = 0
    self % volume = ZERO
    self % centroid = ZERO
    if (allocated(self % edgeIdxs)) deallocate(self % edgeIdxs)
    if (allocated(self % vertexIdxs)) deallocate(self % vertexIdxs)
    if (allocated(self % faceIdxs)) deallocate(self % faceIdxs)
    if (allocated(self % tetrahedronIdxs)) deallocate(self % tetrahedronIdxs)

  end subroutine kill
  
  !! Subroutine 'setIdx'
  !!
  !! Basic description:
  !!   Sets the index of the element.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element.
  !!
  elemental subroutine setIdx(self, idx)
    class(element), intent(inout) :: self
    integer(shortInt), intent(in) :: idx
    
    self % idx = idx

  end subroutine setIdx
  
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
  elemental subroutine split(self, faces, edges, vertices, triangles, pyramids, &
                             lastEdgeIdx, lastTriangleIdx, lastPyramidIdx, lastVertexIdx)
    class(element), intent(inout)                    :: self
    type(faceShelf), intent(in)                      :: faces
    type(edgeShelf), intent(inout)                   :: edges
    type(vertexShelf), intent(inout)                 :: vertices
    type(triangleShelf), intent(inout)               :: triangles
    type(pyramidShelf), intent(inout)                :: pyramids
    integer(shortInt), intent(inout)                 :: lastEdgeIdx, lastTriangleIdx, lastPyramidIdx
    integer(shortInt), intent(in)                    :: lastVertexIdx
    integer(shortInt)                                :: i, j, k, faceIdx, absFaceIdx, nVertices, &
                                                        edgeIdx, commonTriangleIdx
    integer(shortInt), dimension(:), allocatable     :: faceVertices
    integer(shortInt), dimension(3)                  :: testVertices
    real(defReal), dimension(3)                      :: apexCoords
    real(defReal), dimension(3, 3)                   :: array
    logical(defBool)                                 :: createNewEdge
    
    ! Retrieve the apex coordinates and loop over all element faces.
    apexCoords = vertices % getVertexCoordinates(lastVertexIdx)
    array(3, :) = apexCoords
    do i = 1, size(self % faceIdxs)
      ! Increment lastPyramidIdx, retrieve current face index and create absolute face index.
      lastPyramidIdx = lastPyramidIdx + 1
      faceIdx = self % faceIdxs(i)
      absFaceIdx = abs(faceIdx)

      ! Create a new pyramid. Set its index, parent element index (the current element's index),
      ! base face index, vertex indices and centroid.
      faceVertices = faces % getFaceVertexIdxs(absFaceIdx)
      call pyramids % initPyramid(lastPyramidIdx, self % idx, faceIdx, [faceVertices, lastVertexIdx], apexCoords, faces)
      
      ! Compute the number of vertices and loop through all the vertices in the current face to 
      ! see if we need to create new triangles.
      nVertices = size(faceVertices)
      do j = 1, nVertices
        ! Check if there are no triangles already containing the same vertices as this one. Create the three 
        ! vertices in the triangle to be tested.
        testVertices = [faceVertices(j), faceVertices(mod(j, nVertices) + 1), lastVertexIdx]

        ! Note: skip the first face because none of the triangles can already be present. 
        if (i > 1) then
          commonTriangleIdx = vertices % findCommonTriangleIdx(testVertices)
          if (commonTriangleIdx > 0) then
            call pyramids % addTriangleIdxToPyramid(lastPyramidIdx, -commonTriangleIdx)
            cycle

          end if

        end if
        
        ! If the triangle is not already present in the shelf the pyramid owns it.
        ! Create a new triangle. Create vectors pointing to the first and second vertices and
        ! compute the new triangle's centre and normal vector.
        lastTriangleIdx = lastTriangleIdx + 1
        array(1, :) = vertices % getVertexCoordinates(testVertices(1))
        array(2, :) = vertices % getVertexCoordinates(testVertices(2))
        call triangles % buildTriangle(lastTriangleIdx, testVertices, 0, .false., array, &
                                       pyramids % getPyramidCentroid(lastPyramidIdx))
        
        ! Update mesh connectivity information.
        call pyramids % addTriangleIdxToPyramid(lastPyramidIdx, lastTriangleIdx)
        do k = 1, 3
          call vertices % addTriangleIdxToVertex(testVertices(k), lastTriangleIdx)

        end do

        ! Find the common edge between the first two vertices of the new triangle (this edge was
        ! already created during the initial importation process) and update mesh connectivity.
        edgeIdx = vertices % findCommonEdgeIdx(testVertices(1), testVertices(2))
        call edges % addTriangleIdxToEdge(edgeIdx, lastTriangleIdx)
        call triangles % addEdgeIdxToTriangle(lastTriangleIdx, edgeIdx)

        ! Now check if new edges need to be created for the remaining two pairs of vertices in the
        ! new triangle.
        do k = 1, 2
          ! Initialise createNew = .true.; if we are on the very first iteration it is guaranteed that
          ! a new edge must be created so skip below check.
          createNewEdge = .true.
          if (i > 1 .or. j > 1) then
            ! Find the common edge between the two vertices and update createNewEdge depending on
            ! whether this common edge exists.
            edgeIdx = vertices % findCommonEdgeIdx(testVertices(k), lastVertexIdx)
            createNewEdge = edgeIdx == 0

          end if

          if (createNewEdge) then
            ! Create a new edge, set its index and vertices then add this edge to each vertex and
            ! update edgeIdx = lastEdgeIdx.
            lastEdgeIdx = lastEdgeIdx + 1
            call edges % initEdge(lastEdgeIdx, [testVertices(k), lastVertexIdx])
            call vertices % addEdgeIdxToVertex(testVertices(k), lastEdgeIdx)
            call vertices % addEdgeIdxToVertex(lastVertexIdx, lastEdgeIdx)
            edgeIdx = lastEdgeIdx

          end if

          ! Update mesh connectivity information.
          call edges % addTriangleIdxToEdge(edgeIdx, lastTriangleIdx)
          call triangles % addEdgeIdxToTriangle(lastTriangleIdx, edgeIdx)

        end do
      
      end do
    
    end do
  
  end subroutine split
  
  !! Subroutine 'testForInclusion'
  !!
  !! Basic description:
  !!   Tests whether a set of 3-D coordinates is inside the element.
  !!
  !! Detailed description:
  !!   First retrieves the faces making the element up. For each face, the subroutine then checks
  !!   whether the dot product between the face's normal vector and a second vector going from the 
  !!   set of 3-D coordinates to the face's centroid is positive. If it is, then the two vectors 
  !!   point in the same direction. If this test is successful for all faces then the coordinates 
  !!   are inside the element.
  !!
  !! Notes: OpenFOAM always numbers a given face's vertices such that the normal vector to this
  !!        face points from the owner element to the neighbour one. Since neighbour elements
  !!        always have greater indices than owner ones, if a given element neighbours a given face
  !!        then the negative of this face's index is added to the 'faces' component of the
  !!        'element' structure. Therefore, in the function below if a face has a negative index,
  !!        its normal vector is flipped.
  !!
  !! Arguments:
  !!   faces [in]            -> A faceShelf.
  !!   r [in]                -> A set of 3-D coordinates.
  !!   failedFace [out]      -> Index of the last face for which the inclusion test fails.
  !!   surfTolFaceIdxs [out] -> An array listing faces for which the dot product is below
  !!                            SURF_TOL, meaning that the coordinates are on the face. It is
  !!                            used in the main tracking routine to assign an element to the
  !!                            coordinates in case the coordinates are on one or more face(s).
  !!
  pure subroutine testForInclusion(self, faces, r, failedFaceIdx, surfTolFaceIdxs)
    class(element), intent(in)                                :: self
    type(faceShelf), intent(in)                               :: faces
    real(defReal), dimension(3), intent(in)                   :: r
    integer(shortInt), intent(out)                            :: failedFaceIdx
    integer(shortInt), dimension(:), allocatable, intent(out) :: surfTolFaceIdxs
    integer(shortInt)                                         :: i, faceIdx, absFaceIdx
    real(defReal), dimension(3)                               :: normal, centroid
    real(defReal)                                             :: dotProduct
    
    ! Initialise failedFace = 0 and loop over all element faces.
    failedFaceIdx = 0
    do i = 1, size(self % faceIdxs)
      ! Create an absolute face index retrieve the face's normal and centroid vectors.
      faceIdx = self % faceIdxs(i)
      absFaceIdx = abs(faceIdx)
      normal = faces % getFaceNormal(faceIdx)
      centroid = faces % getFaceCentroid(absFaceIdx)
      
      ! Make a vector going from the coordinates to the face's centroid and perform the dot
      ! product between this vector and the face's normal vector.
      dotProduct = dot_product(centroid - r, normal)

      ! Check if the coordinates actually lie on the current face, and if so append zeroDotProductFaces.
      if (areEqual(dotProduct, ZERO)) then
        call append(surfTolFaceIdxs, absFaceIdx)
        cycle

      end if

      ! If dotProduct < ZERO, update failedFace and return early.
      if (dotProduct < ZERO) then 
        failedFaceIdx = absFaceIdx
        return

      end if

    end do

  end subroutine testForInclusion

end module element_class