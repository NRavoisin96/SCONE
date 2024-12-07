module element_class
  
  use face_class,          only : face
  use faceShelf_class,     only : faceShelf
  use genericProcedures,   only : append, areEqual, computePyramidCentre, computePyramidVolume, &
                                  computeTetrahedronCentre, computeTetrahedronVolume, &
                                  findCommon, linFind, swap, fatalError, numToChar
  use numPrecision
  use pyramid_class,       only : pyramid
  use pyramidShelf_class,  only : pyramidShelf
  use tetrahedron_class,   only : tetrahedron
  use triangle_class,      only : triangle
  use triangleShelf_class, only : triangleShelf
  use universalVariables,  only : targetNotFound, SURF_TOL, INF, ZERO
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
    integer(shortInt), dimension(:), allocatable :: faceIdxs, vertexIdxs, tetrahedronIdxs
    real(defReal)                                :: volume = ZERO
    real(defReal), dimension(3)                  :: centroid = ZERO
  contains
    ! Build procedures.
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
    procedure                                    :: getFaces
    procedure                                    :: getIdx
    procedure                                    :: getVertices
    procedure                                    :: getVolume
    procedure                                    :: kill
    procedure                                    :: testForInclusion
  end type element

contains
  
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

    ! If array is not already allocated, allocate it and add the index of the vertex.
    ! Return early in this case.
    if (.not. allocated(self % vertexIdxs)) then
      allocate(self % vertexIdxs(1))
      self % vertexIdxs(1) = vertexIdx
      return

    end if
    
    if (linFind(self % vertexIdxs, vertexIdx) == targetNotFound) self % vertexIdxs = [self % vertexIdxs, vertexIdx]

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
    class(element), intent(in)                               :: self
    real(defReal), dimension(3), intent(in)                  :: startPos, endPos
    integer(shortInt), dimension(:), allocatable, intent(in) :: potentialFaceIdxs
    integer(shortInt), intent(out)                           :: intersectedFaceIdx
    real(defReal), intent(out)                               :: lambda
    type(faceShelf), intent(in)                              :: faces
    integer(shortInt)                                        :: i
    integer(shortInt), dimension(:), allocatable             :: absPotentialFaceIdxs
    real(defReal), dimension(3)                              :: centroid, normal
    real(defReal)                                            :: faceLambda
    
    ! Create absolute triangle indices and initialise lambda = INF.
    absPotentialFaceIdxs = abs(potentialFaceIdxs)
    lambda = INF
    
    ! Loop over all potentially intersected triangles.
    do i = 1, size(potentialFaceIdxs)
      ! Retrieve the current triangle's normal vector and flip it if necessary.
      normal = faces % shelf(absPotentialFaceIdxs(i)) % getNormal(potentialFaceIdxs(i))
      
      ! Retrieve the current triangle's centre and compute triangleLambda.
      centroid = faces % shelf(absPotentialFaceIdxs(i)) % getCentroid()
      faceLambda = dot_product(centroid - startPos, normal) / dot_product(endPos - startPos, normal)
      
      ! If triangleLambda < lambda, update lambda and intersectedTriangleIdx.
      if (faceLambda < lambda) then
        lambda = faceLambda
        intersectedFaceIdx = absPotentialFaceIdxs(i)

      end if

    end do

  end subroutine computeIntersectedFace

  !! Subroutine 'computePotentialFaces'
  !!
  !! Basic description:
  !!   Computes a list of indices of the potentially intersected faces using the element's 
  !!   centroid and the end of a line segment.
  !!
  !! Detailed description:
  !!   See Macpherson, et al. (2009). DOI: 10.1002/cnm.1128.
  !!
  !! Arguments:
  !!   endPos [in]             -> End of the line segment.
  !!   faces [in]              -> A faceShelf.
  !!   potentialFaceIdxs [out] -> Array of indices of the potential faces intersected by the line segment.
  !!
  pure subroutine computePotentialFaces(self, endPos, faces, potentialFaceIdxs)
    class(element), intent(in)                                :: self
    real(defReal), dimension(3), intent(in)                   :: endPos
    type(faceShelf), intent(in)                               :: faces
    integer(shortInt), dimension(:), allocatable, intent(out) :: potentialFaceIdxs
    integer(shortInt)                                         :: i, faceIdx, absFaceIdx
    real(defReal)                                             :: lambda
    real(defReal), dimension(3)                               :: centroid, faceCentroid, normal
    
    ! Retrieve element's centroid then loop over all faces in the tetrahedron.
    centroid = self % centroid
    do i = 1, size(self % faceIdxs)
      ! Retrieve current triangle index and create its absolute index.
      faceIdx = self % faceIdxs(i)
      absFaceIdx = abs(faceIdx)
      
      ! Retrieve the signed normal vector of the current triangle.
      normal = faces % shelf(absFaceIdx) % getNormal(faceIdx)
      
      ! Retrieve the centre of the current triangle and compute lambda.
      faceCentroid = faces % shelf(absFaceIdx) % getCentroid()
      lambda = dot_product(faceCentroid - centroid, normal) / dot_product(endPos - centroid, normal)
      
      ! If ZERO <= lambda <= ONE, append the current triangle to the list of potentially intersected triangles.
      if (ZERO <= lambda .and. lambda <= ONE) call append(potentialFaceIdxs, faceIdx)

    end do

    if (.not. allocated(potentialFaceIdxs)) allocate(potentialFaceIdxs(0))

  end subroutine computePotentialFaces
  
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
    class(element), intent(inout)                                :: self
    type(vertex), dimension(size(self % vertexIdxs)), intent(in) :: vertices
    type(face), dimension(size(self % faceIdxs)), intent(in)     :: faces
    integer(shortInt)                                            :: i, nFaces, nVertices
    real(defReal)                                                :: area, sumAreas, volume, &
                                                                    sumVolumes
    real(defReal), dimension(3)                                  :: sumVolumesCentroid, C
    real(defReal), dimension(:, :), allocatable                  :: array
    character(100), parameter                                    :: Here = 'computeVolumeAndCentroid &
                                                                           &(element_class.f90)'
    
    ! Retrieve the number of vertices in the element.
    nVertices = size(vertices)
    
    ! If the element is a tetrahedron perform a direct computation to minimise round-off errors.
    if (nVertices == 4) then
      allocate(array(4, 3))
      do i = 1, 4
        array(i, :) = vertices(i) % getCoordinates()

      end do
      self % volume = computeTetrahedronVolume(array)
      self % centroid = computeTetrahedronCentre(array)

    else
      ! Retrieve the number of faces in the element and initialise variables.
      nFaces = size(faces)
      C = ZERO
      sumAreas = ZERO
      
      ! Loop through all faces and compute the element's approximate centroid
      ! by performing an area-weighted average of the different faces' centroids.
      do i = 1, nFaces
        ! Retrieve the area of the current face and update 
        area = faces(i) % getArea()
        ! Update the area-weighted centroid and the sum of faces' areas.
        C = C + (faces(i) % getCentroid()) * area
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
        array(1, :) = faces(i) % getNormal() * faces(i) % getArea()
        array(2, :) = C - faces(i) % getCentroid()
        volume = computePyramidVolume(array)

        array(1, :) = 3.0_defReal * faces(i) % getCentroid()
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
  !!   isIt -> A logical which is .true. if the element is convex and .false. otherwise.
  !!
  pure function isConvex(self, vertices, faces) result(isIt)
    class(element), intent(in)                                   :: self
    type(vertex), dimension(size(self % vertexIdxs)), intent(in) :: vertices
    type(face), dimension(size(self % faceIdxs)), intent(in)     :: faces
    integer(shortInt)                                            :: i, j, k
    integer(shortInt), dimension(size(self % faceIdxs))          :: faceIdxs
    integer(shortInt), dimension(size(self % vertexIdxs))        :: vertexIdxs
    integer(shortInt), dimension(:), allocatable                 :: faceVertexIdxs
    logical(defBool)                                             :: isIt
    type(vertex), dimension(:), allocatable                      :: faceVertex
    real(defReal), dimension(3)                                  :: normal, faceVertexCoords
    
    ! Initialise isIt = .false.
    isIt = .false.
    
    ! Retrieve the element's vertices and faces and create an array of absolute face indices.
    vertexIdxs = self % vertexIdxs
    faceIdxs = self % faceIdxs
    
    ! Now loop through all the faces in the element.
    do i = 1, size(faceIdxs)
      ! Retrieve the current face's vertices and signed normal vector.
      faceVertexIdxs = faces(i) % getVertices()
      normal = faces(i) % getNormal(faceIdxs(i))
      
      ! Loop through all the vertices in the current face.
      do j = 1, size(faceVertexIdxs)
        ! Retrieve the coordinates of the current face vertex and loop through all the vertices
        ! in the element.
        faceVertex = pack(vertices, vertices % getIdx() == faceVertexIdxs(j))
        faceVertexCoords = faceVertex(1) % getCoordinates()
        do k = 1, size(vertexIdxs)
          ! Cycle to the next vertex if the current vertex index corresponds to the index of a vertex in the current face.
          if (any(faceVertexIdxs == vertexIdxs(k))) cycle
          
          ! Assemble the test vector and check if normal .dot. testVector > ZERO. If yes, the element
          ! is concave and we can return early.
          if (dot_product(normal, vertices(k) % getCoordinates() - faceVertexCoords) > ZERO) return

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
  !!   vertices [in]           -> A vertexShelf.
  !!   triangles [inout]       -> A triangleShelf.
  !!   pyramids [inout]        -> A pyramidShelf.
  !!   freeTriangleIdx [inout] -> Index of the first free item in the triangleShelf.
  !!   freePyramidIdx [inout]  -> Index of the first free item in the pyramidShelf.
  !!   freeVertexIdx [in]      -> Index of the first free item in the vertexShelf.
  !!
  pure subroutine split(self, faces, vertices, triangles, pyramids, freeTriangleIdx, &
                        freePyramidIdx, freeVertexIdx)
    class(element), intent(inout)                    :: self
    type(faceShelf), intent(in)                      :: faces
    type(vertexShelf), intent(inout)                 :: vertices
    type(triangleShelf), intent(inout)               :: triangles
    type(pyramidShelf), intent(inout)                :: pyramids
    integer(shortInt), intent(inout)                 :: freeTriangleIdx, freePyramidIdx
    integer(shortInt), intent(in)                    :: freeVertexIdx
    integer(shortInt)                                :: i, j, k, faceIdx, absFaceIdx, nVertices, &
                                                        initialFreeTriangleIdx
    integer(shortInt), dimension(:), allocatable     :: faceVertices
    integer(shortInt), dimension(3)                  :: testVertices
    real(defReal), dimension(3)                      :: apexCoords, testCoords, firstVertexCoords, &
                                                        secondVertexCoords, tempCoords
    
    ! Initialise initialFreeTriangleIdx = freeTriangleIdx + 1, retrieve the apex coordinates and
    ! loop over all element faces.
    initialFreeTriangleIdx = freeTriangleIdx + 1
    apexCoords = vertices % shelf(freeVertexIdx) % getCoordinates()
    do i = 1, size(self % faceIdxs)
      ! Increment freePyramidIdx, retrieve current face index and create absolute face index.
      freePyramidIdx = freePyramidIdx + 1
      faceIdx = self % faceIdxs(i)
      absFaceIdx = abs(faceIdx)

      ! Create a new pyramid. Set its index, parent element index (the current element's index),
      ! base face index, vertex indices and centroid.
      faceVertices = faces % shelf(absFaceIdx) % getVertices()
      call pyramids % shelf(freePyramidIdx) % setIdx(freePyramidIdx)
      call pyramids % shelf(freePyramidIdx) % setElement(self % idx)
      call pyramids % shelf(freePyramidIdx) % setFace(faceIdx)
      call pyramids % shelf(freePyramidIdx) % setVertices([faceVertices, freeVertexIdx])
      call pyramids % shelf(freePyramidIdx) % computeCentroid(faces % shelf(absFaceIdx), apexCoords)
      
      ! Compute the number of vertices and loop through all the vertices in the current face to 
      ! see if we need to create new triangles.
      nVertices = size(faceVertices)
      vertexLoop: do j = 1, nVertices
        ! Check if there are no triangles already containing the same vertices as this one. Create the three 
        ! vertices in the triangle to be tested.
        testVertices = [faceVertices(j), faceVertices(mod(j, nVertices) + 1), freeVertexIdx]

        ! Note: skip the first face because none of the triangles can already be present. 
        if (i > 1) then
          
          ! Loop through all generated triangles. TODO: Import edges and refactor this to avoid long loop.
          do k = initialFreeTriangleIdx, freeTriangleIdx
            ! If the triangle already exists then the current pyramid is its neighbour. In this case
            ! we simply update mesh connectivity information and we move onto the next vertex.
            if (size(findCommon(triangles % shelf(k) % getVertices(), testVertices)) == 3) then
                call pyramids % shelf(freePyramidIdx) % addTriangle(-triangles % shelf(k) % getIdx())
                cycle vertexLoop

            end if

          end do

        end if
        
        ! If the triangle is not already present in the shelf the pyramid owns it.
        ! Create a new triangle. Create vectors pointing to the first and second vertices and
        ! compute the new triangle's centre and normal vector.
        freeTriangleIdx = freeTriangleIdx + 1
        call triangles % shelf(freeTriangleIdx) % setIdx(freeTriangleIdx)
        firstVertexCoords = vertices % shelf(testVertices(1)) % getCoordinates()
        secondVertexCoords = vertices % shelf(testVertices(2)) % getCoordinates()
        call triangles % shelf(freeTriangleIdx) % computeCentre(firstVertexCoords, &
                                                                secondVertexCoords, apexCoords)
        call triangles % shelf(freeTriangleIdx) % computeNormal(firstVertexCoords, &
                                                                secondVertexCoords, apexCoords)
        
        ! Now compute the triangle's area and normalise its normal vector. Check that the normal
        ! vector points away from the current face. If not, flip two of its vertices and negate it.
        call triangles % shelf(freeTriangleIdx) % computeArea()
        call triangles % shelf(freeTriangleIdx) % normaliseNormal()
        testCoords = triangles % shelf(freeTriangleIdx) % getCentre() - pyramids % shelf(freePyramidIdx) % getCentroid()
        
        if (dot_product(testCoords, triangles % shelf(freeTriangleIdx) % getNormal()) <= ZERO) then
          call swap(testVertices, 1, 2)
          tempCoords = firstVertexCoords
          firstVertexCoords = secondVertexCoords
          secondVertexCoords = tempCoords
          call triangles % shelf(freeTriangleIdx) % negateNormal()

        end if
        
        ! Set the new triangle's vertices and edge vectors.
        call triangles % shelf(freeTriangleIdx) % setVertices(testVertices)
        
        ! Update mesh connectivity information.
        call pyramids % shelf(freePyramidIdx) % addTriangle(freeTriangleIdx)
        do k = 1, 3
          call vertices % shelf(testVertices(k)) % addTriangleIdx(freeTriangleIdx)

        end do
      
      end do vertexLoop
    
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
  pure subroutine testForInclusion(self, faces, r, failedFace, surfTolFaceIdxs)
    class(element), intent(in)                                :: self
    type(faceShelf), intent(in)                               :: faces
    real(defReal), dimension(3), intent(in)                   :: r
    integer(shortInt), intent(out)                            :: failedFace
    integer(shortInt), dimension(:), allocatable, intent(out) :: surfTolFaceIdxs
    integer(shortInt)                                         :: i, faceIdx, absFaceIdx
    real(defReal), dimension(3)                               :: normal, centroid
    real(defReal)                                             :: dotProduct
    
    ! Initialise failedFace = 0 and loop over all element faces.
    failedFace = 0
    do i = 1, size(self % faceIdxs)
      ! Create an absolute face index retrieve the face's normal and centroid vectors.
      faceIdx = self % faceIdxs(i)
      absFaceIdx = abs(faceIdx)
      normal = faces % shelf(absFaceIdx) % getNormal(faceIdx)
      centroid = faces % shelf(absFaceIdx) % getCentroid()
      
      ! Make a vector going from the coordinates to the face's centroid and perform the dot
      ! product between this vector and the face's normal vector.
      dotProduct = dot_product(centroid - r, normal)

      ! Check if the coordinates actually lie on the current face, and if so append zeroDotProductFaces.
      if (areEqual(dotProduct, ZERO)) then
        call append(surfTolFaceIdxs, faceIdx)
        cycle

      end if

      ! If dotProduct < ZERO, update failedFace and return early.
      if (dotProduct < ZERO) then 
        failedFace = absFaceIdx
        return

      end if

    end do

  end subroutine testForInclusion

end module element_class