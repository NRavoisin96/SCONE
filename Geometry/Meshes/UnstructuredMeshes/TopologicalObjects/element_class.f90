module element_class

  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use coord_class,                  only : coord
  use edge_class,                   only : edgeBox
  use face_class,                   only : faceBox, orientatedFaceBox
  use genericProcedures,            only : append, areEqual, crossProduct, findCommon, fatalError, numToChar
  use numPrecision
  use publicObjects,                only : basicElementInfo
  use topologicalObject_inter,      only : topologicalObject, kill_super => kill
  use universalVariables,           only : FOURTH, INSIDE_ELEMENT, INF, ON_BOUNDARY_ELEMENT, OUTSIDE_ELEMENT, &
                                           SIXTH, SURF_TOL, ZERO
  use vertex_class,                 only : vertexBox
  
  implicit none
  private

  !!
  !!
  !!
  type, public :: buildElementInfo
    integer(shortInt)                                  :: idx = 0, localId = 0, parentIdx = 0
    type(edgeBox), dimension(:), allocatable           :: edges
    type(orientatedFaceBox), dimension(:), allocatable :: orientatedFaces
    type(vertexBox), dimension(:), allocatable         :: vertices
  end type buildElementInfo

  !!
  !! Small, local container to store elements in a single array.
  !!
  !! Public members:
  !!   name -> Name of the mesh.
  !!   ptr  -> Pointer to the mesh.
  !!
  type, public             :: elementBox
    type(element), pointer :: ptr => null()
  end type
  
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
  type, public, extends(topologicalObject)             :: element
    private
    integer(shortInt)                                  :: parentIdx = 0, localId = 0
    type(edgeBox), dimension(:), allocatable           :: edges
    type(orientatedFaceBox), dimension(:), allocatable :: orientatedFaces
    type(vertexBox), dimension(:), allocatable         :: vertices
    integer(shortInt), dimension(:), allocatable       :: tetrahedronIdxs
    real(defReal)                                      :: volume = ZERO
    real(defReal), dimension(3)                        :: centroid = ZERO
    type(axisAlignedBoundingBox)                       :: boundingBox
    logical(defBool)                                   :: isActive = .true., isConvex = .false.
    character(:), allocatable                          :: type
  contains
    ! Build procedures.
    procedure :: addEdge
    procedure :: addFace
    procedure :: addVertex
    procedure :: build
    procedure :: computeConvexity
    procedure :: init
    procedure :: setLocalId
    ! Runtime procedures.
    procedure :: computeIntersectedFace
    procedure :: computePotentialFaces
    procedure :: deactivate
    procedure :: getBoundingBox
    procedure :: getCentroid
    procedure :: getEdges
    procedure :: getIsActive
    procedure :: getIsConvex
    procedure :: getLocalId
    procedure :: getOrientatedFaces
    procedure :: getParentIdx
    procedure :: getType
    procedure :: getVertices
    procedure :: getVolume
    procedure :: kill
    procedure :: pushFromBoundary
    procedure :: isPointInside
  end type element

  !!
  !!
  !!
  type, public        :: inclusionTestResult
    integer(shortInt) :: status = INSIDE_ELEMENT, failedFaceIdx = 0
  end type inclusionTestResult

contains

  !! Subroutine 'addEdgeIdx'
  !!
  !! Basic description:
  !!   Adds the index of an edge sharing the element.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the edge.
  !!
  subroutine addEdge(self, edge)
    class(element), intent(inout)            :: self
    type(edgeBox), intent(in)                :: edge
    integer(shortInt)                        :: nEdges
    type(edgeBox), dimension(:), allocatable :: tempEdges

    if (allocated(self % edges)) then
      nEdges = size(self % edges)
      allocate(tempEdges(nEdges + 1))
      tempEdges(1:nEdges) = self % edges
      tempEdges(nEdges + 1) = edge
      call move_alloc(tempEdges, self % edges)

    else
      allocate(self % edges(1))
      self % edges(1) = edge

    end if

  end subroutine addEdge
  
  !! Subroutine 'addFaceToElement'
  !!
  !! Basic description:
  !!   Adds the index of a face belonging to the element.
  !!
  !! Arguments:
  !!   faceIdx [in] -> Index of the face.
  !!
  subroutine addFace(self, orientatedFace)
    class(element), intent(inout)                      :: self
    type(orientatedFaceBox), intent(in)                :: orientatedFace
    integer(shortInt)                                  :: nFaces
    type(orientatedFaceBox), dimension(:), allocatable :: tempOrientatedFaces
    
    if (allocated(self % orientatedFaces)) then
      nFaces = size(self % orientatedFaces)
      allocate(tempOrientatedFaces(nFaces + 1))
      tempOrientatedFaces(1:nFaces) = self % orientatedFaces
      tempOrientatedFaces(nFaces + 1) = orientatedFace
      call move_alloc(tempOrientatedFaces, self % orientatedFaces)

    else
      allocate(self % orientatedFaces(1))
      self % orientatedFaces(1) = orientatedFace

    end if

  end subroutine addFace
  
  !! Subroutine 'addVertexToElement'
  !!
  !! Basic description:
  !!   Adds the index of a vertex belonging to the element. Only adds it if the index is not already
  !!   present.
  !!
  !! Arguments:
  !!   vertexIdx [in] -> Index of the vertex.
  !!
  subroutine addVertex(self, vertex)
    class(element), intent(inout)              :: self
    type(vertexBox), intent(in)                :: vertex
    integer(shortInt)                          :: nVertices
    type(vertexBox), dimension(:), allocatable :: tempVertices

    if (allocated(self % vertices)) then
      nVertices = size(self % vertices)
      allocate(tempVertices(nVertices + 1))
      tempVertices(1:nVertices) = self % vertices
      tempVertices(nVertices + 1) = vertex
      call move_alloc(tempVertices, self % vertices)

    else
      allocate(self % vertices(1))
      self % vertices(1) = vertex

    end if

  end subroutine addVertex

  !!
  !!
  !!
  subroutine build(self)
    class(element), intent(inout)                      :: self
    integer(shortInt)                                  :: i, nFaces, nVertices
    real(defReal), dimension(3, size(self % vertices)) :: allCoords
    real(defReal)                                      :: faceArea, pyramidVolume, sumVolumes
    real(defReal), dimension(3)                        :: outwardNormal, faceCentroid, geometricCentroid, sumVolumesCentroid
    character(*), parameter                            :: here = 'build (element_class.f90)'

    ! Compute the number of vertices in the element.
    nVertices = size(self % vertices)
    do i = 1, nVertices
      if (.not. associated(self % vertices(i) % ptr)) call fatalError(here, 'Element contains a null vertex pointer.')
      allCoords(:, i) = self % vertices(i) % ptr % getCoordinates()

    end do

    ! If the element is a tetrahedron, perform a direct calculation to avoid round-off errors.
    if (nVertices == 4) then
      self % isConvex = .true.
      self % centroid = FOURTH * sum(allCoords, 2)
      self % volume = SIXTH * abs(dot_product(crossProduct(allCoords(:, 2) - allCoords(:, 1), allCoords(:, 3) - allCoords(:, 1)), &
                                              allCoords(:, 4) - allCoords(:, 1)))

    else
      ! Check if current element is convex and call fatalError if not.
      call self % computeConvexity()
      if (.not. self % isConvex) call fatalError(here, 'Element with index: '//numToChar(self % getIdx())//' is concave.')

      ! Approximate the centroid by taking the arithmetic average of all the vertices in the polyhedron.
      geometricCentroid = sum(allCoords, 2) / nVertices
      
      nFaces = size(self % orientatedFaces)
      sumVolumes = ZERO
      sumVolumesCentroid = ZERO
      
      ! Loop through all faces (pyramids).
      do i = 1, nFaces
        ! Retrieve the volume of the current pyramid and update the volume-weighted centroid and the sum of volumes.
        faceArea = self % orientatedFaces(i) % face % ptr % getArea()
        faceCentroid = self % orientatedFaces(i) % face % ptr % getCentroid()
        outwardNormal = self % orientatedFaces(i) % outwardNormal
        
        pyramidVolume = THIRD * abs(dot_product(faceCentroid - geometricCentroid, outwardNormal * faceArea))
        sumVolumes = sumVolumes + pyramidVolume
        sumVolumesCentroid = sumVolumesCentroid + FOURTH * (3.0_defReal * faceCentroid + geometricCentroid) * pyramidVolume

      end do
      ! The volume of the element is simply the sum of volumes, while the centroid is the average of
      ! the volume-weighted sum.
      self % volume = sumVolumes
      self % centroid = sumVolumesCentroid / sumVolumes

    end if

    ! Initialise bounding box.
    call self % boundingBox % computeBounds(allCoords)

  end subroutine build

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
  subroutine computeConvexity(self)
    class(element), intent(inout)              :: self
    logical(defBool)                           :: isOnFace
    integer(shortInt)                          :: i, j, k
    type(vertexBox), dimension(:), allocatable :: faceVertices
    real(defReal), dimension(3)                :: faceVertexCoords, outwardNormal

    ! Initialise isIt = .false.
    self % isConvex = .false.
    
    ! Now loop through all the faces in the element.
    do i = 1, size(self % orientatedFaces)
      ! Retrieve the current face's vertices and signed normal vector.
      faceVertices = self % orientatedFaces(i) % face % ptr % getVertices()
      faceVertexCoords = faceVertices(1) % ptr % getCoordinates()
      outwardNormal = self % orientatedFaces(i) % outwardNormal

      ! Loop through all vertices in the element.
      do j = 1, size(self % vertices)
        isOnFace = .false.
        do k = 1, size(faceVertices)
          if (.not. associated(faceVertices(k) % ptr)) cycle
          if (associated(self % vertices(j) % ptr, faceVertices(k) % ptr)) then
            isOnFace = .true.
            exit

          end if

        end do

        if (isOnFace) cycle
        if (dot_product(outwardNormal, self % vertices(j) % ptr % getCoordinates() - faceVertexCoords) > ZERO) return

      end do

    end do
    
    ! If reached this point the element is convex. Update isIt = .true.
    self % isConvex = .true.

  end subroutine computeConvexity

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
  subroutine computeIntersectedFace(self, r, rEnd, potentialFaceIdxs, intersectedFaceIdx, lambda)
    class(element), intent(in)                  :: self
    real(defReal), dimension(3), intent(in)     :: r, rEnd
    integer(shortInt), dimension(:), intent(in) :: potentialFaceIdxs
    integer(shortInt), intent(out)              :: intersectedFaceIdx
    real(defReal), intent(out)                  :: lambda
    integer(shortInt)                           :: i
    real(defReal), dimension(3)                 :: normal
    real(defReal)                               :: faceLambda
    
    ! Initialise lambda = INF.
    lambda = INF
    
    ! Loop over all potentially intersected triangles.
    do i = 1, size(potentialFaceIdxs)
      ! Retrieve the current face's signed normal vector and flip it if necessary.
      normal = self % orientatedFaces(potentialFaceIdxs(i)) % outwardNormal
      
      ! Compute lambda.
      faceLambda = dot_product(self % orientatedFaces(potentialFaceIdxs(i)) % face % ptr % getCentroid() - r, normal) / &
                   dot_product(rEnd - r, normal)
      
      ! If triangleLambda < lambda, update lambda and intersectedTriangleIdx.
      if (faceLambda < lambda) then
        lambda = faceLambda
        intersectedFaceIdx = i

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
  function computePotentialFaces(self, rEnd) result(potentialFaceIdxs)
    class(element), intent(in)                   :: self
    real(defReal), dimension(3), intent(in)      :: rEnd
    integer(shortInt), dimension(:), allocatable :: potentialFaceIdxs
    integer(shortInt)                            :: i
    real(defReal)                                :: lambda
    real(defReal), dimension(3)                  :: centroid, faceCentroid, outwardNormal
    
    ! Retrieve element's centroid then loop over all faces in the tetrahedron.
    allocate(potentialFaceIdxs(0))
    centroid = self % centroid
    do i = 1, size(self % orientatedFaces)
      ! Retrieve the signed normal vector of the current face.
      outwardNormal = self % orientatedFaces(i) % outwardNormal
      faceCentroid = self % orientatedFaces(i) % face % ptr % getCentroid()
      
      ! Retrieve the centre of the current face and compute lambda.
      lambda = dot_product(faceCentroid - centroid, outwardNormal) / dot_product(rEnd - centroid, outwardNormal)
      
      ! If ZERO <= lambda <= ONE, append the current face to the list of potentially intersected faces.
      if (ZERO <= lambda .and. lambda <= ONE) call append(potentialFaceIdxs, i)

    end do

  end function computePotentialFaces

  !!
  !!
  !!
  elemental subroutine deactivate(self)
    class(element), intent(inout) :: self

    self % isActive = .false.

  end subroutine deactivate

  !! Function 'getBoundingBox'
  !!
  !! Basic description:
  !!   Returns the bounding box of the element.
  !!
  !! Result:
  !!   boundingBox -> 6-D array representing the bounding box of the element.
  !!
  pure function getBoundingBox(self) result(boundingBox)
    class(element), intent(in)   :: self
    type(axisAlignedBoundingBox) :: boundingBox

    boundingBox = self % boundingBox

  end function getBoundingBox
  
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
  function getEdges(self) result(edges)
    class(element), intent(in)                   :: self
    type(edgeBox), dimension(size(self % edges)) :: edges

    edges = self % edges

  end function getEdges

  !!
  !!
  !!
  elemental function getIsActive(self) result(isActive)
    class(element), intent(in) :: self
    logical(defBool)           :: isActive

    isActive = self % isActive

  end function getIsActive

  !!
  !!
  !!
  elemental function getIsConvex(self) result(isConvex)
    class(element), intent(in) :: self
    logical(defBool)           :: isConvex

    isConvex = self % isConvex

  end function getIsConvex

  !!
  !!
  !!
  elemental function getLocalId(self) result(localId)
    class(element), intent(in) :: self
    integer(shortInt)          :: localId

    localId = self % localId

  end function getLocalId

  !! Function 'getFaces'
  !!
  !! Basic description:
  !!   Returns the indices of the faces in the element.
  !!
  !! Result:
  !!   faceIdxs -> An array listing the indices of the faces in the element.
  !!
  function getOrientatedFaces(self) result(orientatedFaces)
    class(element), intent(in)                         :: self
    type(orientatedFaceBox), dimension(:), allocatable :: orientatedFaces
    
    if (allocated(self % orientatedFaces)) then
      orientatedFaces = self % orientatedFaces

    else
      allocate(orientatedFaces(0))

    end if

  end function getOrientatedFaces

  !! Function 'getParentIdx'
  !!
  !! Basic description:
  !!   Returns the index of the parent element of the element.
  !!
  !! Result:
  !!   parentIdx -> Index of the parent element of the element.
  !!
  elemental function getParentIdx(self) result(parentIdx)
    class(element), intent(in) :: self
    integer(shortInt)          :: parentIdx

    parentIdx = self % parentIdx

  end function getParentIdx

  !!
  !!
  !!
  pure function getType(self) result(type)
    class(element), intent(in) :: self
    character(:), allocatable  :: type

    type = self % type

  end function  getType
  
  !! Function 'getVertices'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in the element.
  !!
  !! Result:
  !!   vertexIdxs -> An array listing indices of the vertices in the element.
  !!
  function getVertices(self) result(vertices)
    class(element), intent(in)                        :: self
    type(vertexBox), dimension(size(self % vertices)) :: vertices
    
    vertices = self % vertices

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

  !!
  !!
  !!
  subroutine init(self, info)
    class(element), intent(inout)      :: self
    type(buildElementInfo), intent(in) :: info
    integer(shortInt)                  :: nFaces, nVertices
    character(*), parameter            :: here = 'init (element_class.f90)'

    ! Catch invalid number of vertices and faces.
    nVertices = size(info % vertices)
    if (nVertices < 4) call fatalError(here, 'An element must have at least 4 vertices. Has: '//numToChar(nVertices)//'.')

    nFaces = size(info % orientatedFaces)
    if (nFaces < 4) call fatalError(here, 'An element must have at least 4 faces. Has: '//numToChar(nFaces)//'.')
    
    if (nVertices == 4) then
      self % type = 'Tetrahedron'

    else
      self % type = 'Polyhedron'

    end if
    
    ! Set everything from payload.
    call self % setIdx(info % idx)
    self % localId = info % localId
    self % parentIdx = info % parentIdx
    self % orientatedFaces = info % orientatedFaces
    self % vertices = info % vertices
    self % edges = info % edges

    ! Build.
    call self % build()

  end subroutine init

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
  function isPointInside(self, r) result(result)
    class(element), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    type(inclusionTestResult)               :: result
    integer(shortInt)                       :: i
    real(defReal)                           :: dotProduct
    logical(defBool)                        :: isOnBoundary
    
    ! Initialise isOnBoundary = .false. and result % status = INSIDE_ELEMENT then loop over all element faces.
    isOnBoundary = .false.
    result % status = INSIDE_ELEMENT
    do i = 1, size(self % orientatedFaces)
      ! Make a vector going from the coordinates to the face's centroid and perform the dot
      ! product between this vector and the face's normal vector.
      dotProduct = dot_product(self % orientatedFaces(i) % face % ptr % getCentroid() - r, &
                               self % orientatedFaces(i) % outwardNormal)

      ! Check if the point is effectively on the plane of this face.
      if (areEqual(dotProduct, ZERO)) then
        isOnBoundary = .true.

        ! Store the first face found and cycle to search other faces.
        if (result % failedFaceIdx == 0) result % failedFaceIdx = self % orientatedFaces(i) % face % ptr % getIdx()
        cycle

      end if

      ! If dotProduct < ZERO, update result and return early.
      if (dotProduct < ZERO) then
        result % status = OUTSIDE_ELEMENT
        result % failedFaceIdx = self % orientatedFaces(i) % face % ptr % getIdx()
        return

      end if

    end do

    ! If point is on boundary, update result % status.
    if (isOnBoundary) result % status = ON_BOUNDARY_ELEMENT

  end function isPointInside
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(element), intent(inout) :: self
    integer(shortInt)             :: i
    
    ! Superclass.
    call kill_super(self)

    ! Local.
    self % parentIdx = 0
    self % localId = 0
    self % volume = ZERO
    self % centroid = ZERO
    self % isActive = .true.
    self % isConvex = .false.
    if (allocated(self % tetrahedronIdxs)) deallocate(self % tetrahedronIdxs)
    if (allocated(self % type)) deallocate(self % type)
    call self % boundingBox % kill()

    if (allocated(self % edges)) then
      do i = 1, size(self % edges)
        nullify(self % edges(i) % ptr)

      end do
      deallocate(self % edges)

    end if

    if (allocated(self % orientatedFaces)) then
      do i = 1, size(self % orientatedFaces)
        nullify(self % orientatedFaces(i) % face % ptr)
        self % orientatedFaces(i) % isOwner = .false.
        self % orientatedFaces(i) % outwardNormal = ZERO

      end do
      deallocate(self % orientatedFaces)

    end if

    if (allocated(self % vertices)) then
      do i = 1, size(self % vertices)
        nullify(self % vertices(i) % ptr)

      end do
      deallocate(self % vertices)

    end if

  end subroutine kill

  !!
  !!
  !!
  subroutine pushFromBoundary(self, coords)
    class(element), intent(in)   :: self
    type(coord), intent(inout)   :: coords
    real(defReal), dimension(3)  :: nudgeDirection, outwardNormal, r, u
    integer(shortInt)            :: i

    ! Initialise nudgeDirection = ZERO then loop over all the faces in the element.
    nudgeDirection = ZERO
    r = coords % getPositionToNudge()
    u = coords % getDirection()
    do i = 1, size(self % orientatedFaces)
      ! Retrieve the normal vector of the current face and test whether the coordinates lie on the face.
      outwardNormal = self % orientatedFaces(i) % outwardNormal
      if (areEqual(dot_product(self % orientatedFaces(i) % face % ptr % getCentroid() - r, outwardNormal), ZERO)) then
        ! If coordinates are parallel to the plane of the current face, append the negative of the normal to
        ! nudgeDirection.
        if (areEqual(dot_product(u, outwardNormal), ZERO)) nudgeDirection = nudgeDirection - outwardNormal

      end if

    end do

    ! Now nudge coordinates with the appropriate direction.
    if (any(nudgeDirection /= ZERO)) then
      call coords % nudgePosition(nudgeDirection / norm2(nudgeDirection))

    else
      call coords % nudgePosition()

    end if

  end subroutine pushFromBoundary

  !!
  !!
  !!
  elemental subroutine setLocalId(self, localId)
    class(element), intent(inout) :: self
    integer(shortInt), intent(in) :: localId

    self % localId = localId

  end subroutine setLocalId

end module element_class