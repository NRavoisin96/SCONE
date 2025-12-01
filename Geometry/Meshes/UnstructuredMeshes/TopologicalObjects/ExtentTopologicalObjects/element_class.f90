module element_class

  use axisAlignedBoundingBox_class,  only : axisAlignedBoundingBox
  use extentTopologicalObject_inter, only : buildExtentTopologicalObjectPayload, extentTopologicalObject, &
                                            intersects_Ray_super => intersects_Ray
  use edge_class,                    only : edgeBox
  use face_class,                    only : faceBox, orientatedFaceBox
  use genericProcedures,             only : append, areEqual, crossProduct, findCommon, fatalError, numToChar
  use numPrecision
  use publicObjects,                 only : basicElementInfo, intersectionTestPayload, intersectionTestResult, &
                                            resetIntersectionTestResult
  use RNG_class,                     only : RNG
  use topologicalObject_inter,       only : buildTopologicalObjectPayload, kill_super => kill, topologicalObjectBox
  use universalVariables,            only : FOURTH, INSIDE_ELEMENT, INF, NUDGE, ON_BOUNDARY_ELEMENT, ONE, OUTSIDE_ELEMENT, &
                                            SIXTH, ZERO
  use vertex_class,                  only : vertexBox
  
  implicit none
  private

  ! Public procedures.
  public :: newElementIntersectionTestPayload, resetElementIntersectionTestResult

  !!
  !!
  !!
  type, public, extends(buildExtentTopologicalObjectPayload) :: buildElementPayload
    integer(shortInt)                                        :: localId = 0, parentIdx = 0
    type(edgeBox), dimension(:), allocatable                 :: edges
    type(orientatedFaceBox), dimension(:), allocatable       :: orientatedFaces
  end type buildElementPayload

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
  type, public, extends(extentTopologicalObject)       :: element
    private
    integer(shortInt)                                  :: parentIdx = 0, localId = 0
    type(edgeBox), dimension(:), allocatable           :: edges
    type(orientatedFaceBox), dimension(:), allocatable :: orientatedFaces
    type(vertexBox), dimension(:), allocatable         :: vertices
    integer(shortInt), dimension(:), allocatable       :: childrenIdxs
    real(defReal)                                      :: volume = ZERO
    logical(defBool)                                   :: isConvex = .false.
    character(:), allocatable                          :: type
  contains
    ! Build procedures.
    procedure          :: addChildIdx
    procedure          :: addEdge
    procedure          :: addFace
    procedure          :: addVertex
    procedure          :: build
    procedure, private :: buildComponents
    procedure          :: computeConvexity
    procedure          :: connectComponents
    procedure          :: setLocalId
    ! Runtime procedures.
    procedure          :: distanceSquared
    procedure          :: getChildrenIdxs
    procedure          :: getEdges
    procedure          :: getSharingElements
    procedure          :: getIsConvex
    procedure          :: getLocalId
    procedure          :: getOrientatedFaces
    procedure          :: getParentIdx
    procedure          :: getType
    procedure          :: getVertices
    procedure          :: getVolume
    procedure          :: intersects_BoundingBox
    procedure          :: intersects_Ray
    procedure          :: isPointInside
    procedure          :: kill
    procedure          :: minimumDistance
    procedure          :: pushFromBoundary
    procedure          :: sampleInitialPosition
  end type element

  !!
  !!
  !!
  type, public        :: inclusionTestResult
    integer(shortInt) :: status = INSIDE_ELEMENT, failedFaceIdx = 0
  end type inclusionTestResult

  !!
  !!
  !!
  type, public, extends(intersectionTestPayload) :: elementIntersectionTestPayload
    logical(defBool)                             :: excludeZeroFaces = .false., skipBoundingBoxIntersectionTest = .false.
  end type elementIntersectionTestPayload

  !!
  !!
  !!
  type, public, extends(intersectionTestResult) :: elementIntersectionTestResult
    type(faceBox)                               :: intersectedFace
  end type elementIntersectionTestResult

contains
  !!
  !!
  !!
  subroutine addChildIdx(self, childIdx)
    class(element), intent(inout)                :: self
    integer(shortInt), intent(in)                :: childIdx
    integer(shortInt)                            :: nChildren
    integer(shortInt), dimension(:), allocatable :: tempChildrenIdxs

    if (allocated(self % childrenIdxs)) then
      nChildren = size(self % childrenIdxs)
      allocate(tempChildrenIdxs(nChildren + 1))
      tempChildrenIdxs(1:nChildren) = self % childrenIdxs
      tempChildrenIdxs(nChildren + 1) = childIdx
      call move_alloc(tempChildrenIdxs, self % childrenIdxs)

    else
      allocate(self % childrenIdxs(1))
      self % childrenIdxs(1) = childIdx

    end if

  end subroutine addChildIdx

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
  subroutine build(self, payload)
    class(element), intent(inout)                       :: self
    class(buildTopologicalObjectPayload), intent(inout) :: payload
    type(buildElementPayload), pointer                  :: payloadPtr
    integer(shortInt)                                   :: nFaces, nVertices
    character(*), parameter                             :: here = 'build (element_class.f90)'

    ! Downcast payload to correct type.
    select type(ptr => payload)
      type is(buildElementPayload)
        payloadPtr => ptr

      class default
        call fatalError(here, 'Invalid payload type.')

    end select

    ! Catch invalid number of vertices and faces.
    nVertices = size(payloadPtr % vertices)
    if (nVertices < 4) call fatalError(here, 'An element must have at least 4 vertices. Has: '//numToChar(nVertices)//'.')

    nFaces = size(payloadPtr % orientatedFaces)
    if (nFaces < 4) call fatalError(here, 'An element must have at least 4 faces. Has: '//numToChar(nFaces)//'.')
    
    if (nVertices == 4) then
      self % type = 'Tetrahedron'

    else
      self % type = 'Polyhedron'

    end if
    
    ! Set everything from payload.
    self % localId = payloadPtr % localId
    self % parentIdx = payloadPtr % parentIdx
    self % orientatedFaces = payloadPtr % orientatedFaces
    self % vertices = payloadPtr % vertices
    self % edges = payloadPtr % edges

    ! Build.
    call self % buildComponents(payloadPtr)

  end subroutine build

  !!
  !!
  !!
  subroutine buildComponents(self, payload)
    class(element), intent(inout)            :: self
    type(buildElementPayload), intent(inout) :: payload
    integer(shortInt)                        :: i, nFaces, nVertices
    real(defReal)                            :: faceArea, pyramidVolume, sumVolumes
    real(defReal), dimension(3)              :: outwardNormal, faceCentroid, geometricCentroid, sumVolumesCentroid
    character(*), parameter                  :: here = 'buildComponents (element_class.f90)'

    ! Compute the number of vertices in the element.
    nVertices = size(self % vertices)
    allocate(payload % allCoords(3, nVertices))
    do i = 1, nVertices
      if (.not. associated(self % vertices(i) % ptr)) call fatalError(here, 'Element contains a null vertex pointer.')
      payload % allCoords(:, i) = self % vertices(i) % ptr % getCoordinates()

    end do

    ! If the element is a tetrahedron, perform a direct calculation to avoid round-off errors.
    if (nVertices == 4) then
      self % isConvex = .true.
      payload % centroid = FOURTH * sum(payload % allCoords, 2)
      self % volume = SIXTH * abs(dot_product(crossProduct(payload % allCoords(:, 2) - payload % allCoords(:, 1), &
                                                           payload % allCoords(:, 3) - payload % allCoords(:, 1)), &
                                              payload % allCoords(:, 4) - payload % allCoords(:, 1)))

    else
      ! Check if current element is convex and call fatalError if not.
      call self % computeConvexity()
      if (.not. self % isConvex) call fatalError(here, 'Element with index: '//numToChar(self % getIdx())//' is concave.')

      ! Approximate the centroid by taking the arithmetic average of all the vertices in the polyhedron.
      geometricCentroid = sum(payload % allCoords, 2) / nVertices
      
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
      payload % centroid = sumVolumesCentroid / sumVolumes

    end if

  end subroutine buildComponents

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

  !!
  !!
  !!
  subroutine connectComponents(self)
    class(element), target, intent(inout) :: self
    type(topologicalObjectBox)            :: box
    integer(shortInt)                     :: i

    box % ptr => self
    do i = 1, size(self % orientatedFaces)
      call self % orientatedFaces(i) % face % ptr % addSharingElement(box)

    end do

    do i = 1, size(self % edges)
      call self % edges(i) % ptr % addSharingElement(box)

    end do

    do i = 1, size(self % vertices)
      call self % vertices(i) % ptr % addSharingElement(box)

    end do

  end subroutine connectComponents

  !!
  !!
  !!
  function distanceSquared(self, r) result(dSquared)
    class(element), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: dSquared
    integer(shortInt)                       :: i

    dSquared = INF
    do i = 1, size(self % orientatedFaces)
      dSquared = min(dSquared, self % orientatedFaces(i) % face % ptr % distanceSquared(r))

    end do

  end function distanceSquared

  !!
  !!
  !!
  pure function getChildrenIdxs(self) result(childrenIdxs)
    class(element), intent(in)                   :: self
    integer(shortInt), dimension(:), allocatable :: childrenIdxs

    if (allocated(self % childrenIdxs)) then
      childrenIdxs = self % childrenIdxs

    else
      allocate(childrenIdxs(0))

    end if

  end function getChildrenIdxs

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
  function getSharingElements(self) result(sharingElements)
    class(element), target, intent(in)                    :: self
    type(topologicalObjectBox), dimension(:), allocatable :: sharingElements

    allocate(sharingElements(1))
    sharingElements(1) % ptr => self

  end function getSharingElements

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
  pure subroutine intersects_BoundingBox(self, boundingBox, doesIt)
    class(element), intent(in)               :: self
    type(axisAlignedBoundingBox), intent(in) :: boundingBox
    logical(defBool), intent(out)            :: doesIt
    integer(shortInt)                        :: i

    ! Initialise doesIt = .false.
    doesIt = .false.
    if (.not. self % intersectsBoundingBox(boundingBox)) return

    ! Loop over all faces in the element and check for intersection with any of them.
    do i = 1, size(self % orientatedFaces)
      call self % orientatedFaces(i) % face % ptr % intersects(boundingBox, doesIt)
      if (doesIt) return

    end do

  end subroutine intersects_BoundingBox

  !!
  !!
  !!
  subroutine intersects_Ray(self, payload, result)
    class(element), intent(in)                    :: self
    class(intersectionTestPayload), intent(in)    :: payload
    class(intersectionTestResult), intent(inout)  :: result
    type(elementIntersectionTestPayload), pointer :: payloadPtr
    type(elementIntersectionTestResult), pointer  :: resultPtr
    real(defReal), dimension(3)                   :: centroid, faceCentroid, outwardNormal, rEnd
    integer(shortInt)                             :: i
    real(defReal)                                 :: centroidLambda, dotProduct, faceLambda, minLambda
    character(*), parameter                       :: here = 'intersects_Ray (element_class.f90)'

    ! Downcast payload to correct type.
    select type(ptr => payload)
      type is(elementIntersectionTestPayload)
        payloadPtr => ptr

      class default
        call fatalError(here, 'Invalid payload type.')

    end select

    ! Allocate result to correct return type then associate pointer.
    select type(ptr => result)
      type is(elementIntersectionTestResult)
        resultPtr => ptr
        call resetElementIntersectionTestResult(resultPtr)

      class default
        ! Should never happen.
        call fatalError(here, 'Failed to downcast result.')

    end select

    ! Check if ray originates from inside the element (skip bounding box intersection in this case.)
    if (payloadPtr % skipBoundingBoxIntersectionTest) then
      ! Retrieve element's centroid then loop over all faces in the element.
      rEnd = payload % r + payload % u * payload % dMax
      centroid = self % getCentroid()
      minLambda = INF
      do i = 1, size(self % orientatedFaces)
        ! Retrieve the signed normal vector of the current face.
        faceCentroid = self % orientatedFaces(i) % face % ptr % getCentroid()
        outwardNormal = self % orientatedFaces(i) % outwardNormal
        
        ! Retrieve the centre of the current face and compute lambda.
        dotProduct = dot_product(rEnd - centroid, outwardNormal)
        if (areEqual(dotProduct, ZERO)) cycle
        centroidLambda = dot_product(faceCentroid - centroid, outwardNormal) / dotProduct
        
        ! If ZERO <= lambda <= ONE, append the current face to the list of potentially intersected faces.
        if (ZERO <= centroidLambda .and. centroidLambda <= ONE) then
          ! Compute lambda for the face using the actual particle coordinates.
          dotProduct = dot_product(rEnd - payload % r, outwardNormal)
          if (areEqual(dotProduct, ZERO)) cycle
          faceLambda = dot_product(faceCentroid - payload % r, outwardNormal) / dotProduct
          if (payloadPtr % excludeZeroFaces .and. faceLambda <= ZERO) cycle
          if (faceLambda < minLambda) then
            minLambda = faceLambda
            resultPtr % intersectedFace = self % orientatedFaces(i) % face

          end if

        else
          ! Check if end point is on the place of the current face.
          if (areEqual(dot_product(rEnd - faceCentroid, outwardNormal), ZERO)) then
            ! End point is on the plane of the face. Check if it is contained inside it.
            if (self % orientatedFaces(i) % face % ptr % isPointInside(rEnd)) then
              resultPtr % intersectedFace = self % orientatedFaces(i) % face
              minLambda = ONE
              exit

            end if

          end if

        end if

      end do

      if (associated(resultPtr % intersectedFace % ptr)) then
        resultPtr % intersects = .true.
        resultPtr % d = norm2(min(ONE, max(ZERO, minLambda)) * (rEnd - payload % r))

      end if

    else
      ! Call fatalError for now.
      call fatalError(here, 'Unsupported procedure.')

    end if

  end subroutine intersects_Ray

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
    self % isConvex = .false.
    if (allocated(self % childrenIdxs)) deallocate(self % childrenIdxs)
    if (allocated(self % type)) deallocate(self % type)

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
  subroutine minimumDistance(self, r, d, orientatedFace)
    class(element), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), intent(out)              :: d
    type(orientatedFaceBox), intent(out)    :: orientatedFace
    integer(shortInt)                       :: i, minIdx
    real(defReal)                           :: dFaceSquared, dSquared

    dSquared = INF
    minIdx = 0
    do i = 1, size(self % orientatedFaces)
      dFaceSquared = self % orientatedFaces(i) % face % ptr % distanceSquared(r)
      minIdx = merge(i, minIdx, dFaceSquared < dSquared)
      dSquared = min(dSquared, dFaceSquared)

    end do
    d = sqrt(dSquared)
    orientatedFace = self % orientatedFaces(minIdx)

  end subroutine minimumDistance

  !!
  !!
  !!
  pure function newElementIntersectionTestPayload(r, u, dMax, skipBoundingBoxIntersectionTest, skipZeroFaces) &
  result(payload)
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal), intent(in)               :: dMax
    logical(defBool), intent(in)            :: skipBoundingBoxIntersectionTest
    logical(defBool), intent(in), optional  :: skipZeroFaces
    type(elementIntersectionTestPayload)    :: payload

    payload % r = r
    payload % u = u
    payload % dMax = dMax
    payload % skipBoundingBoxIntersectionTest = skipBoundingBoxIntersectionTest
    if (present(skipZeroFaces)) payload % excludeZeroFaces = skipZeroFaces

  end function newElementIntersectionTestPayload

  !!
  !!
  !!
  subroutine pushFromBoundary(self, u, r)
    class(element), intent(in)                 :: self
    real(defReal), dimension(3), intent(in)    :: u
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3)                :: nudgeDirection, outwardNormal
    integer(shortInt)                          :: i

    ! Initialise nudgeDirection = ZERO then loop over all the faces in the element.
    nudgeDirection = ZERO
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
      nudgeDirection = nudgeDirection / norm2(nudgeDirection)

    else
      nudgeDirection = u

    end if
    r = r + nudgeDirection * NUDGE

  end subroutine pushFromBoundary

  !!
  !!
  !!
  subroutine sampleInitialPosition(self, rand, localId, r)
    class(element), intent(in)               :: self
    type(RNG), intent(inout)                :: rand
    integer(shortInt), intent(out)           :: localId
    real(defReal), dimension(3), intent(out) :: r
    integer(shortInt)                        :: i
    real(defReal)                            :: factorsProduct, factorsProductTimeRandomNumber3
    real(defReal), dimension(2)              :: factors
    real(defReal), dimension(3)              :: randomNumbers
    real(defReal), dimension(4)              :: barycentricWeights
    real(defReal), dimension(3, 2)           :: boundingBoxBounds
    type(inclusionTestResult)                :: inclusionResult

    localId = self % localId

    ! First check if the element is a tetrahedron and perform a direct sampling using barycentric coordinates if yes.
    if (size(self % vertices) == 4) then
      ! Sample three random numbers.
      call rand % generate(randomNumbers)

      ! Apply transformations to ensure uniform volume sampling.
      factors(1) = randomNumbers(1) ** THIRD
      factors(2) = sqrt(randomNumbers(2))

      ! Calculate barycentric weights.
      factorsProduct = product(factors)
      factorsProductTimeRandomNumber3 = factorsProduct * randomNumbers(3)
      barycentricWeights(1) = ONE - factors(1)
      barycentricWeights(2) = factors(1) - factorsProduct
      barycentricWeights(3) = factorsProduct - factorsProductTimeRandomNumber3
      barycentricWeights(4) = factorsProductTimeRandomNumber3

      ! Sample initial position.
      r = ZERO
      do i = 1, 4
        r = r + barycentricWeights(i) * self % vertices(i) % ptr % getCoordinates()

      end do

    else
      ! Retrieve bounds of element bounding box.
      boundingBoxBounds = self % getBoundingBoxBounds()
      inclusionResult % status = OUTSIDE_ELEMENT

      ! Sample initial position until the point is inside the element.
      do while (.not. inclusionResult % status == INSIDE_ELEMENT)
        ! Sample three random numbers.
        call rand % generate(randomNumbers)
        r = (boundingBoxBounds(:, 2) - boundingBoxBounds(:, 1)) * randomNumbers + boundingBoxBounds(:, 1)
        inclusionResult = self % isPointInside(r)

      end do

    end if

  end subroutine sampleInitialPosition

  !!
  !!
  !!
  subroutine resetElementIntersectionTestResult(result)
    type(elementIntersectionTestResult), intent(inout) :: result

    call resetIntersectionTestResult(result)
    result % intersectedFace % ptr => null()

  end subroutine resetElementIntersectionTestResult

  !!
  !!
  !!
  elemental subroutine setLocalId(self, localId)
    class(element), intent(inout) :: self
    integer(shortInt), intent(in) :: localId

    self % localId = localId

  end subroutine setLocalId

end module element_class