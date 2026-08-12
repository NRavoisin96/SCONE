module element_class

  use axisAlignedBoundingBox_class,  only : axisAlignedBoundingBox
  use extentTopologicalObject_inter, only : buildExtentTopologicalObjectPayload, extentTopologicalObject, &
                                            intersects_Ray_super => intersects_Ray
  use edge_class,                    only : edgeBox
  use errors_mod,                    only : fatalError
  use face_class,                    only : faceBox, orientatedFaceBox, face
  use genericProcedures,             only : append, areEqual, crossProduct, findCommon, numToChar
  use numPrecision
  use publicObjects,                 only : basicElementInfo, intersectionTestPayload, intersectionTestResult, &
                                            resetIntersectionTestResult
  use RNG_class,                     only : RNG
  use topologicalObject_inter,       only : buildTopologicalObjectPayload, kill_super => kill, topologicalObject, &
                                            topologicalObjectBox
  use universalVariables,            only : FOURTH, INSIDE_ELEMENT, INF, NUDGE, ON_BOUNDARY_ELEMENT, ONE, OUTSIDE_ELEMENT, &
                                            SIXTH, ZERO, VALENCE
  use vertex_class,                  only : vertexBox

  use ratint 

  use limb_class



  
  implicit none
  private

  ! Public procedures.
  public :: castElementPtr, newElementIntersectionTestPayload, resetElementIntersectionTestResult

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
    procedure          :: hybridIsPointInside
    procedure          :: hybridIsPointInsideGivenFaces
    procedure          :: hybridIsPointInsideTPO
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
    integer(shortInt), dimension(VALENCE) :: currentFaceIdxs
    integer(shortInt) :: front
    real(defReal), dimension(3) :: intersectionPt
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

  !!
  !!
  !!
  function castElementPtr(source, fatal) result(ptr)
    class(topologicalObject), intent(in)   :: source
    logical(defBool), intent(in), optional :: fatal
    logical(defBool)                       :: throwError
    type(element), pointer                 :: ptr
    character(*), parameter                :: HERE = 'castElementPtr (element_class.f90)'

    ! Downcast.
    select type(temp => source)
      type is(element)
        ptr => temp

      class default
        ptr => null()

    end select

    ! Throw error if requested.
    throwError = .true.
    if(present(fatal)) throwError = fatal
    if(throwError .and. .not. associated(ptr)) call fatalError(HERE, "Topological object is not of type 'element'.")

  end function castElementPtr

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

  ! subroutine intersects_Ray1(self, payload, result)
  !   class(element), intent(in)                    :: self
  !   class(intersectionTestPayload), intent(in)    :: payload
  !   class(intersectionTestResult), intent(inout)  :: result
  !   type(elementIntersectionTestPayload), pointer :: payloadPtr
  !   type(elementIntersectionTestResult), pointer  :: resultPtr
  !   real(defReal), dimension(3)                   :: centroid, faceCentroid, outwardNormal, rEnd
  !   integer(shortInt)                             :: i
  !   real(defReal)                                 :: centroidLambda, dotProduct, faceLambda, minLambda
  !   character(*), parameter                       :: here = 'intersects_Ray (element_class.f90)'

  !   ! Downcast payload to correct type.
  !   select type(ptr => payload)
  !     type is(elementIntersectionTestPayload)
  !       payloadPtr => ptr

  !     class default
  !       call fatalError(here, 'Invalid payload type.')

  !   end select

  !   ! Allocate result to correct return type then associate pointer.
  !   select type(ptr => result)
  !     type is(elementIntersectionTestResult)
  !       resultPtr => ptr
  !       call resetElementIntersectionTestResult(resultPtr)

  !     class default
  !       ! Should never happen.
  !       call fatalError(here, 'Failed to downcast result.')

  !   end select

  !   resultPtr%front = 0

  !   ! Check if ray originates from inside the element (skip bounding box intersection in this case.)
  !   if (payloadPtr % skipBoundingBoxIntersectionTest) then

    
  !     ! Retrieve element's centroid then loop over all faces in the element.
  !     rEnd = payload % r + payload % u * payload % dMax
  !     centroid = self % getCentroid()
  !     minLambda = INF
  !     do i = 1, size(self % orientatedFaces)
  !       ! Retrieve the signed normal vector of the current face.
  !       faceCentroid = self % orientatedFaces(i) % face % ptr % getCentroid()
  !       outwardNormal = self % orientatedFaces(i) % outwardNormal
        
  !       ! Retrieve the centre of the current face and compute lambda.
  !       dotProduct = dot_product(rEnd - centroid, outwardNormal)
  !       if (areEqual(dotProduct, ZERO)) cycle
  !       centroidLambda = dot_product(faceCentroid - centroid, outwardNormal) / dotProduct
        
  !       ! If ZERO <= lambda <= ONE, append the current face to the list of potentially intersected faces.
  !       if (ZERO <= centroidLambda .and. centroidLambda <= ONE) then
  !         ! Compute lambda for the face using the actual particle coordinates.
  !         dotProduct = dot_product(rEnd - payload % r, outwardNormal)
  !         if (areEqual(dotProduct, ZERO)) cycle
  !         faceLambda = dot_product(faceCentroid - payload % r, outwardNormal) / dotProduct
  !         if (payloadPtr % excludeZeroFaces .and. faceLambda <= ZERO) cycle
  !         if (faceLambda < minLambda) then
  !           minLambda = faceLambda
  !           resultPtr % intersectedFace = self % orientatedFaces(i) % face

  !         end if

  !       else
  !         ! Check if end point is on the place of the current face.
  !         if (areEqual(dot_product(rEnd - faceCentroid, outwardNormal), ZERO)) then
  !           ! End point is on the plane of the face. Check if it is contained inside it.
  !           if (self % orientatedFaces(i) % face % ptr % isPointInside(rEnd)) then
  !             resultPtr % intersectedFace = self % orientatedFaces(i) % face
  !             minLambda = ONE
  !             exit

  !           end if

  !         end if

  !       end if

  !     end do

  !     if (associated(resultPtr % intersectedFace % ptr)) then
  !       resultPtr % intersects = .true.
  !       resultPtr % d = norm2(min(ONE, max(ZERO, minLambda)) * (rEnd - payload % r))

  !     end if

  !   else
  !     ! Call fatalError for now.
  !     call fatalError(here, 'Unsupported procedure.')

  !   end if

  ! end subroutine intersects_Ray1


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
    integer(shortInt)                             :: i, zeroFaceCount, faceArrayFront
    real(defReal)                                 :: centroidLambda, dotProduct, faceLambda, minLambda, newMinLambda
    type(faceBox) :: tempFace
    type(orientatedFaceBox), dimension(size(self%orientatedFaces)) :: faceArray
    real(defReal), dimension(size(self%orientatedFaces)) :: faceArrayLambdas
    logical :: intersected
    character(*), parameter                       :: here = 'intersects_Ray (element_class.f90)'

    faceArrayFront = 1

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


      resultPtr%front = 0

      
      centroid = self % getCentroid()
      minLambda = INF

      do i = 1, size(self % orientatedFaces)
        if (payload%front > 0 .and. &
            ANY(payload%currentFaceIdxs(1:payload%front)==self%orientatedFaces(i)%face%ptr%getIdx()) .and. &
              self%orientatedFaces(i)%face%ptr%getIdx() /= 0) then 
          cycle 
        end if
  
        
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


          if (faceLambda >= 0 .and. faceLambda <= 1) then 
            if (faceArrayFront == 1) then 

              faceArray(1) = self % orientatedFaces(i)
              faceArrayLambdas(1) = faceLambda
              faceArrayFront = faceArrayFront + 1 
            else 


              faceArray(faceArrayFront) = self % orientatedFaces(i) 
              faceArrayLambdas(faceArrayFront) = faceLambda
              faceArrayFront = faceArrayFront + 1 
            end if
          end if


          if (faceLambda < minLambda) then
            minLambda = faceLambda
            tempFace = self % orientatedFaces(i) % face

          end if

        else
          ! Check if end point is on the place of the current face.
          if (areEqual(dot_product(rEnd - faceCentroid, outwardNormal), ZERO)) then

            ! End point is on the plane of the face. Check if it is contained inside it.
            if (self % orientatedFaces(i) % face % ptr % isPointInside(rEnd)) then
              
              tempFace = self % orientatedFaces(i) % face
              minLambda = ONE
              exit

            end if

          end if

        end if

      end do



      if (associated(tempFace%ptr)) then

        !!! code added here: check for distance less than epsilon to a vertex or edge
        if (existsEpsilonDistance(self, payload, payloadPtr, tempFace, &
            size(tempFace%ptr%getVertices()), size(tempFace%ptr%getEdges()), &
              (payload % r + (payload%u*payload%dMax)))) then 

          call rescueParticleNearVertexEdge(self, payload, payloadPtr, faceArrayFront-1,  &
                                            faceArray, resultPtr, intersected, newMinLambda)

          if (associated(resultPtr % intersectedFace % ptr) ) then

            resultPtr % intersects = .true.
            resultPtr % intersectionPt = payload % r + min(ONE, max(ZERO, newMinLambda)) * (rEnd - payload % r)
            resultPtr % d = norm2(min(ONE, max(ZERO, newMinLambda)) * (rEnd - payload % r))

          end if

          return

        else 
          do i=1, faceArrayFront-1 

            if (areEqual(minLambda, faceArrayLambdas(i))) then! .and. faceArray(i)%ptr%getFaceIdx() /= tempFace%ptr%getFaceIdx()) then 
              resultPtr%front = resultPtr%front+1
              resultPtr % currentFaceIdxs(resultPtr%front) = faceArray(i)%face%ptr%getIdx() 
            end if 
          end do

          
          resultPtr % intersectedFace = tempFace
          resultPtr % intersects = .true.
          resultPtr % intersectionPt = payload%r + min(ONE, max(ZERO, minLambda)) * (rEnd - payload % r)
          resultPtr % d = norm2(min(ONE, max(ZERO, minLambda)) * (rEnd - payload % r))
        end if

      end if



    else
      ! Call fatalError for now.
      call fatalError(here, 'Unsupported procedure.')

    end if


  end subroutine intersects_Ray


  ! subroutine ratintNormal(face, numVertices, rationalNormal)
  !   type(orientatedFaceBox), intent(in) :: face 
  !   type(ratint_t), dimension(3), intent(inout) :: rationalNormal
  !   integer, intent(in) :: numVertices
  !   type(ratint_t), dimension(3) :: v1, v2, v3, dir1, dir2, centroidDir
  !   type(vertexBox), dimension(numVertices) :: vertices
  !   integer :: i


  !   vertices = face%face%ptr%getVertices()
  !   v1 = vertices(1)%ptr%getRatintCoordinates()
  !   v2 = vertices(2)%ptr%getRatintCoordinates()
  !   v3 = vertices(3)%ptr%getRatintCoordinates()

  !   dir1 = v1 - v2 
  !   dir2 = v1 - v3

  !   print * , 'NORMAL'
  !   rationalNormal = crossProduct_ratint(dir1, dir2)

  !   ! print *, '!!!!!!'
  !   ! do i=1, 3 
  !   !   call printRatInt(rationalNormal(i))
  !   ! end do 

  ! end subroutine ratintNormal


  ! function crossProduct_ratint(a, b) result(c)
  !   type(ratint_t), dimension(3), intent(in) :: a, b
  !   type(ratint_t), dimension(3) :: c
  !   print *, 'CROSSPRODUCT'

  !   print *, 'a2'
  !   call printRatInt(a(1))
  !   print *, 'b3'
  !   call printRatInt(b(2))
  !   print *, 'a2b3'
  !   call printRatInt(a(1) * b(2))

  !   print *, 'a3'
  !   call printRatInt(a(2))
  !   print *, 'b2'
  !   call printRatInt(b(1))
  !   print *, 'a3b2'
  !   call printRatInt(a(2) * b(1))
  !   print *, 'a2b3-a3b2'
  !   call printRatInt((a(1) * b(2)) - (a(2) * b(1)))

  !   c = [a(2)*b(3) - a(3)*b(2), &
  !        a(3)*b(1) - a(1)*b(3), &
  !        a(1)*b(2) - a(2)*b(1)]

  ! end function crossProduct_ratint


!!!NOTE, THROUGH ANOTHER ISSUE, SOMETHING IS WRONG HERE FOR BOUNDARY EXITS?? MAYBE DIRECTION??
  subroutine rescueParticleNearVertexEdge(elementInp, payload, payloadPtr, faceArrayLast, faceArray, res, intersects,newMinLambda)
    class(element), target, intent(in) :: elementInp 
    class(intersectionTestPayload), intent(in) :: payload 
    type(elementIntersectionTestPayload), pointer, intent(in) :: payloadPtr
    integer, intent(in) :: faceArrayLast
    type(elementIntersectionTestResult), intent(inout) :: res
    logical, intent(inout) :: intersects
    real(defReal), intent(inout) :: newMinLambda
    type(orientatedFaceBox), dimension(:), intent(in) :: faceArray
    type(orientatedFaceBox) :: tempFace
    type(orientatedFaceBox), dimension(size(faceArray)) :: faceArrayNew
    type(ratint_t), dimension(size(faceArray)) :: faceArrayNewLambdas
    integer(shortInt) :: faceArrayNewFront
    type(ratint_t), dimension(3) ::  rEnd, rStart, outwardNormal, faceCentroid, elemCentroid
    type(ratint_t) :: minLambda, centroidLambda, faceLambda, dotProduct, tempdot
    real(defReal), dimension(3) ::  rEndReal
    type(topologicalObjectBox), dimension(:), allocatable :: faceElements
    type(vertexBox), dimension(:), allocatable :: faceVertices
    class(element), pointer :: neighbourElem 
    integer :: i, j, k
    type(ratint_t), dimension(3) :: dummy

    intersects = .false.
    
    faceArrayNewFront = 1
    print *, 'rescue'



    call ratintElementCentroid(elementInp, size(elementInp%getVertices()), elemCentroid)

    rEndReal = payload % r + payload % u * payload % dMax

    !call ratintCalcEndpoint(payload,rEnd)

    rEnd = convert_ieee(rEndReal)

    rStart = convert_ieee(payload % r)

    minLambda = def_ratint_large()


    ! Loops through face array
    do i = 1, faceArrayLast


      !call ratintFaceCentroid(faceArray(i)%face, size(faceArray(i)%face % ptr %getVertices()), faceCentroid)


      faceVertices = faceArray(i)%face % ptr %getVertices()


      faceCentroid = faceVertices(1)%ptr%getRatintCoordinates()


      ! call ratintOutwardNormal(faceArray(i)%face, elemCentroid, faceCentroid, &
      !         size(faceArray(i)%face % ptr %getVertices()), outwardNormal)

      outwardNormal = faceArray(i)%ratintOutwardNormal

      

      !call ratintNormal(faceArray(i), size(faceVertices), dummy)

      ! Retrieve the centre of the current face and compute lambda.
      dotProduct = dot_product(rEnd - elemCentroid, outwardNormal)
      

      if (isZero(dotProduct)) cycle
      
      centroidLambda = dot_product(faceCentroid - elemCentroid, outwardNormal) / dotProduct


      ! If ZERO <= lambda <= ONE, append the current face to the list of potentially intersected faces.
      if (centroidLambda >= convert_int(0_8) .and. convert_int(1_8) >= centroidLambda) then
        ! Compute lambda for the face using the actual particle coordinates.
        dotProduct = dot_product(rEnd - rStart, outwardNormal)
        if (isZero(dotProduct)) cycle
        faceLambda = dot_product(faceCentroid - rStart, outwardNormal) / dotProduct


        if (payloadPtr % excludeZeroFaces .and.  convert_int(0_8) >= faceLambda) cycle

      
        if (faceLambda >= convert_int(0_8) .and. convert_int(1_8) >= faceLambda) then 
          !print *, evaluate(faceLambda)
          

          faceArrayNew(faceArrayNewFront)  = faceArray(i) 
          faceArrayNewLambdas(faceArrayNewFront) = faceLambda 
          faceArrayNewFront = faceArrayNewFront + 1

        end if


        if (minLambda > faceLambda) then

          minLambda = faceLambda
          tempFace = faceArray(i)
        end if

      else
        ! Check if end point is on the place of the current face.
        if (isZero(dot_product(rEnd - faceCentroid, outwardNormal))) then
          ! End point is on the plane of the face. Check if it is contained inside it.
          if (hybridIsPointInsideFace(faceArray(i)%face, rEnd, rEndReal, outwardNormal, &
                size(faceArray(i)%face%ptr%getVertices()))) then
            tempFace = faceArray(i)
            minLambda = convert_int(1_8)
            exit
          end if

        end if

      end if

    end do


    if (associated(tempFace%face%ptr)) then
      outwardNormal = tempFace%ratintOutwardNormal
      res%intersectedFace = tempFace%face
      !call ratintFaceCentroid(res%intersectedFace, size(res%intersectedFace % ptr %getVertices()), faceCentroid)

      !call ratintOutwardNormal(res%intersectedFace, elemCentroid, faceCentroid, &
      !      size(res%intersectedFace % ptr %getVertices()), outwardNormal)
      
      intersects = checkIntersected(res%intersectedFace, rStart, rEnd, size(res%intersectedFace%ptr%getVertices()), outwardNormal)
      
      print *, 'minlambda'
      print *, evaluate(minLambda)
      print *, 'intersection lambdas'
      do i =1, faceArrayNewFront-1 
        print *, evaluate(faceArrayNewLambdas(i))
        if (minLambda==faceArrayNewLambdas(i)) then! .and. faceArrayNew(i)%ptr%getFaceIdx()/=res%intersectedFace%ptr%getFaceIdx()) then 
          !print *, evaluate(minLambda)
          res%front = res%front + 1
          res%currentFaceIdxs(res%front) = faceArrayNew(i)%face%ptr%getIdx()
        end if 
      end do

    end if 


    newMinLambda = evaluate(minLambda)
 
  end subroutine rescueParticleNearVertexEdge



  subroutine rescueParticleNearVertexEdgeWrong(elementInp,payload,payloadPtr,faceArrayLast,faceArray,res,intersects,newMinLambda)
    class(element), target, intent(in) :: elementInp 
    class(intersectionTestPayload), intent(in) :: payload 
    type(elementIntersectionTestPayload), pointer, intent(in) :: payloadPtr
    integer, intent(in) :: faceArrayLast
    type(faceBox), intent(inout) :: res
    logical, intent(inout) :: intersects
    real(defReal), intent(inout) :: newMinLambda
    type(faceBox), dimension(:) :: faceArray
    type(ratint_t), dimension(3) ::  rEnd, rStart, outwardNormal, faceCentroid, elemCentroid
    type(ratint_t) :: minLambda, centroidLambda, faceLambda, dotProduct
    real(defReal), dimension(3) ::  rEndReal
    type(topologicalObjectBox), dimension(:), allocatable :: faceElements
    class(element), pointer :: neighbourElem 
    integer(shortInt) :: minElemIdx, minElemIdxLocal

    

    integer :: i, j, k

    call ratintElementCentroid(elementInp, size(elementInp%getVertices()), elemCentroid)
    !! calculate directly
    rEndReal = payload % r + payload % u * payload % dMax
    ! print *, '????'
    ! print *, rEndReal
    !call ratintCalcEndpoint(payload, rEnd)
    rEnd = convert_ieee(rEndReal)
    rStart = convert_ieee(payload % r)
    minLambda = def_ratint_large()

    ! print *, evaluate(rEnd(1))
    ! print *, '?'
    ! call printRatInt(rEnd(1))


    ! Loops through face array
    do i = 1, faceArrayLast
      faceElements = faceArray(i) % ptr % getSharingElements()
      ! Finds the element that isn't the current element
      if (size(faceElements) == 1) then 
        !print *, '?'
        !call printRatInt(minLambda)
        call rescueParticleNearVertexEdgeBoundary(elementInp, payload, payloadPtr, faceArray(i), res, minLambda)
      end if
      do j = 1, size(faceElements)
        ! Downcast element to correct type.
        select type(ptr => faceElements(j) % ptr)
          type is (element)
            if (.not. associated(ptr, elementInp)) then

              ! Checks the lambda values for the element faces that aren't in the intersected faces
              do k = 1, size(ptr%getOrientatedFaces())

                if (associated(faceArray(i) % ptr, ptr%orientatedFaces(k)%face%ptr)) then
                  cycle
                end if


                call ratintFaceCentroid(ptr % orientatedFaces(k) % face, & 
                  size(ptr % orientatedFaces(k) % face % ptr %getVertices()), faceCentroid)


                call ratintOutwardNormal(ptr % orientatedFaces(k) % face, elemCentroid, faceCentroid, &
                  size(ptr % orientatedFaces(k) % face % ptr %getVertices()), outwardNormal)
                
                ! Retrieve the centre of the current face and compute lambda.
                dotProduct = dot_product(rEnd - elemCentroid, outwardNormal)
                if (isZero(dotProduct)) cycle

                centroidLambda = dot_product(faceCentroid - elemCentroid, outwardNormal) / dotProduct

                ! If ZERO <= lambda <= ONE, append the current face to the list of potentially intersected faces.
                if (centroidLambda >= convert_int(0_8) .and. convert_int(1_8) >= centroidLambda) then
                  ! Compute lambda for the face using the actual particle coordinates.
                  dotProduct = dot_product(rEnd - rStart, outwardNormal)
                  if (isZero(dotProduct)) cycle
                  faceLambda = dot_product(faceCentroid - rStart, outwardNormal) / dotProduct
                  if (payloadPtr % excludeZeroFaces .and.  convert_int(0_8) >= faceLambda) cycle

                  if (i == 1) then 
                      minLambda = faceLambda
     
                  else if (minLambda > faceLambda) then
                      minLambda = faceLambda
                      res = ptr % orientatedFaces(k) % face
                      minElemIdx = ptr%getIdx() 
                      minElemIdxLocal = ptr%getLocalId()
                  end if

                else
                  ! Check if end point is on the place of the current face.
                  if (isZero(dot_product(rEnd - faceCentroid, outwardNormal))) then
                    ! End point is on the plane of the face. Check if it is contained inside it.
                    if (hybridIsPointInsideFace(ptr % orientatedFaces(k) % face, rEnd, rEndReal, outwardNormal, &
                          size(ptr % orientatedFaces(k) % face%ptr%getVertices()))) then
                      res = ptr % orientatedFaces(k) % face
                      minLambda = convert_int(1_8)
                      exit
                    end if

                  end if

                end if
              end do

              exit

            end if

          class default

        end select

      end do

      

    end do

    !call printRatInt(minLambda)
    !print *, evaluate(minLambda)

    !call printRatInt(minLambda)
    call ratintOutwardNormal(res, elemCentroid, faceCentroid, size(res % ptr %getVertices()), outwardNormal)
    intersects = checkIntersected(res, rStart, rEnd, size(res%ptr%getVertices()), outwardNormal)
    newMinLambda = evaluate(minLambda)


  end subroutine rescueParticleNearVertexEdgeWrong



  subroutine rescueParticleNearVertexEdgeBoundary(elementInp, payload, payloadPtr, face, res, minLambda)
    class(element), target, intent(in) :: elementInp 
    class(intersectionTestPayload), intent(in) :: payload 
    type(elementIntersectionTestPayload), pointer, intent(in) :: payloadPtr
    type(faceBox), intent(in) :: face
    type(faceBox), intent(inout) :: res
    type(ratint_t), intent(inout) :: minLambda
    type(ratint_t), dimension(3) ::  rEnd, rStart, outwardNormal, faceCentroid, elemCentroid
    type(ratint_t) :: centroidLambda, faceLambda, dotProduct
    real(defReal), dimension(3) ::  rEndReal
    type(topologicalObjectBox), dimension(:), allocatable :: faceElements
    class(element), pointer :: neighbourElem 
    integer :: i, j, k


    call ratintElementCentroid(elementInp, size(elementInp%getVertices()), elemCentroid)
    !! calculate directly
    rEndReal = payload % r + payload % u * payload % dMax
    call ratintCalcEndpoint(payload, rEnd)
    rStart = convert_ieee(payload % r)
    minLambda = def_ratint_large()


    call ratintFaceCentroid(face, size(face % ptr %getVertices()), faceCentroid)


    call ratintOutwardNormal(face, elemCentroid, faceCentroid, size(face % ptr %getVertices()), outwardNormal)
    
    ! Retrieve the centre of the current face and compute lambda.
    dotProduct = dot_product(rEnd - elemCentroid, outwardNormal)
    if (isZero(dotProduct)) then 
      return 
    end if 

    centroidLambda = dot_product(faceCentroid - elemCentroid, outwardNormal) / dotProduct

    ! If ZERO <= lambda <= ONE, append the current face to the list of potentially intersected faces.
    if (centroidLambda >= convert_int(0_8) .and. convert_int(1_8) >= centroidLambda) then
      ! call printRatInt(centroidLambda)
      !print *,' here'
      ! Compute lambda for the face using the actual particle coordinates.
      dotProduct = dot_product(rEnd - rStart, outwardNormal)
      if (isZero(dotProduct)) then 
        !print *,'*'
        return 
      end if
      faceLambda = dot_product(faceCentroid - rStart, outwardNormal) / dotProduct
      if (payloadPtr % excludeZeroFaces .and.  convert_int(0_8) >= faceLambda) then 
        !print *,'**'
        return 
      end if
      !call printRatInt(faceLambda)

      if (minLambda > faceLambda) then
          !print *, 'ummmm'
          minLambda = faceLambda
          res =face
      end if

    else
      ! Check if end point is on the plane of the current face.
      if (isZero(dot_product(rEnd - faceCentroid, outwardNormal))) then
        ! End point is on the plane of the face. Check if it is contained inside it.
        if (hybridIsPointInsideFace(face, rEnd, rEndReal, outwardNormal, size(face%ptr%getVertices()))) then
          res = face
          minLambda = convert_int(1_8)
        end if

      end if

    end if

  
    !call printRatInt(minLambda)

  end subroutine rescueParticleNearVertexEdgeBoundary




  function checkIntersected(face, rStart, rEnd, numVertices, normal) result(intersects)
    type(faceBox), intent(in) :: face 
    type(ratint_t), dimension(3) :: rStart, rEnd, planePt, normal
    logical :: intersects 
    integer :: numVertices, i
    real(defReal), dimension(3) :: vertex 
    type(vertexBox), dimension(numVertices) :: vertices
    type(ratint_t) :: dot1, dot2
    

    intersects = .true.

    vertices = face%ptr%getVertices() 

    vertex = vertices(1)%ptr%getCoordinates()


    planePt = convert_ieee(vertex)

    dot1 = dot_product(rStart - planePt, normal)

    dot2 = dot_product(rEnd - planePt, normal)

    ! call printRatInt(dot1)
    ! call printRatInt(dot2)


    ! Allow equal to?
    if ((dot1 > convert_int(0_8) .and. dot2 > convert_int(0_8)) .or. (convert_int(0_8) > dot1 .and. convert_int(0_8) > dot2)) then
      intersects = .false.
      return 
    end if


  end function checkIntersected


  function existsEpsilonDistance(elementInp, payload, payloadPtr, minFace, numVertices, numEdges, & 
                                  intersectionPt) result(check)
    class(element), intent(in) :: elementInp 
    class(intersectionTestPayload), intent(in) :: payload 
    type(elementIntersectionTestPayload), pointer, intent(in) :: payloadPtr
    type(faceBox), intent(in) :: minFace 
    integer, intent(in) :: numVertices
    integer, intent(in) :: numEdges
    real(real64), dimension(3), intent(in) :: intersectionPt
    real(real64) :: minLambda
    type(vertexBox), dimension(numVertices) :: vertices
    type(edgeBox), dimension(numEdges) :: edges
    real(real64), dimension(3) :: coords
    real(real64), dimension(3) :: dist
    type(vertexBox), dimension(2) :: edgeVertices
    type(vertexBox) :: vertex
    real(real64), dimension(3) :: edgeVertexCoordsA, edgeVertexCoordsB, edgeDirection

    integer :: i 
    logical :: check 

    check = .false.

    vertices = minFace%ptr%getVertices()
    edges = minFace%ptr%getEdges()

    do i=1, numVertices 
      ! print *, 'vertex'
      ! print *,vertices(i)%ptr%getCoordinates()
      dist = norm2(vertices(i)%ptr%getCoordinates() - intersectionPt)
      !print *, dist
      if (areEqual(dist, ZERO)) then 
        check = .true. 
        return 
      end if
    end do

    do i=1, numEdges 
      edgeVertices = (edges(i)%ptr%getVertices())
      edgeVertexCoordsA = edgeVertices(1)%ptr%getCoordinates()
      edgeVertexCoordsB = edgeVertices(2)%ptr%getCoordinates()
      edgeDirection = edgeVertexCoordsB - edgeVertexCoordsA
      
      dist = crossProduct((intersectionPt - edgeVertexCoordsA), edgeDirection)
      dist = dist / norm2(edgeDirection)
      !print *, dist
      if (areEqual(dist, ZERO)) then 
        check = .true. 
        return 
      end if
    end do

  end function existsEpsilonDistance


  subroutine getFaceIDs(faces, numFaces, faceIDs)
    integer(shortInt), intent(in) :: numFaces
    type(topologicalObjectBox), dimension(numFaces), intent(in) :: faces
    integer(shortInt), dimension(numFaces) :: faceIDs
    integer(shortInt) :: i 


    do i=1, numFaces 
      select type(ptr1 => faces(i)%ptr)
        type is(face)
          faceIDs(i) = ptr1%getIdx()
      end select

    end do 



  end subroutine getFaceIDs



  function hybridIsPointInsideTPO(self, r, u, faces, numFaces) result(result)
    class(element), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    integer(shortInt), intent(in) :: numFaces
    type(topologicalObjectBox), dimension(numFaces), intent(in) :: faces
    
    integer(shortInt), dimension(numFaces) :: faceIDs
    type(inclusionTestResult)               :: result
    integer(shortInt)                       :: i
    real(defReal)                           :: dotProduct
    logical(defBool)                        :: isOnBoundary

    isOnBoundary = .false.
    result % status = INSIDE_ELEMENT

    call getFaceIDs(faces, size(faces), faceIDs)

    do i=1, size(self%orientatedFaces)
      if (ANY(faceIDs==self%orientatedFaces(i)%face%ptr%getIdx())) then 

        dotProduct = dot_product(self%orientatedFaces(i)%outwardNormal, u)

        if (areEqual(dotProduct, ZERO)) then 
          if (dotProduct /= 0) then 
            
            result = ratintIsPointInsideGivenFaces(self, r, u, faceIDs)
          end if 

          isOnBoundary = .true.

          ! Store the first face found and cycle to search other faces.
          if (result % failedFaceIdx == 0) then 
            result % failedFaceIdx = self % orientatedFaces(i) % face % ptr % getIdx()
            return 
          end if
      
        end if

        if (dotProduct > 0) then 
          result % status = OUTSIDE_ELEMENT
          result % failedFaceIdx = self % orientatedFaces(i) % face % ptr % getIdx()
          return 
        end if
      
      else 
        cycle
      end if

    end do

    if (isOnBoundary) result % status = ON_BOUNDARY_ELEMENT


  end function hybridIsPointInsideTPO


  function hybridIsPointInsideGivenFaces(self, r, u, faces) result(result)
    class(element), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    integer(shortInt), dimension(:), intent(in) :: faces
    type(inclusionTestResult)               :: result
    integer(shortInt)                       :: i
    real(defReal)                           :: dotProduct
    logical(defBool)                        :: isOnBoundary

    isOnBoundary = .false.
    result % status = INSIDE_ELEMENT

    do i=1, size(self%orientatedFaces)
      if (ANY(faces==self%orientatedFaces(i)%face%ptr%getIdx())) then 

        dotProduct = dot_product(self%orientatedFaces(i)%face%ptr%getNormal(), u)

        if (areEqual(dotProduct, ZERO)) then 
          if (dotProduct /= 0) then 
            
            result = ratintIsPointInsideGivenFaces(self, r, u, faces)
          end if 

          isOnBoundary = .true.

          ! Store the first face found and cycle to search other faces.
          if (result % failedFaceIdx == 0) then 
            result % failedFaceIdx = self % orientatedFaces(i) % face % ptr % getIdx()
            return 
          end if
      
        end if

        if (dotProduct > 0) then 
          result % status = OUTSIDE_ELEMENT
          result % failedFaceIdx = self % orientatedFaces(i) % face % ptr % getIdx()
          return 
        end if
      
      else 
        cycle
      end if

    end do

    if (isOnBoundary) result % status = ON_BOUNDARY_ELEMENT


  end function hybridIsPointInsideGivenFaces


  function ratintIsPointInsideGivenFaces(self, r, u, faces) result(result)
    class(element), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    integer(shortInt), dimension(:), intent(in) :: faces
    type(inclusionTestResult)               :: result
    integer(shortInt)                       :: i
    type(ratint_t)                           :: dotProduct
    logical(defBool)                        :: isOnBoundary
    type(ratint_t), dimension(3) :: ratintNormal, ratintU, rintFaceCentroid, rintElementCentroid


    isOnBoundary = .false.
    result % status = INSIDE_ELEMENT
    ratintU = convert_ieee(u)

    call ratintElementCentroid(self, size(self%getVertices()), rintElementCentroid)

    do i=1, size(self%orientatedFaces)
      if (ANY(faces==self%orientatedFaces(i)%face%ptr%getIdx())) then 

        call ratintFaceCentroid(self % orientatedFaces(i) % face, size(self % orientatedFaces(i) % face%ptr%getVertices()), &
                            rintFaceCentroid)


        call ratintOutwardNormal(self % orientatedFaces(i) % face, rintElementCentroid, rintFaceCentroid, &
              size(self % orientatedFaces(i) % face%ptr%getVertices()), ratintNormal)

        !call ratintElementCentroid(self, size(self%getVertices()), rintElementCentroid)

        !call ratintFaceCentroid(self%orientatedFaces(i)%face, size(self%orientatedFaces(i)%face%ptr%getVertices()), rintFaceCentroid)

        !ratintNormal = self%orientatedFaces(i)%ratintOutwardNormal

        dotProduct = dot_product(ratintNormal, ratintU)
        if (isZero(dotProduct)) then 

          isOnBoundary = .true.

          ! Store the first face found and cycle to search other faces.
          if (result % failedFaceIdx == 0) then 
            result % failedFaceIdx = self % orientatedFaces(i) % face % ptr % getIdx()
            return 
          end if
      
        end if

        if (dotProduct > convert_int(0_8)) then 
          result % status = OUTSIDE_ELEMENT
          result % failedFaceIdx = self % orientatedFaces(i) % face % ptr % getIdx()
          return 
        end if
      
      else 
        cycle
      end if

    end do

    if (isOnBoundary) result % status = ON_BOUNDARY_ELEMENT



  end function ratintIsPointInsideGivenFaces



  function hybridIsPointInside(self, r) result(result)
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
        if (dotProduct /= ZERO) then 
          result = ratintIsPointInsideElement(self, r)
          return
        end if
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


  end function hybridIsPointInside


  function ratintIsPointInsideElement(self, r) result(result)
    class(element), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    type(inclusionTestResult)               :: result
    integer(shortInt)                       :: i
    type(ratint_t)                          :: dotProduct
    logical(defBool)                        :: isOnBoundary
    type(ratint_t), dimension(3) :: ratintR, rintFaceCentroid, ratintNormal, rintElementCentroid
    type(vertexBox), dimension(:), allocatable :: vertices

    ratintR = convert_ieee(r)

    

    call ratintElementCentroid(self, size(self%getVertices()), rintElementCentroid)

  
    isOnBoundary = .false.
    result % status = INSIDE_ELEMENT
    do i = 1, size(self % orientatedFaces)
      call ratintFaceCentroid(self%orientatedFaces(i)%face, size(self%orientatedFaces(i)%face%ptr%getVertices()), rintFaceCentroid)


      vertices = self%orientatedFaces(i)%face%ptr%getVertices()

      !rintFaceCentroid = vertices(1)%ptr%getRatintCoordinates()


      call ratintOutwardNormal(self % orientatedFaces(i) % face, rintElementCentroid, &
        rintFaceCentroid, size(vertices), ratintNormal)

      ratintNormal = self%orientatedFaces(i)%ratintOutwardNormal

      ! Make a vector going from the coordinates to the face's centroid and perform the dot
      ! product between this vector and the face's normal vector.
      dotProduct = dot_product(rintFaceCentroid - ratintR, ratintNormal)

      ! Check if the point is effectively on the plane of this face.
      if (isZero(dotProduct)) then
        isOnBoundary = .true.

        ! Store the first face found and cycle to search other faces.
        if (result % failedFaceIdx == 0) result % failedFaceIdx = self % orientatedFaces(i) % face % ptr % getIdx()
        cycle

      end if

      ! If dotProduct < ZERO, update result and return early.
      if (convert_int(0_8) > dotProduct) then
        result % status = OUTSIDE_ELEMENT
        result % failedFaceIdx = self % orientatedFaces(i) % face % ptr % getIdx()
        return

      end if

    end do

    ! If point is on boundary, update result % status.
    if (isOnBoundary) result % status = ON_BOUNDARY_ELEMENT




  end function ratintIsPointInsideElement


  function hybridIsPointInsideFace(face, endPoint, fpEndPoint, normal, numVertices) result(inside)
     type(faceBox), intent(in) :: face 
      type(ratint_t), dimension(3), intent(in) :: endPoint, normal
      real(defReal),dimension(3), intent(in) :: fpEndPoint
      real(defReal), dimension(3) :: fpNorm
      real(defReal) :: dotProduct
      integer :: numVertices
      logical(defBool) :: inside
      integer(shortInt)                       :: i, nextIdx
      type(vertexBox), dimension(3)             :: vertices
      real(defReal), dimension(numVertices)             :: vertexCoords, nextCoords
      type(ratint_t), dimension(3) :: ratintCoords

      inside = .false.
      vertices = face % ptr % getVertices()
      fpNorm = face%ptr%getNormal()

      ! Loop through all the edges in the face and check if the point lies on the same side
      ! of each edge (note: this assumes a consistent vertex numbering).
      do i = 1, numVertices
        nextIdx = merge(1, i + 1, i == numVertices)
        vertexCoords = vertices(i) % ptr % getCoordinates()
        nextCoords = vertices(nextIdx)%ptr%getCoordinates()
        dotProduct = dot_product(fpNorm, crossProduct(nextCoords - vertexCoords, fpEndPoint - vertexCoords))

        if (areEqual(dotProduct, ZERO) .and. (dotProduct /= ZERO)) then 
          inside = ratintIsPointInsideFace(face, endPoint, normal, numVertices)
          return
        else if (dotProduct < 0) then 
          return
        end if

      end do

      ! If reached here, the point is inside the face.
      inside = .true.



  end function hybridIsPointInsideFace







  function ratintIsPointInsideFace(face, endPoint, normal, numVertices) result(inside)
    type(faceBox), intent(in) :: face 
    type(ratint_t), dimension(3), intent(in) :: endPoint, normal
    integer :: numVertices
    logical(defBool) :: inside
    integer(shortInt)                       :: i, nextIdx
    type(vertexBox), dimension(3)             :: vertices
    real(defReal), dimension(numVertices)             :: vertexCoords
    type(ratint_t), dimension(3) :: ratintCoords, nextRatintCoords
    type(ratint_t) :: dotProduct

    inside = .false.
    vertices = face % ptr % getVertices()

    ! Loop through all the edges in the face and check if the point lies on the same side
    ! of each edge (note: this assumes a consistent vertex numbering).
    do i = 1, numVertices
      nextIdx = merge(1, i + 1, i == numVertices)
      vertexCoords = vertices(i) % ptr % getCoordinates()
      ratintCoords = vertices(i) % ptr % getRatintCoordinates()
      nextRatintCoords = vertices(nextIdx)%ptr%getRatintCoordinates()
      dotProduct = dot_product(normal, crossProduct(nextRatintCoords - ratintCoords, endPoint - ratintCoords))

      if (convert_int(0_8) > dotProduct) then 
        return 
      end if

    end do

    ! If reached here, the point is inside the face.
    inside = .true.



  end function ratintIsPointInsideFace



  subroutine ratintCalcEndpoint(payload, endPoint)
    class(intersectionTestPayload), intent(in) :: payload 
    type(ratint_t), dimension(3), intent(inout) :: endPoint
    type(ratint_t), dimension(3) :: ratintR 
    type(ratint_t), dimension(3) :: ratintU , huhh
    type(ratint_t) :: ratintDMax, huh
    type(limb_t) ::prod1,prod2
    real(defReal), dimension(3) :: ok
    integer :: i

    
    ratintR = convert_ieee(payload % r)
    ratintU = convert_ieee(payload % u)
    ratintDMax = convert_ieee(payload % dmax)


    endPoint = ratintR + (ratintU * ratintDMax)


  end subroutine ratintCalcEndpoint



  subroutine ratintFaceCentroid(face, numVertices, rationalCentroid)
    type(faceBox), intent(in) :: face 
    type(ratint_t), dimension(3), intent(inout) :: rationalCentroid
    integer, intent(in) :: numVertices
    type(vertexBox), dimension(numVertices) :: vertices
    type(ratint_t), dimension(3) :: coords
    integer :: i 

    rationalCentroid = initratint_vector()
    vertices = (face%ptr%getVertices())
    do i = 1, size(face%ptr%getVertices())
      coords = vertices(i)%ptr%getRatintCoordinates()
      rationalCentroid = rationalCentroid + coords
    end do 

    rationalCentroid(1) = rationalCentroid(1) / convert_int(size(face%ptr%getVertices())*1_8)
    rationalCentroid(2) = rationalCentroid(2) / convert_int(size(face%ptr%getVertices())*1_8)
    rationalCentroid(3) = rationalCentroid(3) / convert_int(size(face%ptr%getVertices())*1_8)

  end subroutine ratintFaceCentroid 



  subroutine ratintOutwardNormal(face, elemCentroid, faceCentroid, numVertices, rationalNormal)
    type(faceBox), intent(in) :: face 
    type(ratint_t), dimension(3) :: elemCentroid, faceCentroid
    type(ratint_t), dimension(3), intent(inout) :: rationalNormal
    integer, intent(in) :: numVertices
    type(ratint_t), dimension(3) :: v1, v2, v3, dir1, dir2, centroidDir
    type(ratint_t) :: signTest
    type(vertexBox), dimension(numVertices) :: vertices
    real(defReal), dimension(3) :: coords
    integer :: i

    vertices = face%ptr%getVertices()
    coords = vertices(1)%ptr%getCoordinates()
    v1 = vertices(1)%ptr%getRatintCoordinates()

    coords = vertices(2)%ptr%getCoordinates()
    v2 = vertices(2)%ptr%getRatintCoordinates()

    coords = vertices(3)%ptr%getCoordinates()
    v3 = vertices(3)%ptr%getRatintCoordinates()


    dir1 = v1 - v2 
    dir2 = v1 - v3

    rationalNormal = crossProduct(dir1, dir2)

    centroidDir = faceCentroid - elemCentroid


    signTest = dot_product(rationalNormal, centroidDir)

    if (convert_int(0_8) > signTest) then 
      call swapSign(rationalNormal)
    end if


  end subroutine ratintOutwardNormal
  


  subroutine ratintElementCentroid(elementInp, numVertices, rationalCentroid)
    class(element), intent(in) :: elementInp 
    type(ratint_t), dimension(3), intent(inout) :: rationalCentroid
    integer, intent(in) :: numVertices
    type(vertexBox), dimension(numVertices) :: vertices
    real(real64), dimension(3) :: coords
    integer :: i 


    rationalCentroid = initratint_vector()
    vertices = elementInp%getVertices()
    do i = 1, numVertices
      coords = vertices(i)%ptr%getCoordinates()
      rationalCentroid = rationalCentroid + convert_ieee(coords)
    end do 

    rationalCentroid(1) = rationalCentroid(1) / convert_int(numVertices*1_8)
    rationalCentroid(2) = rationalCentroid(2) / convert_int(numVertices*1_8)
    rationalCentroid(3) = rationalCentroid(3) / convert_int(numVertices*1_8)


  end subroutine ratintElementCentroid


  


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
  pure function newElementIntersectionTestPayload(r, u,dMax, skipBoundingBoxIntersectionTest,skipZeroFaces,currentFaceIdxs,front)&
  result(payload)
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal), intent(in)               :: dMax
    integer(shortInt), dimension(VALENCE), intent(in), optional :: currentFaceIdxs
    integer(shortInt), intent(in), optional :: front
    logical(defBool), intent(in)            :: skipBoundingBoxIntersectionTest
    logical(defBool), intent(in), optional  :: skipZeroFaces
    type(elementIntersectionTestPayload)    :: payload

    payload % r = r
    payload % u = u
    payload % dMax = dMax
    if (present(currentFaceIdxs) .and. present(front)) then
      payload % currentFaceIdxs = currentFaceIdxs
      payload%front = front
    end if
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