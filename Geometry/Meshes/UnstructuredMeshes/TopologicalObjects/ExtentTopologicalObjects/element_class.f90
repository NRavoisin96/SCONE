module element_class

  use axisAlignedBoundingBox_class,  only : axisAlignedBoundingBox
  use extentTopologicalObject_inter, only : buildExtentTopologicalObjectPayload, extentTopologicalObject
  use edge_class,                    only : edgeBox
  use errors_mod,                    only : fatalError
  use face_class,                    only : faceBox, orientatedFaceBox, face
  use genericProcedures,             only : append, areEqual, crossProduct, findCommon, numToChar
  use limb_class
  use numPrecision
  use publicObjects,                 only : basicElementInfo, intersectionTestPayload, intersectionTestResult, &
                                            resetIntersectionTestResult
  use ratint
  use RNG_class,                     only : RNG
  use topologicalObject_inter,       only : buildTopologicalObjectPayload, kill_super => kill, topologicalObject, &
                                            topologicalObjectBox
  use universalVariables,            only : FOURTH, INSIDE_ELEMENT, INF, NUDGE, ON_BOUNDARY_ELEMENT, ONE, OUTSIDE_ELEMENT, &
                                            SIXTH, ZERO, VALENCE
  use vertex_class,                  only : vertexBox
  
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
    character(:), allocatable                          :: type
    integer(shortInt)                                  :: parentIdx = 0, localId = 0
    integer(shortInt), dimension(:), allocatable       :: childrenIdxs
    logical(defBool)                                   :: isConvex = .false.
    real(defReal)                                      :: volume = ZERO
    type(edgeBox), dimension(:), allocatable           :: edges
    type(orientatedFaceBox), dimension(:), allocatable :: orientatedFaces
    type(vertexBox), dimension(:), allocatable         :: vertices
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
    procedure          :: computeExactIntersection
    procedure          :: distanceSquared
    generic            :: entersThroughFaces => entersThroughFaces_defReal, entersThroughFaces_rational
    procedure, private :: entersThroughFaces_defReal
    procedure, private :: entersThroughFaces_rational
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
    generic            :: isPointInside => isPointInside_defReal, isPointInside_rational
    procedure, private :: isPointInside_defReal
    procedure, private :: isPointInside_rational
    procedure          :: kill
    procedure          :: minimumDistance
    procedure          :: pushFromBoundary
    procedure          :: sampleInitialPosition
  end type element

  !!
  !!
  !!
  type, public        :: inclusionTestResult
    integer(shortInt) :: status = INSIDE_ELEMENT
  end type inclusionTestResult

  !!
  !!
  !!
  type, public, extends(intersectionTestPayload) :: elementIntersectionTestPayload
    logical(defBool) :: excludeZeroFaces = .false., skipBoundingBoxIntersectionTest = .false.
  end type elementIntersectionTestPayload

  !!
  !!
  !!
  type, public, extends(intersectionTestResult) :: elementIntersectionTestResult
    integer(shortInt)                     :: front = 0
    integer(shortInt), dimension(VALENCE) :: currentFaceIdxs = 0
    real(defReal), dimension(3)           :: intersectionPt = ZERO
    type(faceBox)                         :: intersectedFace
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
    payload % rationalCentroid = convert_int(0_longInt)
    do i = 1, nVertices
      if (.not. associated(self % vertices(i) % ptr)) call fatalError(here, 'Element contains a null vertex pointer.')
      payload % allCoords(:, i) = self % vertices(i) % ptr % getCoordinates()
      payload % rationalCentroid = payload % rationalCentroid + self % vertices(i) % ptr % getRatintCoordinates()

    end do

    ! Now average exact centroid.
    payload % rationalCentroid = payload % rationalCentroid / int(nVertices, longInt)

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

  !!!NOTE, THROUGH ANOTHER ISSUE, SOMETHING IS WRONG HERE FOR BOUNDARY EXITS?? MAYBE DIRECTION??
  subroutine computeExactIntersection(self, excludedFaceIdxs, dMax, r, rEnd, u, res, minLambda)
    class(element), intent(in)                                       :: self
    integer(shortInt), dimension(:), intent(in)                      :: excludedFaceIdxs
    real(defReal), intent(in)                                        :: dMax
    real(defReal), dimension(3), intent(in)                          :: r, rEnd, u
    type(elementIntersectionTestResult), intent(inout)               :: res
    real(defReal), intent(out)                                       :: minLambda
    integer(shortInt)                                                :: i, nPotentiallyIntersectedFaces
    integer(shortInt), dimension(size(self % orientatedFaces))       :: potentiallyIntersectedFaceIdxs
    type(faceBox)                                                    :: tempFace
    type(ratint_t)                                                   :: dotProduct_ratint, faceLambda, minLambda_ratint, &
                                                                        ONE_ratint, ZERO_ratint
    type(ratint_t), dimension(3)                                     :: firstVertexCoordinates_ratint, outwardNormal_ratint, &
                                                                        r_ratint, rEnd_ratint
    type(ratint_t), dimension(size(self % orientatedFaces))          :: potentiallyIntersectedFaceLambdas
    type(vertexBox), dimension(:), allocatable                       :: faceVertices

    ! Reset result % front = 0.
    res % front = 0

    ! Pre-compute ONE and ZERO represented as a ratint.
    ONE_ratint = convert_int(1_longInt)
    ZERO_ratint = convert_int(0_longInt)
    
    nPotentiallyIntersectedFaces = 0

    r_ratint = convert_ieee(r)
    rEnd_ratint = r_ratint + convert_ieee(u) * convert_ieee(dMax)
    minLambda_ratint = def_ratint_large()

    ! Loops through all the faces in the element.
    do i = 1, size(self % orientatedFaces)
      ! Check if the current face is excluded from the computation and cycle if so.
      if(any(excludedFaceIdxs == self % orientatedFaces(i) % face % ptr % getIdx())) cycle

      ! Retrieve the coordinates of the first vertex in the face.
      faceVertices = self % orientatedFaces(i) % face % ptr % getVertices()
      firstVertexCoordinates_ratint = faceVertices(1) % ptr % getRatintCoordinates()

      ! Retrieve the outward normal of the current face then compute lambda.
      outwardNormal_ratint = self % orientatedFaces(i) % ratintOutwardNormal
      
      ! Compute denominator and cycle to next face if it is ZERO.
      dotProduct_ratint = dot_product(rEnd_ratint - r_ratint, outwardNormal_ratint)
      if(isZero(dotProduct_ratint)) cycle

      ! Compute lambda for the current face.
      faceLambda = dot_product(firstVertexCoordinates_ratint - r_ratint, outwardNormal_ratint) / dotProduct_ratint
      
      ! Update the array of faces intersected by the ray.
      if(faceLambda >= ZERO_ratint .and. ONE_ratint >= faceLambda) then
        nPotentiallyIntersectedFaces = nPotentiallyIntersectedFaces + 1
        potentiallyIntersectedFaceIdxs(nPotentiallyIntersectedFaces) = self % orientatedFaces(i) % face % ptr % getIdx()
        potentiallyIntersectedFaceLambdas(nPotentiallyIntersectedFaces) = faceLambda

        ! Update intersected face.
        if(minLambda_ratint > faceLambda) then
          minLambda_ratint = faceLambda
          tempFace = self % orientatedFaces(i) % face

        end if

      end if

    end do

    if(associated(tempFace % ptr)) then
      res % intersectedFace = tempFace
      do i = 1, nPotentiallyIntersectedFaces
        if(minLambda_ratint == potentiallyIntersectedFaceLambdas(i)) then
          res % front = res % front + 1
          res % currentFaceIdxs(res % front) = potentiallyIntersectedFaceIdxs(i)

        end if

      end do
      ! Update minLambda here.
      minLambda = evaluate(minLambda_ratint)

    else
      ! Return infinite distance.
      minLambda = INF

    end if
 
  end subroutine computeExactIntersection

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
  pure function entersThroughFaces_defReal(self, faceIdxs, u) result(doesIt)
    class(element), intent(in)                  :: self
    integer(shortInt), dimension(:), intent(in) :: faceIdxs
    real(defReal), dimension(3), intent(in)     :: u
    integer(shortInt)                           :: i
    logical(defBool)                            :: doesIt
    real(defReal)                               :: dotProduct

    ! Initialise doesIt = .false.
    doesIt = .false.

    do i = 1, size(self % orientatedFaces)
      ! Skip current face if it is not in the set of faces to be tested against.
      if(.not. any(faceIdxs == self % orientatedFaces(i) % face % ptr % getIdx())) cycle

      ! Compute dot product between face outward normal and direction.
      dotProduct = dot_product(u, self % orientatedFaces(i) % outwardNormal)

      if(areEqual(dotProduct, ZERO)) then
        ! If dot product is below tolerance, call exact procedure and return.
        doesIt = self % entersThroughFaces(faceIdxs, convert_ieee(u))
        return

      elseif(ZERO < dotProduct) then
        ! If dot product is positive, return immediately.
        return

      end if

    end do

    ! If reached here, update doesIt = .true.
    doesIt = .true.

  end function entersThroughFaces_defReal

  !!
  !!
  !!
  pure function entersThroughFaces_rational(self, faceIdxs, u) result(doesIt)
    class(element), intent(in)                  :: self
    integer(shortInt), dimension(:), intent(in) :: faceIdxs
    type(ratint_t), dimension(3), intent(in)    :: u
    integer(shortInt)                           :: i
    logical(defBool)                            :: doesIt
    type(ratint_t)                              :: dotProduct, ZERO_rational

    ! Initialise doesIt = .false. and precompute ZERO_rational.
    doesIt = .false.
    ZERO_rational = convert_int(0_longInt)

    do i = 1, size(self % orientatedFaces)
      ! Skip current face if it is not in the set of faces to be tested against.
      if(.not. any(faceIdxs == self % orientatedFaces(i) % face % ptr % getIdx())) cycle

      ! Compute dot product between face outward normal and direction.
      dotProduct = dot_product(u, self % orientatedFaces(i) % ratintOutwardNormal)

      ! If dot product is positive, return immediately.
      if(dotProduct > ZERO_rational) return

    end do

    ! If reached here, update doesIt = .true.
    doesIt = .true.

  end function entersThroughFaces_rational

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
      call self % orientatedFaces(i) % face % ptr % intersects_BoundingBox(boundingBox, doesIt)
      if (doesIt) return

    end do

  end subroutine intersects_BoundingBox

  !!
  !!
  !!
  subroutine intersects_Ray(self, payload, result)
    class(element), intent(in)                                       :: self
    class(intersectionTestPayload), intent(in)                       :: payload
    class(intersectionTestResult), intent(inout)                     :: result
    integer(shortInt)                                                :: i, faceArrayFront
    logical(defBool)                                                 :: needsRescue
    real(defReal)                                                    :: dotProduct, faceLambda, minLambda, newMinLambda
    real(defReal), dimension(3)                                      :: centroid, firstVertexCoordinates, outwardNormal, &
                                                                        rEnd
    real(defReal), dimension(size(self % orientatedFaces))           :: faceArrayLambdas
    type(elementIntersectionTestPayload), pointer                    :: payloadPtr
    type(elementIntersectionTestResult), pointer                     :: resultPtr
    type(faceBox)                                                    :: tempFace
    type(orientatedFaceBox), dimension(size(self % orientatedFaces)) :: faceArray
    type(vertexBox), dimension(:), allocatable                       :: faceVertices
    character(*), parameter                                          :: HERE = 'intersects_Ray (element_class.f90)'

    ! Initialise number of intersected faces to 0 and needsRescue = .false.
    faceArrayFront = 0
    needsRescue = .false.

    ! Downcast payload to correct type.
    select type(ptr => payload)
      type is(elementIntersectionTestPayload)
        payloadPtr => ptr

      class default
        call fatalError(HERE, 'Invalid payload type.')

    end select

    ! Allocate result to correct return type then associate pointer.
    select type(ptr => result)
      type is(elementIntersectionTestResult)
        resultPtr => ptr
        call resetElementIntersectionTestResult(resultPtr)

      class default
        ! Should never happen.
        call fatalError(HERE, 'Failed to downcast result.')

    end select

    ! Check if ray originates from inside the element (skip bounding box intersection in this case.)
    if(payloadPtr % skipBoundingBoxIntersectionTest) then

      ! Compute end position of the ray and initialise resultPtr % front = 0.
      rEnd = payload % r + payload % u * payload % dMax
      resultPtr % front = 0
      
      ! Retrieve element's centroid then loop over all faces in the element.
      centroid = self % getCentroid()
      minLambda = INF

      do i = 1, size(self % orientatedFaces)
        ! Check if the current face is excluded from the computation and cycle if so.
        if(0 < payload % front) then
          if(any(payload % currentFaceIdxs(1:payload % front) == self % orientatedFaces(i) % face % ptr % getIdx())) cycle

        end if
  
        ! Retrieve the centroid and outward normal vector of the current face.
        faceVertices = self % orientatedFaces(i) % face % ptr % getVertices()
        firstVertexCoordinates = faceVertices(1) % ptr % getCoordinates()
        outwardNormal = self % orientatedFaces(i) % outwardNormal

        ! Check if the starting point or end point of the ray are within tolerance of the face's plane and 
        ! set needsRescue = .true. if so.
        if(areEqual(dot_product(payload % r - firstVertexCoordinates, outwardNormal), ZERO) .or. &
           areEqual(dot_product(rEnd - firstVertexCoordinates, outwardNormal), ZERO)) needsRescue = .true.

        ! Compute denominator and cycle to next face if it is ZERO.
        dotProduct = dot_product(rEnd - payload % r, outwardNormal)
        if(areEqual(dotProduct, ZERO)) cycle

        ! Compute lambda for the current face.
        faceLambda = dot_product(firstVertexCoordinates - payload % r, outwardNormal) / dotProduct

        ! Update the array of faces intersected by the ray.
        if(ZERO <= faceLambda .and. faceLambda <= ONE) then
          faceArrayFront = faceArrayFront + 1
          faceArray(faceArrayFront) = self % orientatedFaces(i) 
          faceArrayLambdas(faceArrayFront) = faceLambda

          ! Update intersected face.
          if(faceLambda < minLambda) then
            minLambda = faceLambda
            tempFace = self % orientatedFaces(i) % face

          end if

        end if

      end do

      if(associated(tempFace % ptr)) then
        !!! code added here: check for distance less than epsilon to a vertex or edge
        if(tempFace % ptr % isPointNearEdgeOrVertex(payloadPtr % r + minLambda * (rEnd - payloadPtr % r))) &
          needsRescue = .true.

      end if

      ! If rescue is needed, compute intersection point exactly.
      if(needsRescue) then
        ! If intersection point lands within epsilon tolerance of an edge or a vertex of the intersected face, begin rescue operation.
        call self % computeExactIntersection(payloadPtr % currentFaceIdxs(1:payload % front), payload % dMax, &
                                             payload % r, rEnd, payloadPtr % u, resultPtr, newMinLambda)

        if(associated(resultPtr % intersectedFace % ptr) ) then
          resultPtr % intersects = .true.
          resultPtr % intersectionPt = payloadPtr % r + min(ONE, max(ZERO, newMinLambda)) * (rEnd - payloadPtr % r)
          resultPtr % d = norm2(min(ONE, max(ZERO, newMinLambda)) * (rEnd - payloadPtr % r))

        end if

      elseif(associated(tempFace % ptr)) then
        do i = 1, faceArrayFront
          if(areEqual(minLambda, faceArrayLambdas(i))) then
            resultPtr % front = resultPtr % front + 1
            resultPtr % currentFaceIdxs(resultPtr % front) = faceArray(i) % face % ptr % getIdx()

          end if

        end do
        resultPtr % intersectedFace = tempFace
        resultPtr % intersects = .true.
        resultPtr % intersectionPt = payload % r + min(ONE, max(ZERO, minLambda)) * (rEnd - payload % r)
        resultPtr % d = norm2(min(ONE, max(ZERO, minLambda)) * (rEnd - payload % r))

      end if

    else
      ! Call fatalError for now.
      call fatalError(HERE, 'Unsupported procedure.')

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
  pure function isPointInside_defReal(self, r) result(result)
    class(element), intent(in)                 :: self
    real(defReal), dimension(3), intent(in)    :: r
    type(inclusionTestResult)                  :: result
    integer(shortInt)                          :: i
    real(defReal)                              :: dotProduct
    
    ! Initialise result % status = INSIDE_ELEMENT then loop over all element faces.
    result % status = INSIDE_ELEMENT
    do i = 1, size(self % orientatedFaces)
      ! Make a vector going from the coordinates to the face's first vertex and perform the dot
      ! product between this vector and the face's normal vector.
      dotProduct = dot_product(self % orientatedFaces(i) % face % ptr % getFirstVertexCoordinates() - r, &
                               self % orientatedFaces(i) % outwardNormal)

      ! Check if the point is effectively on the plane of this face.
      if(areEqual(dotProduct, ZERO)) then
        ! If the point is within tolerance of the plane of this face, call exact procedure and return.
        result = self % isPointInside(convert_ieee(r))
        return

      elseif(dotProduct < ZERO) then
        ! If dot product is negative, the query point is outside the element so return immediately.
        result % status = OUTSIDE_ELEMENT
        return

      end if

    end do

  end function isPointInside_defReal

  !!
  !!
  !!
  pure function isPointInside_rational(self, r) result(res)
    class(element), intent(in)                 :: self
    type(ratint_t), dimension(3), intent(in)   :: r
    type(inclusionTestResult)                  :: res
    integer(shortInt)                          :: i
    logical(defBool)                           :: isOnBoundary
    type(ratint_t)                             :: dotProduct, ZERO_rational

    ! Initialise isOnBoundary = .false. and pre-compute ZERO_rational.
    isOnBoundary = .false.
    ZERO_rational = convert_int(0_longInt)

    ! Initialise res % status = INSIDE_ELEMENT then loop over all element faces.
    res % status = INSIDE_ELEMENT
    do i = 1, size(self % orientatedFaces)
      ! Make a vector going from the coordinates to the face's first vertex and perform the dot
      ! product between this vector and the face's normal vector.
      dotProduct = dot_product(self % orientatedFaces(i) % face % ptr % getFirstVertexRationalCoordinates() - r, &
                               self % orientatedFaces(i) % ratintOutwardNormal)

      ! Check if the point is effectively on the plane of this face.
      if(isZero(dotProduct)) then
        isOnBoundary = .true.

      elseif(ZERO_rational > dotProduct) then
        ! If dot product is negative, the query point is outside the element so return immediately.
        res % status = OUTSIDE_ELEMENT
        return

      end if

    end do

    ! If query point is on boundary, flag it here.
    if(isOnBoundary) res % status = ON_BOUNDARY_ELEMENT

  end function isPointInside_rational
  
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
  pure function newElementIntersectionTestPayload(r, u, dMax, skipBoundingBoxIntersectionTest, front, currentFaceIdxs, &
                                                  skipZeroFaces) result(payload)
    real(defReal), dimension(3), intent(in)                     :: r, u
    real(defReal), intent(in)                                   :: dMax
    logical(defBool), intent(in)                                :: skipBoundingBoxIntersectionTest
    integer(shortInt), intent(in), optional                     :: front
    integer(shortInt), dimension(VALENCE), intent(in), optional :: currentFaceIdxs
    logical(defBool), intent(in), optional                      :: skipZeroFaces
    type(elementIntersectionTestPayload)                        :: payload

    payload % r = r
    payload % u = u
    payload % dMax = dMax
    payload % skipBoundingBoxIntersectionTest = skipBoundingBoxIntersectionTest
    if(present(front)) payload % front = front
    if(present(currentFaceIdxs)) payload % currentFaceIdxs = currentFaceIdxs
    if(present(skipZeroFaces)) payload % excludeZeroFaces = skipZeroFaces

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