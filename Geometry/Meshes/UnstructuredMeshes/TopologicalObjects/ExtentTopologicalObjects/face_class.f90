module face_class
  
  use axisAlignedBoundingBox_class,  only : axisAlignedBoundingBox
  use extentTopologicalObject_inter, only : buildExtentTopologicalObjectPayload, extentTopologicalObject
  use edge_class,                    only : edgeBox
  use genericProcedures,             only : append, areEqual, crossProduct, fatalError, numToChar
  use numPrecision
  use publicObjects,                 only : intersectionTestPayload, intersectionTestResult, meshBoundaryConditionInfo, &
                                            rationalIntersectionTestPayload, resetIntersectionTestResult
  use topologicalObject_inter,       only : buildTopologicalObjectPayload, kill_super => kill, topologicalObjectBox
  use universalVariables
  use vertex_class,                  only : vertexBox
  use ratint
  
  implicit none
  private

  !!
  !!
  !!
  type, public, extends(buildExtentTopologicalObjectPayload) :: buildFacePayload
    integer(shortInt)                        :: parentIdx = 0
    integer(shortInt), dimension(N_BC_TYPES) :: boundaryConditions = [INTERNAL_TRANSPORT_BC, INTERNAL_TEMPERATURE_BC]
    logical(defBool)                         :: isBoundary = .false., testNormal = .false.
    real(defReal), dimension(N_BC_TYPES)     :: boundaryValues = ZERO
    real(defReal), dimension(3)              :: testCentroid = ZERO
    type(edgeBox), dimension(:), allocatable :: edges
  end type buildFacePayload

  !!
  !! Small, local container to store polymorphic faces in a single array.
  !!
  !! Public members:
  !!   name -> Name of the mesh.
  !!   ptr  -> Pointer to the mesh.
  !!
  type, public          :: faceBox
    type(face), pointer :: ptr => null()
  end type

  !!
  !!
  !!
  type, public                   :: orientatedFaceBox
    logical(defBool)             :: isOwner = .false.
    real(defReal), dimension(3)  :: outwardNormal = ZERO
    type(faceBox)                :: face
    type(ratint_t), dimension(3) :: ratintOutwardNormal
  end type orientatedFaceBox
  
  !! Face of an unstructured mesh. Consists of a list of vertices indices making the face up and 
  !! face-to-element connectivity information.
  !!
  !! Private members:
  !!   idx            -> Index of the face.
  !!   vertices       -> Array of vertices indices making the face up.
  !!   faceToElements -> Array listing the owner and neighbour elements for the face.
  !!   triangles      -> Array of triangles indices into which the face is decomposed.
  !!   boundaryFace   -> Is the face a boundary face?
  !!   area           -> Area of the face.
  !!   centroid       -> Vector pointing to the centroid of the face.
  !!   normal         -> Normal vector of the face.
  !!
  type, public, extends(extentTopologicalObject)          :: face
    private
    character(:), allocatable                             :: type
    integer(shortInt)                                     :: parentIdx = 0
    integer(shortInt), dimension(N_BC_TYPES)              :: boundaryConditions = [INTERNAL_TRANSPORT_BC, INTERNAL_TEMPERATURE_BC]
    integer(shortInt), dimension(:), allocatable          :: childrenIdxs
    logical(defBool)                                      :: isBoundary = .false.
    real(defReal)                                         :: area = ZERO
    real(defReal), dimension(N_BC_TYPES)                  :: boundaryValues = ZERO
    real(defReal), dimension(3)                           :: normal = ZERO
    type(ratint_t), dimension(3)                          :: ratintNormal
    type(edgeBox), dimension(:), allocatable              :: edges
    type(topologicalObjectBox), dimension(:), allocatable :: sharingElements
    type(vertexBox), dimension(:), allocatable            :: vertices
  contains
    procedure          :: addChildIdx
    procedure          :: addEdge
    procedure          :: addSharingElement
    procedure          :: addVertex
    procedure          :: build
    procedure, private :: buildComponents
    procedure          :: connectComponents
    procedure          :: distanceSquared
    procedure          :: flipDirection
    procedure          :: getArea
    procedure          :: getBoundaryCondition
    procedure          :: getBoundaryConditions
    procedure          :: getBoundaryValue
    procedure          :: getBoundaryValues
    procedure          :: getChildrenIdxs
    procedure          :: getEdgeIdxs
    procedure          :: getEdges
    procedure          :: getFaceIdx
    procedure          :: getFirstVertexCoordinates
    procedure          :: getFirstVertexRationalCoordinates
    procedure          :: getIsBoundary
    procedure          :: getNormal
    procedure          :: getRatintNormal
    procedure          :: getSharingElements
    procedure          :: getSharingFaceIdxs
    procedure          :: getSharingFaces
    procedure          :: getType
    procedure          :: getVertexIdxs
    procedure          :: getVertices
    procedure          :: intersects_BoundingBox
    procedure          :: intersects_Ray
    procedure          :: intersects_Ray_rational
    procedure          :: isPointInside
    procedure          :: isPointInside_rational
    procedure          :: isPointNearEdgeOrVertex
    procedure          :: kill
    procedure          :: setArea
    procedure          :: setBoundaryConditions
    procedure          :: setIsBoundary
    procedure          :: setNormal
    procedure          :: setVertices
  end type face

contains

  !! Subroutine 'addTriangleIdx'
  !!
  !! Basic description:
  !!   Adds the index of a triangle in the face.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the triangle.
  !!
  elemental subroutine addChildIdx(self, idx)
    class(face), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx
    
    call append(self % childrenIdxs, idx)

  end subroutine addChildIdx

  !! Subroutine 'addEdgeIdx'
  !!
  !! Basic description:
  !!   Adds the index of an edge sharing the face.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the edge sharing the face.
  !!
  subroutine addEdge(self, edge)
    class(face), intent(inout)               :: self
    type(edgeBox), intent(in)                :: edge
    type(edgeBox), dimension(:), allocatable :: tempEdges
    integer(shortInt)                        :: nEdges
    
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
  
  !! Subroutine 'addElementIdx'
  !!
  !! Basic description:
  !!   Adds the index of an element containing the face.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element containing the face.
  !!
  subroutine addSharingElement(self, box)
    class(face), intent(inout)                            :: self
    type(topologicalObjectBox), intent(in)                :: box
    integer(shortInt)                                     :: nSharingElements
    type(topologicalObjectBox), dimension(:), allocatable :: tempSharingElements
    
    if (allocated(self % sharingElements)) then
      nSharingElements = size(self % sharingElements)
      allocate(tempSharingElements(nSharingElements + 1))
      tempSharingElements(1:nSharingElements) = self % sharingElements
      tempSharingElements(nSharingElements + 1) = box
      call move_alloc(tempSharingElements, self % sharingElements)

    else
      allocate(self % sharingElements(1))
      self % sharingElements(1) = box

    end if

  end subroutine addSharingElement
  
  !! Subroutine 'addVertexIdx'
  !!
  !! Basic description:
  !!   Adds the index of a vertex in the face.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the vertex.
  !!
  subroutine addVertex(self, vertex)
    class(face), intent(inout)                 :: self
    type(vertexBox), intent(in)                :: vertex
    type(vertexBox), dimension(:), allocatable :: tempVertices
    integer(shortInt)                          :: nVertices
    
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
    class(face), intent(inout)                          :: self
    class(buildTopologicalObjectPayload), intent(inout) :: payload
    type(buildFacePayload), pointer                     :: payloadPtr
    integer(shortInt)                                   :: nVertices
    character(*), parameter                             :: here = 'init (face_class.f90)'

    ! Downcast payload to correct type.
    select type(ptr => payload)
      type is(buildFacePayload)
        payloadPtr => ptr

      class default
        call fatalError(here, 'Invalid payload type.')

    end select

    ! Catch invalid number of vertices.
    nVertices = size(payloadPtr % vertices)
    if (nVertices < 3) then
      call fatalError(here, 'A face must have at least three vertices. Has: '//numToChar(nVertices)//'.')

    elseif(nVertices == 3) then
      self % type = 'Triangle'

    else
      self % type = 'Polygon'

    end if

    ! Set everything from payload.
    self % parentIdx = payloadPtr % parentIdx
    self % isBoundary = payloadPtr % isBoundary
    self % boundaryConditions = payloadPtr % boundaryConditions
    self % boundaryValues = payloadPtr % boundaryValues
    self % vertices = payloadPtr % vertices
    self % edges = payloadPtr % edges

    ! Build components.
    call self % buildComponents(payloadPtr)

    ! Check if normal test was requested.
    if (payloadPtr % testNormal) then
      if (dot_product(payloadPtr % centroid - payloadPtr % testCentroid, self % normal) < ZERO) then
        self % vertices(1) = payloadPtr % vertices(2)
        self % vertices(2) = payloadPtr % vertices(1)
        self % normal = -self % normal
        self % ratintNormal = (-1_8) * (self % ratintNormal)

      end if

    end if

  end subroutine build

  !!
  !!
  !!
  subroutine buildComponents(self, payload)
    class(face), intent(inout)            :: self
    type(buildFacePayload), intent(inout) :: payload
    integer(shortInt)                     :: i, nVertices
    real(defReal), dimension(3, 3)        :: triangleCoordsArray
    real(defReal), dimension(3)           :: normal, sumAreasCentroid, sumNormals
    real(defReal)                         :: normalNorm, sumAreas
    character(*), parameter               :: here = 'buildComponents (face_class.f90)'

    ! First retrieve the coordinates of all the vertices in the face.
    nVertices = size(self % vertices)
    allocate(payload % allCoords(3, nVertices))
    payload % rationalCentroid = convert_int(0_longInt)
    do i = 1, nVertices
      if (.not. associated(self % vertices(i) % ptr)) &
      call fatalError(here, 'Face with index '//numToChar(self % getIdx())//' contains a null vertex pointer.')
      payload % allCoords(:, i) = self % vertices(i) % ptr % getCoordinates()
      payload % rationalCentroid = payload % rationalCentroid + self % vertices(i) % ptr % getRatintCoordinates()

    end do

    ! Now average the exact face centroid.
    payload % rationalCentroid = payload % rationalCentroid / int(nVertices, longInt)

    ! Check if the face is a triangle. If so, perform a direct computation to avoid round-off errors.
    if (nVertices == 3) then
      normal = computeTriangleNormal(payload % allCoords)
      normalNorm = norm2(normal)
      self % area = HALF * normalNorm
      payload % centroid = THIRD * sum(payload % allCoords, 2)
      self % normal = normal / normalNorm

    else
      ! Calculate the polygon's geometric centroid.
      triangleCoordsArray(:, 3) = sum(payload % allCoords, 2) / nVertices
      sumAreas = ZERO
      sumAreasCentroid = ZERO
      sumNormals = ZERO
      do i = 1, nVertices
        triangleCoordsArray(:, 1) = payload % allCoords(:, i)
        triangleCoordsArray(:, 2) = payload % allCoords(:, merge(1, i + 1, i == nVertices))

        normal = computeTriangleNormal(triangleCoordsArray)
        sumNormals = sumNormals + normal

        normalNorm = norm2(normal)
        sumAreas = sumAreas + normalNorm
        sumAreasCentroid = sumAreasCentroid + normalNorm * sum(triangleCoordsArray, 2)

      end do
      self % area = HALF * sumAreas
      payload % centroid = THIRD * sumAreasCentroid / sumAreas

      normal = computeTriangleNormal(payload % allCoords(:, 1:3))
      self % normal = normal / norm2(normal)

    end if

    self % ratintNormal = crossProduct(self % vertices(1) % ptr % getRatintCoordinates() - &
                                       self % vertices(2) % ptr % getRatintCoordinates(), &
                                       self % vertices(1) % ptr % getRatintCoordinates() - &
                                       self % vertices(3) % ptr % getRatintCoordinates())

  contains
    !!
    !!
    !!
    pure function computeTriangleNormal(array) result(n)
      real(defReal), dimension(3, 3), intent(in) :: array
      real(defReal), dimension(3)                :: n

      n = crossProduct(array(:, 2) - array(:, 1), array(:, 3) - array(:, 1))

    end function computeTriangleNormal
    
  end subroutine buildComponents

  !!
  !!
  !!
  subroutine connectComponents(self)
    class(face), target, intent(inout) :: self
    type(topologicalObjectBox)         :: box
    integer(shortInt)                  :: i

    box % ptr => self
    do i = 1, size(self % vertices)
      call self % vertices(i) % ptr % addSharingFace(box)
      call self % edges(i) % ptr % addSharingFace(box)

    end do

  end subroutine connectComponents

  !!
  !!
  !!
  function distanceSquared(self, r) result(dSquared)
    class(face), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: d, dSquared, inverseNormalSquared
    real(defReal), dimension(3)             :: diff, proj
    integer(shortInt)                       :: i

    ! First compute the distance between the point and the plane of the face.
    diff = r - self % getCentroid()
    d = dot_product(diff, self % normal)

    ! Now project the point on the plane of the face and check if the projection lies inside the face.
    inverseNormalSquared = ONE / dot_product(self % normal, self % normal)
    proj = r - self % normal * d * inverseNormalSquared

    ! If projection is inside the face, compute dSquared and return.
    if (self % isPointInside(proj)) then
      dSquared = d * d * inverseNormalSquared
      return

    end if

    ! If projection is outside the face, we need to compute the distance to each edge of the face and
    ! retain the mininum distance.
    dSquared = INF
    do i = 1, size(self % edges)
      dSquared = min(dSquared, self % edges(i) % ptr % distanceSquared(r))

    end do

  end function distanceSquared

  !!
  !!
  !!
  subroutine flipDirection(self, u)
    class(face), intent(in)                    :: self
    real(defReal), dimension(3), intent(inout) :: u

    u = u - TWO * dot_product(u, self % normal) * self % normal

  end subroutine flipDirection
  
  !! Function 'getArea'
  !!
  !! Basic description:
  !!   Returns the area of the face.
  !!
  !! Result:
  !!   area -> Area of the face.
  !!
  elemental function getArea(self) result(area)
    class(face), intent(in) :: self
    real(defReal)           :: area
    
    area = self % area

  end function getArea

  !!
  !!
  !!
  function getBoundaryCondition(self, boundaryConditionType) result(boundaryCondition)
    class(face), intent(in)       :: self
    integer(shortInt), intent(in) :: boundaryConditionType
    integer(shortInt)             :: boundaryCondition
    character(*), parameter       :: here = 'getBoundaryCondition (face_class.f90)'

    select case(boundaryConditionType)
      case(TRANSPORT_BCs)
        boundaryCondition = self % boundaryConditions(TRANSPORT_BCs)

      case(TEMPERATURE_BCs)
        boundaryCondition = self % boundaryConditions(TEMPERATURE_BCs)

      case default
        call fatalError(here, 'Invalid boundary condition type: '//numToChar(boundaryConditionType)//'.')

    end select

  end function getBoundaryCondition

  !!
  !!
  !!
  pure function getBoundaryConditions(self) result(boundaryConditions)
    class(face), intent(in)                  :: self
    integer(shortInt), dimension(N_BC_TYPES) :: boundaryConditions

    boundaryConditions = self % boundaryConditions

  end function getBoundaryConditions

  !!
  !!
  !!
  function getBoundaryValue(self, boundaryConditionType) result(boundaryValue)
    class(face), intent(in)       :: self
    integer(shortInt), intent(in) :: boundaryConditionType
    real(defReal)                 :: boundaryValue
    character(*), parameter       :: here = 'getBoundaryValue (face_class.f90)'

    select case(boundaryConditionType)
      case(TRANSPORT_BCs)
        boundaryValue = self % boundaryValues(TRANSPORT_BCs)

      case(TEMPERATURE_BCs)
        boundaryValue = self % boundaryValues(TEMPERATURE_BCs)

      case default
        call fatalError(here, 'Invalid boundary condition type: '//numToChar(boundaryConditionType)//'.')

    end select

  end function getBoundaryValue

  !!
  !!
  !!
  pure function getBoundaryValues(self) result(boundaryValues)
    class(face), intent(in)              :: self
    real(defReal), dimension(N_BC_TYPES) :: boundaryValues

    boundaryValues = self % boundaryValues

  end function getBoundaryValues

  !! Function 'getTriangleIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the triangles in the face.
  !!
  !! Result:
  !!   trianglesIdxs -> Indices of the triangles in the face.
  !!
  pure function getChildrenIdxs(self) result(childrenIdxs)
    class(face), intent(in)                      :: self
    integer(shortInt), dimension(:), allocatable :: childrenIdxs
    
    if (allocated(self % childrenIdxs)) then
      childrenIdxs = self % childrenIdxs

    else
      allocate(childrenIdxs(0))

    end if

  end function getChildrenIdxs

  !!
  !!
  !!
  pure function getEdgeIdxs(self) result(edgeIdxs)
    class(face), intent(in)                      :: self
    integer(shortInt)                            :: i, nEdges
    integer(shortInt), dimension(:), allocatable :: edgeIdxs

    if(allocated(self % edges)) then
      nEdges = size(self % edges)
      allocate(edgeIdxs(nEdges))

      do i = 1, nEdges
        edgeIdxs(i) = self % edges(i) % ptr % getIdx()

      end do

    else
      allocate(edgeIdxs(0))

    end if

  end function getEdgeIdxs

  !! Function 'getEdgeIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the edges in the face.
  !!
  !! Result:
  !!   edgeIdxs -> Indices of the edges in the face.
  !!
  function getEdges(self) result(edges)
    class(face), intent(in)                      :: self
    type(edgeBox), dimension(size(self % edges)) :: edges

    edges = self % edges

  end function getEdges
  
  !! Function 'getElementIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the elements containing the face.
  !!
  !! Result:
  !!   elementIdxs -> Indices of the elements containing the face.
  !!
  function getSharingElements(self) result(sharingElements)
    class(face), target, intent(in)                       :: self
    type(topologicalObjectBox), dimension(:), allocatable :: sharingElements
    
    if (allocated(self % sharingElements)) then
      sharingElements = self % sharingElements

    else
      allocate(sharingElements(0))

    end if

  end function getSharingElements

  !!
  !!
  !!
  pure function getSharingFaceIdxs(self) result(sharingFaceIdxs)
    class(face), intent(in)                      :: self
    integer(shortInt), dimension(:), allocatable :: sharingFaceIdxs

    allocate(sharingFaceIdxs(1))
    sharingFaceIdxs(1) = self % getIdx()

  end function getSharingFaceIdxs

  !!
  !!
  !!
  function getSharingFaces(self) result(sharingFaces)
    class(face), target, intent(in)                       :: self
    type(topologicalObjectBox), dimension(:), allocatable :: sharingFaces

    allocate(sharingFaces(1))
    sharingFaces(1) % ptr => self

  end function getSharingFaces

  !! Function 'getFaceIdx'
  !!
  !! Basic description:
  !!   Returns the index of the face from which the face originates.
  !!
  !! Result:
  !!   faceIdx -> Index of the face from which the face originates.
  !!
  elemental function getFaceIdx(self) result(faceIdx)
    class(face), intent(in) :: self
    integer(shortInt)       :: faceIdx
    
    faceIdx = self % parentIdx

  end function getFaceIdx

  !!
  !!
  !!
  pure function getFirstVertexCoordinates(self) result(firstVertexCoordinates)
    class(face), intent(in)     :: self
    real(defReal), dimension(3) :: firstVertexCoordinates

    firstVertexCoordinates = self % vertices(1) % ptr % getCoordinates()

  end function getFirstVertexCoordinates

  !!
  !!
  !!
  pure function getFirstVertexRationalCoordinates(self) result(firstVertexRationalCoordinates)
    class(face), intent(in)      :: self
    type(ratint_t), dimension(3) :: firstVertexRationalCoordinates

    firstVertexRationalCoordinates = self % vertices(1) % ptr % getRatintCoordinates()

  end function getFirstVertexRationalCoordinates

  !! Function 'getIsBoundary'
  !!
  !! Basic description:
  !!   Returns .true. if the face is a boundary face.
  !!
  !! Result:
  !!   isBoundary -> .true. if the face is a boundary face.
  !!
  elemental function getIsBoundary(self) result(isBoundary)
    class(face), intent(in) :: self
    logical(defBool)        :: isBoundary

    isBoundary = self % isBoundary

  end function getIsBoundary 
  
  !! Function 'getNormal'
  !!
  !! Basic description:
  !!   Returns the normal vector of the face.
  !!
  !! Result:
  !!   normal -> Normal vector of the face.
  !!
  pure function getNormal(self, idx) result(normal)
    class(face), intent(in)                 :: self
    real(defReal), dimension(3)             :: normal
    integer(shortInt), intent(in), optional :: idx
    
    normal = self % normal
    
    if (.not. present(idx)) return
    if (idx < 0) normal = -normal

  end function getNormal



  pure function getRatintNormal(self, idx) result(ratintNormal)
    class(face), intent(in)                 :: self
    type(ratint_t), dimension(3)             :: ratintNormal
    integer(shortInt), intent(in), optional :: idx
    
    ratintNormal = self % ratintNormal
    
    if (.not. present(idx)) return
    if (idx < 0) ratintNormal = (-1_8)*ratintNormal

  end function getRatintNormal

  !! Function 'getTriangleIdxs'
  !!
  !! Basic description:
  !!   Returns the type of the face.
  !!
  !! Result:
  !!   type -> Type of the face.
  !!
  pure function getType(self) result(type)
    class(face), intent(in)      :: self
    character(len(self % type)) :: type
    
    type = self % type

  end function getType

  !!
  !!
  !!
  pure function getVertexIdxs(self) result(vertexidxs)
    class(face), intent(in)                      :: self
    integer(shortInt)                            :: i, nVertices
    integer(shortInt), dimension(:), allocatable :: vertexidxs

    if(allocated(self % vertices)) then
      nVertices = size(self % vertices)
      allocate(vertexidxs(nVertices))

      do i = 1, nVertices
        vertexidxs(i) = self % vertices(i) % ptr % getIdx()

      end do

    else
      allocate(vertexidxs(0))

    end if

  end function getVertexIdxs
  
  !! Function 'getVertexIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in the face.
  !!
  !! Result:
  !!   vertexIdxs -> Array of vertex indices in the face.
  !!
  function getVertices(self) result(vertices)
    class(face), intent(in)                           :: self
    type(vertexBox), dimension(size(self % vertices)) :: vertices
    
    vertices = self % vertices

  end function getVertices

  !!
  !!
  !!
  pure subroutine intersects_BoundingBox(self, boundingBox, doesIt)
    class(face), intent(in)                            :: self
    type(axisAlignedBoundingBox), intent(in)           :: boundingBox
    logical(defBool), intent(out)                      :: doesIt
    real(defReal), dimension(3)                        :: boundingBoxCentre, halfwidths, axis, edgeVector, boxAxis
    real(defReal), dimension(3, size(self % vertices)) :: centredVertexCoords
    integer(shortInt)                                  :: i, j, nVertices

    ! Initialise doesIt = .false., retrieve the centre and halfwidths of the boundingBox.
    doesIt = .false.

    ! First check if the bounding box intersects the face's bounding box.
    if (.not. self % intersectsBoundingBox(boundingBox)) return

    boundingBoxCentre = boundingBox % getCentre()
    halfwidths = boundingBox % getHalfwidths()
    ! Offset the coordinates of the face vertices with respect to the box centre.
    nVertices = size(self % vertices)
    do i = 1, nVertices
      centredVertexCoords(:, i) = self % vertices(i) % ptr % getCoordinates() - boundingBoxCentre

    end do

    ! First test for intersection along the three bounding box axes.
    do i = 1, 3
      axis = ZERO
      axis(i) = ONE
      if (.not. overlaps(halfwidths, centredVertexCoords, axis, nVertices)) return

    end do

    ! Now test the face's normal vector.
    if (.not. overlaps(halfwidths, centredVertexCoords, self % normal, nVertices)) return

    ! Finally, test cross products between the face's edges and the bounding box's edges.
    do i = 1, size(self % edges)
      edgeVector = self % edges(i) % ptr % getEdgeVector()
      do j = 1, 3
        boxAxis = ZERO
        boxAxis(j) = ONE
        axis = crossProduct(edgeVector, boxAxis)
        if (.not. overlaps(halfwidths, centredVertexCoords, axis, nVertices)) return

      end do

    end do

    ! If reached here, the face and the bounding box intersect so update doesIt = .true.
    doesIt = .true.

  contains
    !!
    !!
    !!
    pure function overlaps(h, coords, ax, n) result(isOverlapping)
      real(defReal), dimension(3), intent(in)                        :: h, ax
      real(defReal), dimension(3, size(self % vertices)), intent(in) :: coords
      integer(shortInt), intent(in)                                  :: n
      logical(defBool)                                               :: isOverlapping
      real(defReal)                                                  :: radius, minProjection, maxProjection, d
      integer(shortInt)                                              :: k

      ! Compute the box radius.
      radius = dot_product(h, abs(ax))

      ! Compute d and initialise minProjection and maxProjection.
      d = dot_product(coords(:, 1), ax)
      minProjection = d
      maxProjection = d

      do k = 2, n
        d = dot_product(coords(:, k), ax)
        minProjection = min(minProjection, d)
        maxProjection = max(maxProjection, d)

      end do

      ! Check if overlap between projections.
      isOverlapping = minProjection <= radius .and. -radius <= maxProjection

    end function overlaps

  end subroutine intersects_BoundingBox

  !!
  !!
  !!
  subroutine intersects_Ray(self, payload, result)
    class(face), intent(in)                    :: self
    class(intersectionTestPayload), intent(in) :: payload
    class(intersectionTestResult), intent(out) :: result
    real(defReal)                              :: denominator, numerator, t
    real(defReal), dimension(3)                :: firstVertexCoordinates, rIntersection

    ! First check if ray intersects the face's bounding box and return early if not.
    !call intersects_Ray_super(self, payload, result)
    !if(.not. result % intersects) return

    ! Retrieve the coordinates of the first vertex in the face.
    firstVertexCoordinates = self % getFirstVertexCoordinates()

    ! Compute numerator and denominator.
    numerator = dot_product(firstVertexCoordinates - payload % r, self % normal)
    denominator = dot_product(self % normal, payload % u)
    if(areEqual(denominator, ZERO)) then
      ! If denominator is nearly equal to zero, and the ray lies in the plane of the face, escalate to exact computation.
      if(areEqual(numerator, ZERO)) result % needsRescue = .true.
      return

    end if
    
    ! Compute distance along the ray to intersection.
    t = numerator / denominator

    ! Escalate to exact computation if the segment either starts or ends within tolerance of the face.
    if(areEqual(t, ZERO) .or. areEqual(t, payload % dMax)) result % needsRescue = .true.

    ! Return early if intersection is not possible and the result is unambiguous.
    if((t < ZERO .or. payload % dMax < t) .and. .not. result % needsRescue) return

    ! Compute intersection point.
    rIntersection = payload % r + payload % u * t

    ! Escalate to exact computation if the intersection point lands within tolerance of any edges or vertices.
    if(self % isPointNearEdgeOrVertex(rIntersection)) result % needsRescue = .true.

    ! Check if the intersection point coordinates are inside the face.
    if(self % isPointInside(rIntersection)) then
      result % intersects = .true.
      result % d = t

    end if

  end subroutine intersects_Ray

  !!
  !!
  !!
  subroutine intersects_Ray_rational(self, payload, result)
    class(face), intent(in)                            :: self
    class(rationalIntersectionTestPayload), intent(in) :: payload
    class(intersectionTestResult), intent(out)         :: result
    type(ratint_t)                                     :: denominator, t
    type(ratint_t), dimension(3)                       :: firstVertexCoordinates

    ! Retrieve the coordinates of the first vertex in the face.
    firstVertexCoordinates = self % getFirstVertexRationalCoordinates()

    ! Compute denominator.
    denominator = dot_product(self % ratintNormal, payload % u)
    if(isZero(denominator)) return

    ! Compute distance along the ray to intersection.
    t = dot_product(firstVertexCoordinates - payload % r, self % ratintNormal) / denominator

    ! Return early if intersection is not possible.
    if((convert_int(0_longInt) > t .or. t > payload % dMax)) return

    ! Check if the intersection point coordinates are inside the face.
    if(self % isPointInside_rational(payload % r + payload % u * t)) then
      result % intersects = .true.
      result % d_rational = t

    end if

  end subroutine intersects_Ray_rational

  !!
  !!
  !!
  function isPointInside(self, r) result(isIt)
    class(face), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    logical(defBool)                        :: isIt
    integer(shortInt)                       :: i, nVertices
    real(defReal)                           :: dotProduct
    real(defReal), dimension(3)             :: firstVertexCoords, nextVertexCoords, vertexCoords

    ! Initialise isIt = .false. and compute the number of vertices in the face.
    isIt = .false.
    nVertices = size(self % vertices)

    ! Retrieve the coordinates of the first vertex and initialise vertexCoords = firstVertexCoords.
    firstVertexCoords = self % vertices(1) % ptr % getCoordinates()
    vertexCoords = firstVertexCoords

    ! Loop through all the edges in the face and check if the point lies on the same side
    ! of each edge (note: this assumes a consistent vertex numbering).
    do i = 1, nVertices
      if(i < nVertices) then
        nextVertexCoords = self % vertices(i + 1) % ptr % getCoordinates()

      else
        nextVertexCoords = firstVertexCoords

      end if
      dotProduct = dot_product(self % normal, crossProduct(nextVertexCoords - vertexCoords, r - vertexCoords))
      if (dotProduct < ZERO) return

      ! Update vertexCoords.
      vertexCoords = nextVertexCoords

    end do

    ! If reached here, the point is inside the face.
    isIt = .true.

  end function isPointInside

  !!
  !!
  !!
  function isPointInside_rational(self, r) result(isIt)
    class(face), intent(in)                  :: self
    type(ratint_t), dimension(3), intent(in) :: r
    logical(defBool)                         :: isIt
    integer(shortInt)                        :: i, nVertices
    type(ratint_t)                           :: dotProduct, ZERO_rational
    type(ratint_t), dimension(3)             :: firstVertexCoords, nextVertexCoords, vertexCoords

    ! Pre-compute ZERO_rational.
    ZERO_rational = convert_int(0_longInt)

    ! Initialise isIt = .false. and compute the number of vertices in the face.
    isIt = .false.
    nVertices = size(self % vertices)

    ! Retrieve the coordinates of the first vertex and initialise vertexCoords = firstVertexCoords.
    firstVertexCoords = self % vertices(1) % ptr % getRatintCoordinates()
    vertexCoords = firstVertexCoords

    ! Loop through all the edges in the face and check if the point lies on the same side
    ! of each edge (note: this assumes a consistent vertex numbering).
    do i = 1, nVertices
      if(i < nVertices) then
        nextVertexCoords = self % vertices(i + 1) % ptr % getRatintCoordinates()

      else
        nextVertexCoords = firstVertexCoords

      end if
      dotProduct = dot_product(self % ratintNormal, crossProduct(nextVertexCoords - vertexCoords, r - vertexCoords))
      if(ZERO_rational > dotProduct) return

      ! Update vertexCoords.
      vertexCoords = nextVertexCoords

    end do

    ! If reached here, the point is inside the face.
    isIt = .true.

  end function isPointInside_rational

  !!
  !!
  !!
  pure function isPointNearEdgeOrVertex(self, r) result(isIt)
    class(face), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    integer(shortInt)                       :: i
    logical(defBool)                        :: isIt

    ! Initialise isIt = .true.
    isIt = .true.

    ! Compute distance to each vertex and immediately return if point is within tolerance to any of them.
    do i = 1, size(self % vertices) 
      if(areEqual(sqrt(self % vertices(i) % ptr % distanceSquared(r)), ZERO)) return

    end do

    ! Compute distance to each edge and immediately return if point is within tolerance to any of them.
    do i = 1, size(self % edges) 
      if(areEqual(sqrt(self % edges(i) % ptr % distanceSquared(r)), ZERO)) return

    end do

    ! If reached here, the point is not within epsilon tolerance of any edges or vertices.
    isIt = .false.

  end function isPointNearEdgeOrVertex
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(face), intent(inout) :: self
    integer(shortInt)          :: i
    
    ! Superclass.
    call kill_super(self)

    ! Local.
    self % parentIdx = 0
    self % isBoundary = .false.
    self % area = ZERO
    self % boundaryValues = ZERO
    self % normal = ZERO
    self % boundaryConditions = [INTERNAL_TRANSPORT_BC, INTERNAL_TEMPERATURE_BC]
    if (allocated(self % childrenIdxs)) deallocate(self % childrenIdxs)
    if (allocated(self % edges)) then
      do i = 1, size(self % edges)
        nullify(self % edges(i) % ptr)

      end do
      deallocate(self % edges)

    end if

    if (allocated(self % sharingElements)) then
      do i = 1, size(self % sharingElements)
        nullify(self % sharingElements(i) % ptr)

      end do
      deallocate(self % sharingElements)

    end if

    if (allocated(self % vertices)) then
      do i = 1, size(self % vertices)
        nullify(self % vertices(i) % ptr)

      end do
      deallocate(self % vertices)

    end if

  end subroutine kill

  !! Subroutine 'setArea'
  !!
  !! Basic description:
  !!   Sets the area of the face.
  !!
  !! Arguments:
  !!   area [in] -> Area of the face.
  !!
  elemental subroutine setArea(self, area)
    class(face), intent(inout) :: self
    real(defReal), intent(in)  :: area

    self % area = area

  end subroutine setArea

  !!
  !!
  !!
  elemental subroutine setBoundaryConditions(self, info)
    class(face), intent(inout)                  :: self
    type(meshBoundaryConditionInfo), intent(in) :: info

    self % boundaryConditions = info % boundaryConditions
    self % boundaryValues = info % boundaryValues

  end subroutine setBoundaryConditions
  
  !! Subroutine 'setBoundaryFace'
  !!
  !! Basic description:
  !!   Sets the face as a boundary face.
  !!
  elemental subroutine setIsBoundary(self, isBoundary)
    class(face), intent(inout)   :: self
    logical(defBool), intent(in) :: isBoundary
    
    self % isBoundary = isBoundary

  end subroutine setIsBoundary

  !! Subroutine 'setNormal'
  !!
  !! Basic description:
  !!   Sets the normal vector of the face.
  !!
  !! Arguments:
  !!   normal [in] -> 3-D coordinates of the normal vector of the face.
  !!
  pure subroutine setNormal(self, normal)
    class(face), intent(inout)               :: self
    real(defReal), dimension(3), intent(in)  :: normal

    self % normal = normal

  end subroutine setNormal

  !! Subroutine 'setVertexIdxs'
  !!
  !! Basic description:
  !!   Sets the indices of the vertices in the face.
  !!
  !! Arguments:
  !!   vertexIdxs [in] -> Indices of the vertices in the face.
  !!
  subroutine setVertices(self, vertices)
    class(face), intent(inout)                :: self
    type(vertexBox), dimension(:), intent(in) :: vertices

    self % vertices = vertices

  end subroutine setVertices
  
end module face_class