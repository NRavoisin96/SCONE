module face_class
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use coord_class,                  only : coord
  use edge_class,                   only : edgeBox
  use genericProcedures,            only : append, areEqual, crossProduct, fatalError, numToChar
  use topologicalObject_inter,      only : topologicalObject, kill_super => kill
  use numPrecision
  use universalVariables,           only : HALF, INF, ONE, SURF_TOL, THIRD, ZERO
  use vertex_class,                 only : vertexBox
  
  implicit none
  private

  !!
  !!
  !!
  type, public :: buildFaceInfo
    integer(shortInt)                          :: idx = 0, parentIdx = 0
    logical(defBool)                           :: isBoundary = .false., testNormal = .false.
    type(vertexBox), dimension(:), allocatable :: vertices
    type(edgeBox), dimension(:), allocatable   :: edges
    real(defReal), dimension(3)                :: testCentroid = ZERO
  end type buildFaceInfo

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
  type, public                  :: orientatedFaceBox
    type(faceBox)               :: face
    logical(defBool)            :: isOwner = .false.
    real(defReal), dimension(3) :: outwardNormal = ZERO
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
  type, public, extends(topologicalObject)       :: face
    private
    integer(shortInt)                            :: parentIdx = 0
    type(edgeBox), dimension(:), allocatable     :: edges
    type(vertexBox), dimension(:), allocatable   :: vertices
    integer(shortInt), dimension(:), allocatable :: childrenIdxs, elementIdxs
    logical(defBool)                             :: isActive = .true., isBoundary = .false.
    real(defReal)                                :: area = ZERO
    real(defReal), dimension(3)                  :: centroid = ZERO, normal = ZERO
    type(axisAlignedBoundingBox)                 :: boundingBox
    character(:), allocatable                    :: type
  contains
    procedure          :: addChildIdx
    procedure          :: addEdge
    procedure          :: addElementIdx
    procedure          :: addVertex
    procedure, private :: build
    procedure          :: computeIntersection
    procedure          :: deactivate
    procedure          :: distanceSquared
    procedure          :: getArea
    procedure          :: getBoundingBox
    procedure          :: getCentroid
    procedure          :: getChildrenIdxs
    procedure          :: getEdges
    procedure          :: getElementIdxs
    procedure          :: getFaceIdx
    procedure          :: getHasElements
    procedure          :: getIsActive
    procedure          :: getIsBoundary
    procedure          :: getNormal
    procedure          :: getType
    procedure          :: getVertices
    procedure          :: init
    generic            :: intersects => intersects_BoundingBox
    procedure, private :: intersects_BoundingBox
    generic            :: intersectsBoundingBox => intersectsBoundingBox_BoundingBox
    procedure, private :: intersectsBoundingBox_BoundingBox
    procedure          :: isPointInside
    procedure          :: kill
    procedure          :: setArea
    procedure          :: setCentroid
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
  elemental subroutine addElementIdx(self, idx)
    class(face), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx
    
    call append(self % elementIdxs, idx)

  end subroutine addElementIdx
  
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
  subroutine build(self)
    class(face), intent(inout)                         :: self
    integer(shortInt)                                  :: i, nVertices
    real(defReal), dimension(3, 3)                     :: triangleCoordsArray
    real(defReal), dimension(3, size(self % vertices)) :: allCoords
    real(defReal), dimension(3)                        :: normal, sumAreasCentroid, sumNormals
    real(defReal)                                      :: normalNorm, sumAreas
    character(*), parameter                            :: here = 'build (face_class.f90)'

    ! First retrieve the coordinates of all the vertices in the face.
    nVertices = size(self % vertices)
    do i = 1, nVertices
      if (.not. associated(self % vertices(i) % ptr)) &
      call fatalError(here, 'Face with index '//numToChar(self % getIdx())//' contains a null vertex pointer.')
      allCoords(:, i) = self % vertices(i) % ptr % getCoordinates()

    end do

    ! Check if the face is a triangle. If so, perform a direct computation to avoid round-off errors.
    if (nVertices == 3) then
      normal = computeTriangleNormal(allCoords)
      normalNorm = norm2(normal)
      self % area = HALF * normalNorm
      self % centroid = THIRD * sum(allCoords, 2)
      self % normal = normal / normalNorm

    else
      ! Calculate the polygon's geometric centroid.
      triangleCoordsArray(:, 3) = sum(allCoords, 2) / nVertices
      sumAreas = ZERO
      sumAreasCentroid = ZERO
      sumNormals = ZERO
      do i = 1, nVertices
        triangleCoordsArray(:, 1) = allCoords(:, i)
        triangleCoordsArray(:, 2) = allCoords(:, merge(1, i + 1, i == nVertices))

        normal = computeTriangleNormal(triangleCoordsArray)
        sumNormals = sumNormals + normal

        normalNorm = norm2(normal)
        sumAreas = sumAreas + normalNorm
        sumAreasCentroid = sumAreasCentroid + normalNorm * sum(triangleCoordsArray, 2)

      end do
      self % area = HALF * sumAreas
      self % centroid = THIRD * sumAreasCentroid / sumAreas
      self % normal = sumNormals / norm2(sumNormals)

    end if
    ! Initialise face bounding box.
    call self % boundingBox % computeBounds(allCoords)

  contains
    !!
    !!
    !!
    pure function computeTriangleNormal(array) result(n)
      real(defReal), dimension(3, 3), intent(in) :: array
      real(defReal), dimension(3)                :: n

      n = crossProduct(array(:, 2) - array(:, 1), array(:, 3) - array(:, 1))

    end function computeTriangleNormal
    
  end subroutine build

  !! Subroutine 'computeIntersection'
  !!
  !! Basic description:
  !!   Checks whether a line segment intersects the face and if so, computes the distance from
  !!   the segment's origin to the point of intersection.
  !!
  !! Detailed description:
  !!
  !! Arguments:
  !!   startPos [in]          -> 3-D coordinates of the line segment's origin.
  !!   endPos [in]            -> 3-D coordinates of the line segment's end.
  !!   firstVertexCoords [in] -> 3-D coordinates of the face's first vertex.
  !!   isIntersecting [out]   -> .true. if the line segment intersects the triangle.
  !!   d [out]                -> Distance from the line segment's origin to the point of intersection.
  !!
  subroutine computeIntersection(self, coords, d)
    class(face), intent(in)     :: self
    type(coord), intent(in)     :: coords
    real(defReal), intent(out)  :: d
    real(defReal), dimension(3) :: diff, normal, r, rIntersection, u
    real(defReal)               :: denominator, s

    ! Initialise d = INF.
    d = INF
    
    ! Retrieve the face's normal vector and pre-compute the difference between the line segment's end
    ! and beginning positions.
    normal = self % normal
    r = coords % getPosition()
    u = coords % getDirection()
    diff = coords % getEndPosition() - r
    denominator = dot_product(normal, diff)

    ! If the denominator is ZERO, return early since the line segment is parallel to the face's plane.
    ! Else, compute the fraction of the line segment required to intersect the face's plane, s.
    if (areEqual(denominator, ZERO)) return
    s = dot_product(normal, self % centroid - r) / denominator
    
    ! If s is ZERO, the line segment's origin is on the face. In this case return early if the segment
    ! points in the same direction as the face's normal.
    if (areEqual(s, ZERO) .and. dot_product(normal, u) >= ZERO) return
    
    ! If s < ZERO or s > ONE, return early since the intersection is outside the line segment.
    if (s < ZERO .or. s > ONE) return

    ! Compute the coordinates of the intersection point.
    diff = s * diff
    rIntersection = r + diff

    ! Check if the intersection point coordinates are inside the face.
    if (self % isPointInside(rIntersection)) d = norm2(diff)

  end subroutine computeIntersection

  !!
  !!
  !!
  elemental subroutine deactivate(self)
    class(face), intent(inout) :: self

    self % isActive = .false.

  end subroutine deactivate

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
    diff = r - self % centroid
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

  !! Function 'getBoundingBox'
  !!
  !! Basic description:
  !!   Returns the bounding box of the face.
  !!
  !! Result:
  !!   boundingBox -> 6-D coordinates of the bounding box of the face.
  !!
  pure function getBoundingBox(self) result(boundingBox)
    class(face), intent(in)      :: self
    type(axisAlignedBoundingBox) :: boundingBox
    
    boundingBox = self % boundingBox

  end function getBoundingBox
  
  !! Function 'getCentroid'
  !!
  !! Basic description:
  !!   Returns the centroid of the face.
  !!
  !! Result:
  !!   centroid -> 3-D coordinates of the centroid of the face.
  !!
  pure function getCentroid(self) result(centroid)
    class(face), intent(in)     :: self
    real(defReal), dimension(3) :: centroid
    
    centroid = self % centroid

  end function getCentroid

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
  pure function getElementIdxs(self) result(elementIdxs)
    class(face), intent(in)                                :: self
    integer(shortInt), dimension(size(self % elementIdxs)) :: elementIdxs
    
    elementIdxs = self % elementIdxs

  end function getElementIdxs

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

  !! Function 'getHasElements'
  !!
  !! Basic description:
  !!   Returns .true. if elementIdxs is allocated.
  !!
  !! Result:
  !!   hasElements -> .true. if elementIdxs is allocated.
  !!
  elemental function getHasElements(self) result(hasElements)
    class(face), intent(in) :: self
    logical(defBool)        :: hasElements

    hasElements = allocated(self % elementIdxs)

  end function getHasElements

  !!
  !!
  !!
  elemental function getIsActive(self) result(isActive)
    class(face), intent(in) :: self
    logical(defBool)        :: isActive

    isActive = self % isActive

  end function getIsActive

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
  subroutine init(self, info)
    class(face), intent(inout)      :: self
    type(buildFaceInfo), intent(in) :: info
    integer(shortInt)               :: nVertices
    character(*), parameter         :: here = 'init (face_class.f90)'

    ! Catch invalid number of vertices.
    nVertices = size(info % vertices)
    if (nVertices < 3) then
      call fatalError(here, 'A face must have at least three vertices. Has: '//numToChar(nVertices)//'.')

    elseif(nVertices == 3) then
      self % type = 'Triangle'

    else
      self % type = 'Polygon'

    end if

    ! Set everything from payload.
    call self % setIdx(info % idx)
    self % parentIdx = info % parentIdx
    self % isBoundary = info % isBoundary
    self % vertices = info % vertices
    self % edges = info % edges

    ! Build.
    call self % build()

    ! Check if normal test was requested.
    if (info % testNormal) then
      if (dot_product(self % centroid - info % testCentroid, self % normal) < ZERO) then
        self % vertices(1) = info % vertices(2)
        self % vertices(2) = info % vertices(1)
        self % normal = -self % normal

      end if

    end if

  end subroutine init

  !!
  !!
  !!
  subroutine intersects_BoundingBox(self, boundingBox, doesIt)
    class(face), intent(in)                            :: self
    type(axisAlignedBoundingBox), intent(in)           :: boundingBox
    logical(defBool), intent(out)                      :: doesIt
    real(defReal), dimension(3)                        :: boundingBoxCentre, halfwidths, axis, edgeVector, boxAxis
    real(defReal), dimension(3, size(self % vertices)) :: centredVertexCoords
    integer(shortInt)                                  :: i, j, nVertices

    ! Initialise doesIt = .false., retrieve the centre and halfwidths of the boundingBox.
    doesIt = .false.
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

      ! Compute d and initialise minProjection and maxProjections.
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
  elemental subroutine intersectsBoundingBox_BoundingBox(self, boundingBox, doesIt)
    class(face), intent(in)                  :: self
    type(axisAlignedBoundingBox), intent(in) :: boundingBox
    logical(defBool), intent(out)            :: doesIt

    doesIt = self % boundingBox % intersects(boundingbox)

  end subroutine intersectsBoundingBox_BoundingBox

  !!
  !!
  !!
  function isPointInside(self, r) result(isIt)
    class(face), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    logical(defBool)                        :: isIt
    integer(shortInt)                       :: i, nextIdx, nVertices
    real(defReal)                           :: dotProduct
    real(defReal), dimension(3)             :: vertexCoords

    ! Initialise isIt = .false. and compute the number of vertices in the face.
    isIt = .false.
    nVertices = size(self % vertices)

    ! Loop through all the edges in the face and check if the point lies on the same side
    ! of each edge (note: this assumes a consistent vertex numbering).
    do i = 1, nVertices
      nextIdx = merge(1, i + 1, i == nVertices)
      vertexCoords = self % vertices(i) % ptr % getCoordinates()
      dotProduct = dot_product(self % normal, &
                               crossProduct(self % vertices(nextIdx) % ptr % getCoordinates() - vertexCoords, &
                                            r - vertexCoords))
      
      if (dotProduct < ZERO) return

    end do

    ! If reached here, the point is inside the face.
    isIt = .true.

  end function isPointInside
  
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
    self % isActive = .true.
    self % isBoundary = .false.
    self % area = ZERO
    self % centroid = ZERO
    self % normal = ZERO
    call self % boundingBox % kill()
    if (allocated(self % childrenIdxs)) deallocate(self % childrenIdxs)
    if (allocated(self % elementIdxs)) deallocate(self % elementIdxs)
    if (allocated(self % edges)) then
      do i = 1, size(self % edges)
        nullify(self % edges(i) % ptr)

      end do

    end if

    if (allocated(self % vertices)) then
      do i = 1, size(self % vertices)
        nullify(self % vertices(i) % ptr)

      end do

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

  !! Subroutine 'setCentroid'
  !!
  !! Basic description:
  !!   Sets the centroid of the face.
  !!
  !! Arguments:
  !!   centroid [in] -> 3-D coordinates of the centroid of the face.
  !!
  pure subroutine setCentroid(self, centroid)
    class(face), intent(inout)               :: self
    real(defReal), dimension(3), intent(in)  :: centroid

    self % centroid = centroid

  end subroutine setCentroid

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