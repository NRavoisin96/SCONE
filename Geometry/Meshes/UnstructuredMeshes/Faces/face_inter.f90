module face_inter
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use coord_class,                  only : coord
  use edgeShelf_class,              only : edgeShelf
  use genericProcedures,            only : append, areEqual, crossProduct, fatalError, findCommon, numToChar, swap
  use numPrecision
  use universalVariables,           only : HALF, INF, ONE, SURF_TOL, THIRD, ZERO
  use vertexShelf_class,            only : vertexShelf
  
  implicit none
  private

  ! Extendable procedures.
  public :: kill
  
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
  type, public, abstract                         :: face
    private
    integer(shortInt)                            :: idx = 0, parentIdx = 0
    integer(shortInt), dimension(:), allocatable :: edgeIdxs, elementIdxs, triangleIdxs, vertexIdxs, &
                                                    normalSigns
    logical(defBool)                             :: isBoundary = .false.
    real(defReal)                                :: area = ZERO, const = ZERO, extraDistance = ZERO
    real(defReal), dimension(3)                  :: centroid = ZERO, normal = ZERO, AB = ZERO, AC = ZERO
    type(axisAlignedBoundingBox)                 :: boundingBox
    character(:), allocatable                    :: type
    real(defReal), dimension(:), allocatable     :: extraDistanceArr
  contains
    procedure, non_overridable                   :: addEdgeIdx
    procedure, non_overridable                   :: addElementIdx
    procedure, non_overridable                   :: swapElementIdxsOrder
    procedure, non_overridable                   :: addTriangleIdx
    procedure, non_overridable                   :: addVertexIdx
    procedure, non_overridable                   :: build
    procedure(computeComponents), deferred       :: computeComponents
    procedure, non_overridable                   :: computeIntersection
    procedure(createTriangle), deferred          :: createTriangle
    procedure, non_overridable                   :: distanceSquared
    procedure, non_overridable                   :: distanceSquaredToEdge
    procedure, non_overridable                   :: getAB
    procedure, non_overridable                   :: getAC
    procedure, non_overridable                   :: getArea
    procedure, non_overridable                   :: getConst
    procedure, non_overridable                   :: getExtraDistance
    procedure, non_overridable                   :: getExtraDistanceArr
    procedure, non_overridable                   :: getNormalSigns
    procedure, non_overridable                   :: getBoundingBox
    procedure, non_overridable                   :: getCentroid
    procedure, non_overridable                   :: getEdgeIdxs
    procedure, non_overridable                   :: getElementIdxs
    procedure, non_overridable                   :: getFaceIdx
    procedure, non_overridable                   :: getHasElements
    procedure, non_overridable                   :: getIdx
    procedure, non_overridable                   :: getIsBoundary
    procedure, non_overridable                   :: getNormal
    procedure, non_overridable                   :: getTriangleIdxs
    procedure, non_overridable                   :: getType
    procedure, non_overridable                   :: getVertexIdxs
    procedure, non_overridable                   :: init
    generic                                      :: intersects => intersects_BoundingBox
    procedure, private, non_overridable          :: intersects_BoundingBox
    generic                                      :: intersectsBoundingBox => intersectsBoundingBox_BoundingBox
    procedure, private, non_overridable          :: intersectsBoundingBox_BoundingBox
    procedure, non_overridable                   :: isPointInside
    procedure                                    :: kill
    procedure, non_overridable                   :: setArea
    procedure, non_overridable                   :: setConst
    procedure, non_overridable                   :: setExtraDistance
    procedure, non_overridable                   :: setExtraDistanceArr
    procedure, non_overridable                   :: deallocateExtraDistanceArr
    procedure, non_overridable                   :: setNormalSigns
    procedure, non_overridable                   :: setCentroid
    procedure, non_overridable                   :: setIsBoundary
    procedure, non_overridable                   :: setIdx
    procedure, non_overridable                   :: setNormal
    procedure, non_overridable                   :: setVertexIdxs
    procedure                                    :: split
  end type face

  !!
  !! Small, local container to store polymorphic faces in a single array.
  !!
  !! Public members:
  !!   name -> Name of the mesh.
  !!   ptr  -> Pointer to the mesh.
  !!
  type, public               :: faceBox
    class(face), allocatable :: item
  end type

  abstract interface

    !! Subroutine 'computeAreaAndNormal'
    !!
    !! Basic description:
    !!   Computes the area and the normal vector of the face.
    !!
    !! Detailed description:
    !!   First estimates the centroid of the face by performing a simple arithmetic average of the 
    !!   vertices' coordinates. Using this centroid, the face is then decomposed into a number of 
    !!   triangles equal to the number of edges in the face where each triangle shares the face 
    !!   centroid as a common point. The normal vector of each triangle is then computed and the 
    !!   normal vector for the face is then obtained by summing the normal vectors for each triangle.
    !!   The overall face area is obtained by normalising the overall normal vector. The face centroid 
    !!   is then recomputed by performing an area-weighted average of the triangles' centres.
    !!
    !! Arguments:
    !!   vertices -> A vertexShelf.
    !!
    !! Errors:
    !!   fatalError if area is negative.
    !!   fatalError if area is infinite.
    !!
    pure subroutine computeComponents(self, vertexIdxs, vertices, centroid, normal, area)
      import                                      :: face, shortInt, vertexShelf, defReal
      class(face), intent(inout)                  :: self
      integer(shortInt), dimension(:), intent(in) :: vertexIdxs
      type(vertexShelf), intent(in)               :: vertices
      real(defReal), dimension(3), intent(out)    :: centroid, normal
      real(defReal), intent(out)                  :: area

    end subroutine computeComponents

    !!
    !!
    !!
    pure subroutine createTriangle(self, lastNewFaceIdx, edgeIdxs, newVertices, newTriangle, vertexIdxs, boundingBox)
      import                                         :: axisAlignedBoundingBox, face, shortInt, vertexShelf, faceBox, &
                                                        defReal
      class(face), intent(in)                        :: self
      integer(shortInt), intent(in)                  :: lastNewFaceIdx
      integer(shortInt), dimension(3), intent(in)    :: edgeIdxs
      type(vertexShelf), intent(in)                  :: newVertices
      type(faceBox), intent(inout)                   :: newTriangle
      integer(shortInt), dimension(3), intent(inout) :: vertexIdxs
      type(axisAlignedBoundingBox), intent(in)       :: boundingBox

    end subroutine createTriangle

  end interface

contains

  !! Subroutine 'addEdgeIdx'
  !!
  !! Basic description:
  !!   Adds the index of an edge sharing the face.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the edge sharing the face.
  !!
  elemental subroutine addEdgeIdx(self, idx)
    class(face), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx

    call append(self % edgeIdxs, idx)

  end subroutine addEdgeIdx
  
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

  !! Changes the order of element indices of a give face so that
  !! the order is [owner element, non-owner element]
  !! For (usual OpenFoam mesh format) boundary faces, the only element attached is
  !! usually set as owner element so [owner element]. However, if this format is violated,
  !! [0, the only element attached] where "0" indicates the outside of the mesh domain.
  subroutine swapElementIdxsOrder(self)
    class(face), intent(inout)                   :: self
    integer(shortInt), dimension(:), allocatable :: idxsOld

    if (allocated(idxsOld)) deallocate(idxsold)

    idxsOld = self % elementIdxs
    if (size(self % elementIdxs) == 1) then
      self % elementIdxs = 0_shortInt
      call append(self % elementIdxs, idxsOld)
    else
      self % elementIdxs(1) = idxsOld(2)
      self % elementIdxs(2) = idxsOld(1)
    end if

  end subroutine swapElementIdxsOrder


  !! Subroutine 'addTriangleIdx'
  !!
  !! Basic description:
  !!   Adds the index of a triangle in the face.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the triangle.
  !!
  elemental subroutine addTriangleIdx(self, idx)
    class(face), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx
    
    call append(self % triangleIdxs, idx)

  end subroutine addTriangleIdx
  
  !! Subroutine 'addVertexIdx'
  !!
  !! Basic description:
  !!   Adds the index of a vertex in the face.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the vertex.
  !!
  elemental subroutine addVertexIdx(self, idx)
    class(face), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx
    
    call append(self % vertexIdxs, idx)

  end subroutine addVertexIdx

  !!
  !!
  !!
  pure subroutine build(self, idx, faceIdx, isBoundary, vertexIdxs, vertices, type, boundingBox, testCentroid, edgeIdxs)
    class(face), intent(inout)                            :: self
    integer(shortInt), intent(in)                         :: idx, faceIdx
    logical(defBool), intent(in)                          :: isBoundary
    integer(shortInt), dimension(:), intent(inout)        :: vertexIdxs
    type(vertexShelf), intent(in)                         :: vertices
    character(*), intent(in)                              :: type
    type(axisAlignedBoundingBox), intent(in)              :: boundingBox
    real(defReal), dimension(3), intent(in), optional     :: testCentroid
    integer(shortInt), dimension(:), intent(in), optional :: edgeIdxs
    real(defReal), dimension(3)                           :: centroid, normal, firstVertexCoords, AB, AC
    real(defReal)                                         :: area

    ! Compute centroid, normal vector and area.
    call self % computeComponents(vertexIdxs, vertices, centroid, normal, area)

    ! Check normal orientation.
    if (present(testCentroid)) then
      ! Check that the triangle's normal vector points in the correct direction (outward). If not,
      ! swap two of the new triangle's vertices and negate its normal vector.
      if (dot_product(centroid - testCentroid, normal) <= ZERO) then
        call swap(vertexIdxs, 1, 2)
        normal = -normal

      end if

    end if

    ! Compute the two edge vectors of the face then set everything.
    firstVertexCoords = vertices % getVertexCoordinates(vertexIdxs(1))
    AB = vertices % getVertexCoordinates(vertexIdxs(2)) - firstVertexCoords
    AC = vertices % getVertexCoordinates(vertexIdxs(3)) - firstVertexCoords
    call self % init(idx, faceIdx, isBoundary, area, centroid, normal, AB, AC, vertexIdxs, type, boundingBox, edgeIdxs)

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
  pure subroutine computeIntersection(self, coords, vertices, d)
    class(face), intent(in)                            :: self
    type(coord), intent(in)                            :: coords
    type(vertexShelf), intent(in)                      :: vertices
    real(defReal), intent(out)                         :: d
    real(defReal), dimension(3)                        :: diff, normal, r, rIntersection, u
    real(defReal)                                      :: denominator, s

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
    if (self % isPointInside(rIntersection, vertices)) d = norm2(diff)

  end subroutine computeIntersection

  !!
  !!
  !!
  pure function distanceSquared(self, r, vertices) result(dSquared)
    class(face), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    type(vertexShelf), intent(in)           :: vertices
    real(defReal)                           :: d, dSquared, inverseNormalSquared, temp
    real(defReal), dimension(3)             :: diff, proj
    integer(shortInt)                       :: i, nextIdx, nVertices

    ! First compute the distance between the point and the plane of the face.
    diff = r - self % centroid
    d = dot_product(diff, self % normal)

    ! Now project the point on the plane of the face and check if the projection lies inside the face.
    inverseNormalSquared = ONE / dot_product(self % normal, self % normal)
    proj = r - self % normal * d * inverseNormalSquared

    ! If projection is inside the face, compute dSquared and return.
    if (self % isPointInside(proj, vertices)) then
      dSquared = d * d * inverseNormalSquared
      return

    end if

    ! If projection is outside the face, we need to compute the distance to each edge of the face and
    ! retain the mininum distance.
    dSquared = INF
    nVertices = size(self % vertexIdxs)
    do i = 1, nVertices
      nextIdx = merge(1, i + 1, i == nVertices)
      temp = self % distanceSquaredToEdge(r, vertices, i, nextIdx)
      dSquared = min(dSquared, temp)

    end do

  end function distanceSquared

  !!
  !!
  !!
  pure function distanceSquaredToEdge(self, r, vertices, idx, nextIdx) result(dSquared)
    class(face), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    type(vertexShelf), intent(in)           :: vertices
    integer(shortInt), intent(in)           :: idx, nextIdx
    real(defReal), dimension(3)             :: edgeVector, pointVector, vertexCoords
    real(defReal)                           :: dSquared, lSquared, t

    ! First compute edgeVector and pointVector.
    vertexCoords = vertices % getVertexCoordinates(idx)
    edgeVector = vertices % getVertexCoordinates(nextIdx) - vertexCoords
    pointVector = r - vertexCoords

    ! Compute the square of the edge length.
    lSquared = dot_product(edgeVector, edgeVector)

    ! Handle the case of a zero-length segment.
    if (areEqual(lSquared, ZERO)) then
        dSquared = dot_product(pointVector, pointVector)
        return
        
    end if

    ! Compute the normalisation parameter t by projecting pointVector onto edgeVector and
    ! snap it to the range [0, 1].
    t = max(ZERO, min(ONE, dot_product(pointVector, edgeVector) / lSquared))

    ! Now compute dSquared.
    pointVector = pointVector - edgeVector * t
    dSquared = dot_product(pointVector, pointVector)

  end function distanceSquaredToEdge

  !!
  !!
  !!
  pure function getAB(self) result(AB)
    class(face), intent(in)     :: self
    real(defReal), dimension(3) :: AB

    AB = self % AB

  end function getAB

  !!
  !!
  !!
  pure function getAC(self) result(AC)
    class(face), intent(in)     :: self
    real(defReal), dimension(3) :: AC

    AC = self % AC

  end function getAC
  
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
  elemental function getConst(self, idx) result(const)
    class(face), intent(in)                 :: self
    integer(shortInt), intent(in), optional :: idx
    real(defReal)                           :: const
    
    const = self % const

    if (.not. present(idx)) return
    if (idx < 0) const = -const

  end function getConst

  !!
  !!
  !!
  elemental function getExtraDistance(self) result(extraDistance)
    class(face), intent(in) :: self
    real(defReal)           :: extraDistance
    
    extraDistance = self % extraDistance

  end function getExtraDistance

  !!
  !!
  !!
  elemental function getExtraDistanceArr(self, currLayer) result(extraDistanceArr)
    class(face), intent(in)         :: self
    integer(shortInt), intent(in)   :: currLayer
    real(defReal)                   :: extraDistanceArr
    
    extraDistanceArr = self % extraDistanceArr(currLayer)

  end function getExtraDistanceArr

  !!
  !!
  !!
  pure function getNormalSigns(self, idx) result(normalSigns)
    class(face), intent(in)                 :: self
    integer(shortInt), intent(in), optional :: idx
    integer(shortInt), dimension(3)         :: normalSigns
    
    normalSigns = self % normalSigns

    if (.not. present(idx)) return
    if (idx < 0) normalSigns = -normalSigns

  end function getNormalSigns

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

  !! Function 'getEdgeIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the edges in the face.
  !!
  !! Result:
  !!   edgeIdxs -> Indices of the edges in the face.
  !!
  pure function getEdgeIdxs(self) result(edgeIdxs)
    class(face), intent(in)                             :: self
    integer(shortInt), dimension(size(self % edgeIdxs)) :: edgeIdxs

    edgeIdxs = self % edgeIdxs

  end function getEdgeIdxs
  
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
  
  !! Function 'getIdx'
  !!
  !! Basic description:
  !!   Returns the index of the face.
  !!
  !! Result:
  !!   idx -> Index of the face.
  !!
  elemental function getIdx(self) result(idx)
    class(face), intent(in) :: self
    integer(shortInt)       :: idx
    
    idx = self % idx

  end function getIdx

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
  !!   Returns the indices of the triangles in the face.
  !!
  !! Result:
  !!   trianglesIdxs -> Indices of the triangles in the face.
  !!
  pure function getTriangleIdxs(self) result(triangleIdxs)
    class(face), intent(in)                                 :: self
    integer(shortInt), dimension(size(self % triangleIdxs)) :: triangleIdxs
    
    triangleIdxs = self % triangleIdxs

  end function getTriangleIdxs

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
  pure function getVertexIdxs(self) result(vertexIdxs)
    class(face), intent(in)                               :: self
    integer(shortInt), dimension(size(self % vertexIdxs)) :: vertexIdxs
    
    vertexIdxs = self % vertexIdxs

  end function getVertexIdxs

  !!
  !!
  !!
  pure subroutine init(self, idx, parentIdx, isBoundary, area, centroid, normal, AB, AC, vertexIdxs, type, boundingBox, edgeIdxs)
    class(face), intent(inout)                            :: self
    integer(shortInt), intent(in)                         :: idx, parentIdx
    logical(defBool), intent(in)                          :: isBoundary
    real(defReal), intent(in)                             :: area
    real(defReal), dimension(3), intent(in)               :: centroid, normal, AB, AC
    integer(shortInt), dimension(:), intent(in)           :: vertexIdxs
    character(*), intent(in)                              :: type
    type(axisAlignedBoundingBox), intent(in)              :: boundingBox
    integer(shortInt), dimension(:), intent(in), optional :: edgeIdxs

    ! Set everything.
    self % idx = idx
    self % parentIdx = parentIdx
    self % isBoundary = isBoundary
    self % area = area
    self % centroid = centroid
    self % normal = normal
    self % AB = AB
    self % AC = AC
    self % vertexIdxs = vertexIdxs
    self % type = type
    self % boundingBox = boundingBox
    if (present(edgeIdxs)) self % edgeIdxs = edgeIdxs

  end subroutine init

  !!
  !!
  !!
  elemental subroutine intersects_BoundingBox(self, vertices, boundingBox, doesIt)
    class(face), intent(in)                              :: self
    type(vertexShelf), intent(in)                        :: vertices
    type(axisAlignedBoundingBox), intent(in)             :: boundingBox
    logical(defBool), intent(out)                        :: doesIt
    real(defReal), dimension(3)                          :: boundingBoxCentre, halfwidths, axis, edge, boxAxis
    real(defReal), dimension(3, size(self % vertexIdxs)) :: centredVertexCoords
    integer(shortInt)                                    :: i, j, nextIdx, nVertices

    ! Initialise doesIt = .false., retrieve the centre and halfwidths of the boundingBox.
    doesIt = .false.
    boundingBoxCentre = boundingBox % getCentre()
    halfwidths = boundingBox % getHalfwidths()

    ! Offset the coordinates of the face vertices with respect to the box centre.
    nVertices = size(self % vertexIdxs)
    centredVertexCoords = vertices % getVertexCoordinates(self % vertexIdxs) - spread(boundingBoxCentre, 2, nVertices)

    ! First test for intersection along the three bounding box's axes.
    do i = 1, 3
      axis = ZERO
      axis(i) = ONE
      if (.not. overlaps(halfwidths, centredVertexCoords, axis, nVertices)) return

    end do

    ! Now test the face's normal vector.
    if (.not. overlaps(halfwidths, centredVertexCoords, self % normal, nVertices)) return

    ! Finally, test cross products between the face's edges and the bounding box's edges.
    do i = 1, nVertices
      nextIdx = merge(1, i + 1, i == nVertices)
      edge = centredVertexCoords(:, nextIdx) - centredVertexCoords(:, i)
      do j = 1, 3
        boxAxis = ZERO
        boxAxis(j) = ONE
        axis = crossProduct(edge, boxAxis)
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
      real(defReal), dimension(3), intent(in)                          :: h, ax
      real(defReal), dimension(3, size(self % vertexIdxs)), intent(in) :: coords
      integer(shortInt), intent(in)                                    :: n
      logical(defBool)                                                 :: isOverlapping
      real(defReal)                                                    :: radius, minProjection, maxProjection, d
      integer(shortInt)                                                :: k

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
      isOverlapping = minProjection <= radius .and. maxProjection >= -radius

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
  pure function isPointInside(self, r, vertices) result(isIt)
    class(face), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    type(vertexShelf), intent(in)           :: vertices
    logical(defBool)                        :: isIt
    integer(shortInt)                       :: i, nextIdx, nVertices
    real(defReal)                           :: dotProduct
    real(defReal), dimension(3)             :: vertexCoords

    ! Initialise isIt = .false. and compute the number of vertices in the face.
    isIt = .false.
    nVertices = size(self % vertexIdxs)

    ! Loop through all the edges in the face and check if the point lies on the same side
    ! of each edge (note: this assumes a consistent vertex numbering).
    do i = 1, nVertices
      nextIdx = merge(1, i + 1, i == nVertices)
      vertexCoords = vertices % getVertexCoordinates(self % vertexIdxs(i))
      dotProduct = dot_product(self % normal, &
                               crossProduct(vertices % getVertexCoordinates(self % vertexIdxs(nextIdx)) - vertexCoords, &
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
    
    self % idx = 0
    self % parentIdx = 0
    self % isBoundary = .false.
    self % area = ZERO
    self % centroid = ZERO
    self % normal = ZERO
    self % AB = ZERO
    self % AC = ZERO
    call self % boundingBox % kill()
    if (allocated(self % edgeIdxs)) deallocate(self % edgeIdxs)
    if (allocated(self % elementIdxs)) deallocate(self % elementIdxs)
    if (allocated(self % vertexIdxs)) deallocate(self % vertexIdxs)

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
  elemental subroutine setExtraDistance(self, extraDistance)
    class(face), intent(inout) :: self
    real(defReal), intent(in)  :: extraDistance

    self % extraDistance = extraDistance

  end subroutine setExtraDistance

  !!
  !!
  !!
  pure subroutine setExtraDistanceArr(self, extraDistanceArr, n_layers)
    class(face), intent(inout)              :: self
    real(defReal), dimension(:), intent(in) :: extraDistanceArr
    integer(shortInt), intent(in)           :: n_layers

    allocate(self % extraDistanceArr(n_layers))
    self % extraDistanceArr = extraDistanceArr

  end subroutine setExtraDistanceArr

  !!
  !!
  !!
  elemental subroutine deallocateExtraDistanceArr(self)
    class(face), intent(inout)      :: self

    deallocate(self % extraDistanceArr)

  end subroutine deallocateExtraDistanceArr

  !!
  !!
  !!
  pure subroutine setNormalSigns(self, normalSigns)
    class(face), intent(inout)                  :: self
    integer(shortInt), dimension(3), intent(in) :: normalSigns

    if (.NOT. allocated(self % normalSigns)) allocate(self % normalSigns(3))
    self % normalSigns = normalSigns

  end subroutine setNormalSigns
  !!
  !!
  !!
  elemental subroutine setConst(self, const)
    class(face), intent(inout) :: self
    real(defReal), intent(in)  :: const

    self % const = const

  end subroutine setConst
  
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
  
  !! Subroutine 'setIdx'
  !!
  !! Basic description:
  !!   Sets the index of the face.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the face.
  !!
  !! Error:
  !!   fatalError if idx < 1.
  !!
  subroutine setIdx(self, idx)
    class(face), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx
    character(100), parameter     :: Here = 'setIdx (face_class.f90)'
    
    ! Catch invalid index.
    if (idx < 1) call fatalError(Here, 'Face index must be +ve. Is: '//numToChar(idx)//'.')
    self % idx = idx

  end subroutine setIdx

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
  pure subroutine setVertexIdxs(self, vertexIdxs)
    class(face), intent(inout)                  :: self
    integer(shortInt), dimension(:), intent(in) :: vertexIdxs

    self % vertexIdxs = vertexIdxs

  end subroutine setVertexIdxs

  !!
  !!
  !!
  subroutine split(self, newEdges, newVertices, lastNewEdgeIdx, lastNewFaceIdx, triangles)
    class(face), intent(inout)                 :: self
    type(edgeShelf), intent(inout)             :: newEdges
    type(vertexShelf), intent(inout)           :: newVertices
    integer(shortInt), intent(inout)           :: lastNewEdgeIdx, lastNewFaceIdx
    type(faceBox), dimension(:), intent(inout) :: triangles
    integer(shortInt)                          :: i, j, k, minVertexLoc, nTriangles, nVertices
    integer(shortInt), dimension(3)            :: edgeIdxs, vertexIdxs
    real(defReal), dimension(3, 3)             :: vertexCoords
    type(axisAlignedBoundingBox)               :: boundingBox

    ! First compute the number of vertices and the location of the vertex of minimum index in the current face.
    nVertices = size(self % vertexIdxs)
    minVertexLoc = minloc(self % vertexIdxs, 1)
    vertexIdxs(1) = self % vertexIdxs(minVertexLoc)

    ! Compute the number of triangles to be generated, allocate memory and loop through all new triangles.
    nTriangles = nVertices - 2
    allocate(self % triangleIdxs(nTriangles))
    do i = 1, nTriangles
      ! Increment lastNewFaceIdx and add the new triangle to the list of face's triangle indices.
      lastNewFaceIdx = lastNewFaceIdx + 1
      self % triangleIdxs(i) = lastNewFaceIdx

      ! Compute the locations of the second and third vertices in the new triangle then set their indices.
      if (nVertices == 3) then
        edgeIdxs = self % edgeIdxs
        vertexIdxs = self % vertexIdxs
        boundingBox = self % getBoundingBox()

      else
        vertexIdxs(2:3) = self % vertexIdxs([mod(minVertexLoc + i - 1, nVertices) + 1, mod(minVertexLoc + i, nVertices) + 1])

        ! If we are not at the last iteration of the loop, create a new edge.
        if (i < nTriangles) then
          lastNewEdgeIdx = lastNewEdgeIdx + 1
          call newEdges % initEdge(lastNewEdgeIdx, vertexIdxs([1, 3]))
          call newVertices % addEdgeIdxToVertex(vertexIdxs(1), lastNewEdgeIdx)
          call newVertices % addEdgeIdxToVertex(vertexIdxs(3), lastNewEdgeIdx)

        end if

        ! Loop through all vertices (and edges) in the new triangle.
        do j = 1, 3
          edgeIdxs(j) = newVertices % findCommonEdgeIdx(vertexIdxs(j), vertexIdxs(mod(j, 3) + 1))
          vertexCoords(:, j) = newVertices % getVertexCoordinates(vertexIdxs(j))

        end do
        call boundingBox % computeBounds(vertexCoords)

      end if

      ! Create a new triangle.
      call self % createTriangle(lastNewFaceIdx, edgeIdxs, newVertices, triangles(lastNewFaceIdx), vertexIdxs, boundingBox)

      ! Update mesh connectivity information.
      do j = 1, 3
        call newEdges % addFaceIdxToEdge(edgeIdxs(j), lastNewFaceIdx)
        call newVertices % addFaceIdxToVertex(vertexIdxs(j), lastNewFaceIdx)

      end do

    end do

  end subroutine split
  
end module face_inter