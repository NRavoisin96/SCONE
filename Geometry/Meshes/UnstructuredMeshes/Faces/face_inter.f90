module face_inter
  
  use edgeShelf_class,     only : edgeShelf
  use genericProcedures,   only : append, areEqual, fatalError, findCommon, numToChar, swap
  use numPrecision
  use universalVariables,  only : INF, HALF, THIRD, SURF_TOL, ZERO
  use vertexShelf_class,   only : vertexShelf
  
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
    integer(shortInt), dimension(:), allocatable :: edgeIdxs, elementIdxs, triangleIdxs, vertexIdxs
    logical(defBool)                             :: isBoundary = .false.
    real(defReal)                                :: area = ZERO
    real(defReal), dimension(3)                  :: centroid = ZERO, normal = ZERO, AB = ZERO, AC = ZERO
    character(:), allocatable                    :: type
  contains
    procedure, non_overridable                   :: addEdgeIdx
    procedure, non_overridable                   :: addElementIdx
    procedure, non_overridable                   :: addTriangleIdx
    procedure, non_overridable                   :: addVertexIdx
    procedure, non_overridable                   :: build
    procedure(computeComponents), deferred       :: computeComponents
    procedure, non_overridable                   :: computeIntersection
    procedure(createTriangle), deferred          :: createTriangle
    procedure, non_overridable                   :: getAB
    procedure, non_overridable                   :: getAC
    procedure, non_overridable                   :: getArea
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
    procedure                                    :: kill
    procedure, non_overridable                   :: setArea
    procedure, non_overridable                   :: setCentroid
    procedure, non_overridable                   :: setIsBoundary
    procedure, non_overridable                   :: setIdx
    procedure, non_overridable                   :: setNormal
    procedure, non_overridable                   :: setVertexIdxs
    procedure                                    :: split
    procedure(testForInclusion), deferred        :: testForInclusion
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
    pure subroutine createTriangle(self, lastNewFaceIdx, edgeIdxs, newVertices, newTriangle, vertexIdxs)
      import                                         :: face, shortInt, vertexShelf, faceBox
      class(face), intent(in)                        :: self
      integer(shortInt), intent(in)                  :: lastNewFaceIdx
      integer(shortInt), dimension(3), intent(in)    :: edgeIdxs
      type(vertexShelf), intent(in)                  :: newVertices
      type(faceBox), intent(inout)                   :: newTriangle
      integer(shortInt), dimension(3), intent(inout) :: vertexIdxs

    end subroutine createTriangle

    !!
    !!
    !!
    pure subroutine testForInclusion(self, vertices, intersectionCoords, diff, d, edgeIdx, vertexIdx)
      import                                  :: face, defReal, shortInt, vertexShelf
      class(face), intent(in)                 :: self
      type(vertexShelf), intent(in)           :: vertices
      real(defReal), dimension(3), intent(in) :: intersectionCoords, diff
      real(defReal), intent(inout)            :: d
      integer(shortInt), intent(inout)        :: edgeIdx, vertexIdx

    end subroutine testForInclusion

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
  pure subroutine build(self, idx, faceIdx, isBoundary, vertexIdxs, vertices, type, testCentroid, edgeIdxs)
    class(face), intent(inout)                            :: self
    integer(shortInt), intent(in)                         :: idx, faceIdx
    logical(defBool), intent(in)                          :: isBoundary
    integer(shortInt), dimension(:), intent(inout)        :: vertexIdxs
    type(vertexShelf), intent(in)                         :: vertices
    character(*), intent(in)                              :: type
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
    call self % init(idx, faceIdx, isBoundary, area, centroid, normal, AB, AC, vertexIdxs, type, edgeIdxs)

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
  pure subroutine computeIntersection(self, r, rEnd, u, vertices, d, edgeIdx, vertexIdx)
    class(face), intent(in)                            :: self
    real(defReal), dimension(3), intent(in)            :: r, rEnd, u
    type(vertexShelf), intent(in)                      :: vertices
    real(defReal), intent(out)                         :: d
    integer(shortInt), intent(out)                     :: edgeIdx, vertexIdx
    real(defReal), dimension(3)                        :: normal, diff, intersectionCoords
    real(defReal)                                      :: denominator, s

    ! Initialise d = INF, edgeIdx = 0 and vertexIdx = 0.
    d = INF
    edgeIdx = 0
    vertexIdx = 0
    
    ! Retrieve the face's normal vector and pre-compute the difference between the line segment's end
    ! and beginning positions.
    normal = self % normal
    diff = rEnd - r
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
    intersectionCoords = r + diff

    ! Check if the intersection point coordinates are inside the face.
    call self % testForInclusion(vertices, intersectionCoords, diff, d, edgeIdx, vertexIdx)

  end subroutine computeIntersection

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
  pure subroutine init(self, idx, parentIdx, isBoundary, area, centroid, normal, AB, AC, vertexIdxs, type, edgeIdxs)
    class(face), intent(inout)                            :: self
    integer(shortInt), intent(in)                         :: idx, parentIdx
    logical(defBool), intent(in)                          :: isBoundary
    real(defReal), intent(in)                             :: area
    real(defReal), dimension(3), intent(in)               :: centroid, normal, AB, AC
    integer(shortInt), dimension(:), intent(in)           :: vertexIdxs
    character(*), intent(in)                              :: type
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
    if (present(edgeIdxs)) self % edgeIdxs = edgeIdxs

  end subroutine init
  
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
    integer(shortInt)                          :: i, j, minVertexLoc, nTriangles, nVertices
    integer(shortInt), dimension(3)            :: edgeIdxs, vertexIdxs

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

      else
        vertexIdxs(2:3) = self % vertexIdxs([mod(minVertexLoc + i - 1, nVertices) + 1, mod(minVertexLoc + i, nVertices) + 1])

        ! If we are not at the last iteration of the loop, create a new edge.
        if (i < nTriangles) then
          lastNewEdgeIdx = lastNewEdgeIdx + 1
          call newEdges % initEdge(lastNewEdgeIdx, vertexIdxs([1, 3]))
          call newVertices % addEdgeIdxToVertex(vertexIdxs(1), lastNewEdgeIdx)
          call newVertices % addEdgeIdxToVertex(vertexIdxs(3), lastNewEdgeIdx)

        end if

        do j = 1, 3
          edgeIdxs(j) = newVertices % findCommonEdgeIdx(vertexIdxs(j), vertexIdxs(mod(j, 3) + 1))

        end do

      end if

      ! Create a new triangle.
      call self % createTriangle(lastNewFaceIdx, edgeIdxs, newVertices, triangles(lastNewFaceIdx), vertexIdxs)

      ! Update mesh connectivity information.
      do j = 1, 3
        call newEdges % addFaceIdxToEdge(edgeIdxs(j), lastNewFaceIdx)
        call newVertices % addFaceIdxToVertex(vertexIdxs(j), lastNewFaceIdx)

      end do

    end do

  end subroutine split
  
end module face_inter