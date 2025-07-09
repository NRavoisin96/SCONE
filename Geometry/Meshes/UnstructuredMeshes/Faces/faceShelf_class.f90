module faceShelf_class
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use coord_class,                  only : coord
  use edgeShelf_class,              only : edgeShelf
  use numPrecision
  use genericProcedures,            only : fatalError, findCommon, removeDuplicates, quickSort
  use face_inter,                   only : face, faceBox
  use polygon_class,                only : polygon
  use triangle_class,               only : triangle
  use vertexShelf_class,            only : vertexShelf
  
  implicit none
  private
  
  !!
  !! Storage space for faces of an OpenFOAM mesh.
  !!
  !! Private members:
  !!   shelf -> Array to store faces.
  !!
  type, public                               :: faceShelf
    private
    type(faceBox), dimension(:), allocatable :: shelf
  contains
    procedure                                :: addEdgeIdxToFace
    procedure                                :: addElementIdxToFace
    procedure                                :: addFace
    procedure                                :: addVertexIdxToFace
    procedure                                :: allocateFace
    procedure                                :: allocateShelf
    procedure                                :: buildFace
    procedure                                :: computeFaceIntersection
    procedure                                :: distanceSquaredFromFace
    procedure                                :: findCommonEdgeIdx
    procedure                                :: findCommonVertexIdx
    procedure                                :: getAllFaceBoundingBoxes
    procedure                                :: getAllFaceCentroids
    procedure                                :: getFaceArea
    procedure                                :: getFaceConst
    procedure                                :: setFaceConst
    procedure                                :: getFaceBoundingBox
    procedure                                :: getFaceCentroid
    procedure                                :: getFaceEdgeIdxs
    generic                                  :: getFaceElementIdxs => getFaceElementIdxs_shortInt, &
                                                                      getFaceElementIdxs_shortIntArray
    procedure, private                       :: getFaceElementIdxs_shortInt
    procedure, private                       :: getFaceElementIdxs_shortIntArray
    procedure                                :: getFaceHasElements
    procedure                                :: getFaceIsBoundary
    procedure                                :: getFaceNormal
    procedure                                :: getFaceTriangleIdxs
    procedure                                :: getFaceType
    procedure                                :: getFaceVertexIdxs
    procedure                                :: getSize
    procedure                                :: initFace
    generic                                  :: intersectsFace => intersectsFace_BoundingBox
    procedure, private                       :: intersectsFace_BoundingBox
    generic                                  :: intersectsFaceBoundingBox => intersectsFaceBoundingBox_BoundingBox
    procedure, private                       :: intersectsFaceBoundingBox_BoundingBox
    procedure                                :: kill
    procedure                                :: splitFace
  end type

contains

  !! Subroutine 'addEdgeIdxToFace'
  !!
  !! Basic description:
  !!   Adds the index of an edge to a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the face in the shelf.
  !!   edgeIdx [in] -> Index of the edge in the face.
  !!
  elemental subroutine addEdgeIdxToFace(self, idx, edgeIdx)
    class(faceShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx, edgeIdx

    call self % shelf(idx) % item % addEdgeIdx(edgeIdx)

  end subroutine addEdgeIdxToFace

  !! Subroutine 'addElementIdxToFace'
  !!
  !! Basic description:
  !!   Adds the index of an element to a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in]        -> Index of the face in the shelf.
  !!   elementIdx [in] -> Index of the element containing the face.
  !!
  elemental subroutine addElementIdxToFace(self, idx, elementIdx)
    class(faceShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx, elementIdx

    call self % shelf(idx) % item % addElementIdx(elementIdx)

  end subroutine addElementIdxToFace

  !!
  !!
  !!
  elemental subroutine addFace(self, idx, newFace)
    class(faceShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx
    type(faceBox), intent(in)       :: newFace

    self % shelf(idx) = newFace

  end subroutine addFace

  !! Subroutine 'addVertexIdxToFace'
  !!
  !! Basic description:
  !!   Adds the index of a vertex to a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in]       -> Index of the face in the shelf.
  !!   vertexIdx [in] -> Index of the vertex in the face.
  !!
  elemental subroutine addVertexIdxToFace(self, idx, vertexIdx)
    class(faceShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx, vertexIdx

    call self % shelf(idx) % item % addVertexIdx(vertexIdx)

  end subroutine addVertexIdxToFace

  !!
  !!
  !!
  elemental subroutine allocateFace(self, idx, type)
    class(faceShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx
    character(*), intent(in)        :: type

    select case(type)
      case('Polygon')
        allocate(polygon :: self % shelf(idx) % item)
      case('Triangle')
        allocate(triangle :: self % shelf(idx) % item)
      case default

    end select

  end subroutine allocateFace

  !! Subroutine 'allocateShelf'
  !!
  !! Basic description:
  !!   Allocates memory in the shelf.
  !!
  !! Arguments:
  !!   nFaces [in] -> Number of faces to be included in the shelf.
  !!
  elemental subroutine allocateShelf(self, nFaces)
    class(faceShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: nFaces

    allocate(self % shelf(nFaces))

  end subroutine allocateShelf

  !!
  !!
  !!
  pure subroutine buildFace(self, idx, faceIdx, isBoundary, vertexIdxs, vertices, type, boundingBox, testCentroid)
    class(faceShelf), intent(inout)                   :: self
    integer(shortInt), intent(in)                     :: idx, faceIdx
    logical(defBool), intent(in)                      :: isBoundary
    integer(shortInt), dimension(:), intent(inout)    :: vertexIdxs
    type(vertexShelf), intent(in)                     :: vertices
    character(*), intent(in)                          :: type
    type(axisAlignedBoundingBox), intent(in)          :: boundingBox
    real(defReal), dimension(3), intent(in), optional :: testCentroid

    ! Allocate the face in the shelf.
    call self % allocateFace(idx, type)

    ! Set face properties and compute its centroid, normal vector and area.
    call self % shelf(idx) % item % build(idx, faceIdx, isBoundary, vertexIdxs, vertices, type, boundingBox, testCentroid)

  end subroutine buildFace

  !! Subroutine 'computeFaceIntersection'
  !!
  !! Basic description:
  !!   Computes the intersection of a line segment of origin r, end rEnd and direction u with a
  !!   face in the shelf.
  !!
  !! Arguments:
  !!   idx [in]        -> Index of the face in the shelf.
  !!   r [in]          -> Line segment's origin coordinates.
  !!   rEnd [in]       -> Line segment's end coordinates.
  !!   u [in]          -> Line segment's direction
  !!   vertices [in]   -> A vertexShelf.
  !!   d [out]         -> Distance to intersection with the face.
  !!   edgeIdx [out]   -> Used in case the line segment intersects the face at one of its edges.
  !!   vertexIdx [out] -> Used in case the line segment intersects the face at one of its vertices.
  !!
  pure subroutine computeFaceIntersection(self, idx, coords, vertices, d)
    class(faceShelf), intent(in)            :: self
    integer(shortInt), intent(in)           :: idx
    type(coord), intent(in)                 :: coords
    type(vertexShelf), intent(in)           :: vertices
    real(defReal), intent(out)              :: d

    call self % shelf(idx) % item % computeIntersection(coords, vertices, d)

  end subroutine computeFaceIntersection

  !! Function 'distanceSquaredFromFace'
  !!
  !!
  pure function distanceSquaredFromFace(self, idx, r, vertices) result(dSquared)
    class(faceShelf), intent(in)            :: self
    integer(shortInt), intent(in)           :: idx
    real(defReal), dimension(3), intent(in) :: r
    type(vertexShelf), intent(in)           :: vertices
    real(defReal)                           :: dSquared

    dSquared = self % shelf(idx) % item % distanceSquared(r, vertices)

  end function distanceSquaredFromFace

  !! Function 'findCommonEdgeIdx'
  !!
  !! Basic description:
  !!   Returns the index of the common edge between two faces.
  !!
  !! Arguments:
  !!   firstFaceIdx [in]  -> Index of the first face.
  !!   secondFaceIdx [in] -> Index of the second face.
  !!
  !! Result:
  !!   edgeIdx            -> Index of the common edge between the two faces.
  !!
  elemental function findCommonEdgeIdx(self, firstFaceIdx, secondFaceIdx) result(edgeIdx)
    class(faceShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: firstFaceIdx, secondFaceIdx
    integer(shortInt)                            :: edgeIdx
    integer(shortInt), dimension(:), allocatable :: commonIdxs

    ! Initialise edgeIdx = 0 then find common edge indices between the two faces.
    edgeIdx = 0
    commonIdxs = findCommon(self % shelf(firstFaceIdx) % item % getEdgeIdxs(), self % shelf(secondFaceIdx) % item % getEdgeIdxs())

    ! If a common edge has been found update edgeIdx.
    if (size(commonIdxs) > 0) edgeIdx = commonIdxs(1)

  end function findCommonEdgeIdx

  !! Function 'findCommonVertexIdx'
  !!
  !! Basic description:
  !!   Returns the index of the common vertex between a set of faces in the shelf.
  !!
  !! Arguments:
  !!   idxs [in] -> Indices of the faces in the shelf.
  !!
  !! Result:
  !!   vertexIdx -> Index of the common vertex between the faces.
  !!
  pure function findCommonVertexIdx(self, idxs) result(vertexIdx)
    class(faceShelf), intent(in)                 :: self
    integer(shortInt), dimension(:), intent(in)  :: idxs
    integer(shortInt)                            :: vertexIdx, i
    integer(shortInt), dimension(:), allocatable :: commonIdxs

    ! Initialise vertexIdx = 0 then find common vertex indices between the faces.
    vertexIdx = 0
    commonIdxs = self % shelf(idxs(1)) % item % getVertexIdxs()

    do i = 2, size(idxs)
      commonIdxs = findCommon(commonIdxs, self % shelf(idxs(i)) % item % getVertexIdxs())

    end do

    ! If a common vertex has been found update vertexIdx.
    if (size(commonIdxs) > 0) vertexIdx = commonIdxs(1)

  end function findCommonVertexIdx

  !! Function 'getAllFaceBoundingBoxes'
  !!
  !! Basic description:
  !!   Returns the bounding boxes of all the faces in the shelf.
  !!
  !! Results:
  !!   boundingBoxes -> A defReal array containing the bounding boxes of all the faces in the shelf.
  !!
  pure function getAllFaceBoundingBoxes(self) result(boundingBoxes)
    class(faceShelf), intent(in)                              :: self
    type(axisAlignedBoundingBox), dimension(self % getSize()) :: boundingBoxes
    integer(shortInt)                                         :: i

    do i = 1, self % getSize()
      boundingBoxes(i) = self % shelf(i) % item % getBoundingBox()

    end do

  end function getAllFaceBoundingBoxes

  !! Function 'getAllFaceCentroids'
  !!
  !! Basic description:
  !!   Returns the centroids of all the faces in the shelf.
  !!
  !! Results:
  !!   centroids -> A defReal array containing the centroids of all the faces in the shelf.
  !!
  pure function getAllFaceCentroids(self) result(centroids)
    class(faceShelf), intent(in)                  :: self
    real(defReal), dimension(3, self % getSize()) :: centroids
    integer(shortInt)                             :: i

    do i = 1, self % getSize()
      centroids(:, i) = self % getFaceCentroid(i)

    end do

  end function getAllFaceCentroids

  !! Function 'getFaceArea'
  !!
  !! Basic description:
  !!   Returns the area of a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the face in the shelf.
  !!
  !! Result:
  !!   area     -> Area of the face's centroid.
  !!
  elemental function getFaceArea(self, idx) result(area)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    real(defReal)                 :: area

    area = self % shelf(idx) % item % getArea()

  end function getFaceArea

  !!
  !!
  !!
  elemental function getFaceConst(self, idx) result(const)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    real(defReal)                 :: const

    const = self % shelf(abs(idx)) % item % getConst(idx)

  end function getFaceConst

  !!
  !!
  !!
  elemental subroutine setFaceConst(self, idx, const)
    class(faceShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx
    real(defReal), intent(in)       :: const

    call self % shelf(idx) % item % setConst(const)

  end subroutine setFaceConst

  !! Function 'getFaceBoundingBox'
  !!
  !! Basic description:
  !!   Returns the bounding box of a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in]    -> Index of the face in the shelf.
  !!
  !! Result:
  !!   boundingBox -> 6-D coordinates of the face's bounding box.
  !!
  pure function getFaceBoundingBox(self, idx) result(boundingBox)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    type(axisAlignedBoundingBox)  :: boundingBox

    boundingBox = self % shelf(idx) % item % getBoundingBox()

  end function getFaceBoundingBox

  !! Function 'getFaceCentroid'
  !!
  !! Basic description:
  !!   Returns the 3-D coordinates of the centroid of a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the face in the shelf.
  !!
  !! Result:
  !!   centroid -> 3-D coordinates of the face's centroid.
  !!
  pure function getFaceCentroid(self, idx) result(centroid)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    real(defReal), dimension(3)   :: centroid

    centroid = self % shelf(idx) % item % getCentroid()

  end function getFaceCentroid

  !! Function 'getFaceEdgeIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the edges in a face of the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the face in the shelf.
  !!
  !! Result:
  !!   edgeIdxs -> Indices of the edges in the face.
  !!
  pure function getFaceEdgeIdxs(self, idx) result(edgeIdxs)
    class(faceShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: edgeIdxs

    edgeIdxs = self % shelf(idx) % item % getEdgeIdxs()

  end function getFaceEdgeIdxs

  !! Function 'getFaceElementIdxs_shortInt'
  !!
  !! Basic description:
  !!   Returns the indices of the elements sharing a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in]    -> Index of the face in the shelf.
  !!
  !! Result:
  !!   elementIdxs -> Indices of the elements sharing the face.
  !!
  pure function getFaceElementIdxs_shortInt(self, idx) result(elementIdxs)
    class(faceShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: elementIdxs

    elementIdxs = self % shelf(idx) % item % getElementIdxs()

  end function getFaceElementIdxs_shortInt

  !! Function 'getFaceElementIdxs_shortIntArray'
  !!
  !! Basic description:
  !!   Returns the unique indices of the elements sharing faces in the shelf.
  !!
  !! Arguments:
  !!   idxs [in]   -> Indices of the faces in the shelf.
  !!
  !! Result:
  !!   elementIdxs -> Indices of the elements sharing the faces.
  !!
  pure function getFaceElementIdxs_shortIntArray(self, idxs) result(elementIdxs)
    class(faceShelf), intent(in)                 :: self
    integer(shortInt), dimension(:), intent(in)  :: idxs
    integer(shortInt), dimension(:), allocatable :: elementIdxs, faceElementIdxs, tempIdxs
    integer(shortInt)                            :: i, idx, j, nIdxs, nTempIdxs

    ! Compute nIdxs and initialise nTempIdxs = 0
    nIdxs = size(idxs)
    nTempIdxs = 0
    
    ! Do a first pass and count the number of elements sharing each face.
    do i = 1, nIdxs
      nTempIdxs = nTempIdxs + size(self % shelf(idxs(i)) % item % getElementIdxs())

    end do
    allocate(tempIdxs(nTempIdxs))

    ! Do a second pass and populate tempIdxs.
    idx = 0
    do i = 1, nIdxs
      faceElementIdxs = self % shelf(idxs(i)) % item % getElementIdxs()
      do j = 1, size(faceElementIdxs)
        idx = idx + 1
        tempIdxs(idx) = faceElementIdxs(j)

      end do

    end do

    ! Now remove potential duplicates from the tempIdxs array.
    elementIdxs = removeDuplicates(tempIdxs)

  end function getFaceElementIdxs_shortIntArray

  !! Function 'getFaceHasElements'
  !!
  !! Basic description:
  !!   Returns .true. if a face in the shelf is already associated to elements.
  !!
  !! Arguments:
  !!   idx [in]    -> Index of the face in the shelf.
  !!
  !! Result:
  !!   hasElements -> .true. if the face in the shelf is associated to elements.
  !!
  elemental function getFaceHasElements(self, idx) result(hasElements)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    logical(defBool)              :: hasElements

    hasElements = self % shelf(idx) % item % getHasElements()

  end function getFaceHasElements

  !! Function 'getFaceIsBoundary'
  !!
  !! Basic description:
  !!   Returns .true. if a face in the shelf is a boundary face.
  !!
  !! Arguments:
  !!   idx [in]   -> Index of the face in the shelf.
  !!
  !! Result:
  !!   isBoundary -> .true. if the face in the shelf is a boundary face.
  !!
  elemental function getFaceIsBoundary(self, idx) result(isBoundary)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    logical(defBool)              :: isBoundary

    isBoundary = self % shelf(idx) % item % getIsBoundary()

  end function getFaceIsBoundary

  !! Function 'getFaceNormal'
  !!
  !! Basic description:
  !!   Returns the signed normal vector of a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the face in the shelf. Can be negative if the normal vector needs to be flipped.
  !!
  !! Result:
  !!   normal   -> 3-D coordinates of the face's signed normal vector.
  !!
  pure function getFaceNormal(self, idx) result(centroid)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    real(defReal), dimension(3)   :: centroid

    centroid = self % shelf(abs(idx)) % item % getNormal(idx)

  end function getFaceNormal

  !! Function 'getFaceTriangleIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the triangles in a face of the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the face in the shelf.
  !!
  !! Result:
  !!   triangleIdxs -> Indices of the triangles in the face.
  !!
  pure function getFaceTriangleIdxs(self, idx) result(triangleIdxs)
    class(faceShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: triangleIdxs

    triangleIdxs = self % shelf(idx) % item % getTriangleIdxs()

  end function getFaceTriangleIdxs

  !! Function 'getFaceType'
  !!
  !! Basic description:
  !!   Returns the type of a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the face in the shelf.
  !!
  !! Result:
  !!   type     -> Type of the face.
  !!
  pure function getFaceType(self, idx) result(type)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    character(:), allocatable     :: type

    type = self % shelf(idx) % item % getType()

  end function getFaceType

  !! Function 'getFaceVertexIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in a face of the shelf.
  !!
  !! Arguments:
  !!   idx [in]   -> Index of the face in the shelf.
  !!
  !! Result:
  !!   vertexIdxs -> Indices of the vertices in the face.
  !!
  pure function getFaceVertexIdxs(self, idx) result(vertexIdxs)
    class(faceShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: vertexIdxs

    vertexIdxs = self % shelf(idx) % item % getVertexIdxs()

  end function getFaceVertexIdxs
  
  !! Function 'getSize'
  !!
  !! Basic description:
  !!   Returns the size of the faceShelf.
  !!
  !! Result:
  !!   size -> Size of the faceShelf.
  !!
  elemental function getSize(self) result(nFaces)
    class(faceShelf), intent(in) :: self
    integer(shortInt)            :: nFaces
    
    nFaces = size(self % shelf)
    
  end function getSize

  !! Subroutine 'initFace'
  !!
  !! Basic description:
  !!   Initialises a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in]            -> Index of the face.
  !!   nInternalFaces [in] -> Number of internal faces in the shelf.
  !!
  subroutine initFace(self, idx, faceIdx, isBoundary, vertexIdxs, AB, AC, centroid, normal, area, type, boundingBox)
    class(faceShelf), intent(inout)             :: self
    integer(shortInt), intent(in)               :: idx, faceIdx
    logical(defBool), intent(in)                :: isBoundary
    integer(shortInt), dimension(:), intent(in) :: vertexIdxs
    real(defReal), dimension(3), intent(in)     :: AB, AC, centroid, normal
    real(defReal), intent(in)                   :: area
    character(*), intent(in)                    :: type
    type(axisAlignedBoundingBox), intent(in)    :: boundingBox

    ! Allocate face in the shelf and set everything.
    call self % allocateFace(idx, type)
    call self % shelf(idx) % item % init(idx, faceIdx, isBoundary, area, centroid, normal, AB, AC, vertexIdxs, type, boundingBox)

  end subroutine initFace

  !!
  !!
  !!
  elemental function intersectsFace_BoundingBox(self, idx, vertices, boundingBox) result(doesIt)
    class(faceShelf), intent(in)             :: self
    integer(shortInt), intent(in)            :: idx
    type(vertexShelf), intent(in)            :: vertices
    type(axisAlignedBoundingBox), intent(in) :: boundingBox
    logical(defBool)                         :: doesIt

    call self % shelf(idx) % item % intersects(vertices, boundingBox, doesIt)

  end function intersectsFace_BoundingBox

  !!
  !!
  !!
  elemental function intersectsFaceBoundingBox_BoundingBox(self, idx, boundingBox) result(doesIt)
    class(faceShelf), intent(in)             :: self
    integer(shortInt), intent(in)            :: idx
    type(axisAlignedBoundingBox), intent(in) :: boundingBox
    logical(defBool)                         :: doesIt
    class(face), allocatable                 :: item

    call self % shelf(idx) % item % intersectsBoundingBox(boundingbox, doesIt)

  end function intersectsFaceBoundingBox_BoundingBox
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(faceShelf), intent(inout) :: self
    integer(shortInt)               :: i
    
    if (allocated(self % shelf)) then
      do i = 1, size(self % shelf)
        if (allocated(self % shelf(i) % item)) deallocate(self % shelf(i) % item)

      end do
      deallocate(self % shelf)

    end if

  end subroutine kill

  !! Subroutine 'splitFace'
  !!
  !! Basic description:
  !!   Splits a face in the shelf into triangles.
  !!
  !! Arguments:
  !!   idx [in]                -> Index of the face in the shelf.
  !!   edges [inout]           -> An edgeShelf.
  !!   triangles [inout]       -> A triangleShelf.
  !!   vertices [inout]        -> A vertexShelf.
  !!   lastEdgeIdx [inout]     -> Index of the last edge in the edgeShelf.
  !!   lastTriangleIdx [inout] -> Index of the last triangle in the triangleShelf.
  !!
  subroutine splitFace(self, idx, edges, vertices, lastEdgeIdx, lastTriangleIdx, triangles)
    class(faceShelf), intent(inout)            :: self
    integer(shortInt), intent(in)              :: idx
    type(edgeShelf), intent(inout)             :: edges
    type(vertexShelf), intent(inout)           :: vertices
    integer(shortInt), intent(inout)           :: lastEdgeIdx, lastTriangleIdx
    type(faceBox), dimension(:), intent(inout) :: triangles

    call self % shelf(idx) % item % split(edges, vertices, lastEdgeIdx, lastTriangleIdx, triangles)

  end subroutine splitFace
  
end module faceShelf_class