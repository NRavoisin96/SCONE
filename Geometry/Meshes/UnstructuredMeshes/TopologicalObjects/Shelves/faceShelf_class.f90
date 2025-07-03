module faceShelf_class
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use coord_class,                  only : coord
  use edge_class,                   only : edgeBox
  use face_class,                   only : buildFaceInfo, face, faceBox
  use faceFactory_func,             only : newFaceBox
  use genericProcedures,            only : fatalError, numToChar, removeDuplicates
  use numPrecision
  use topologicalObject_inter,      only : topologicalObjectBox
  use topologicalObjectShelf_inter, only : topologicalObjectShelf
  use vertex_class,                 only : vertexBox
  
  implicit none
  private
  
  !!
  !! Storage space for faces of an OpenFOAM mesh.
  !!
  !! Private members:
  !!   shelf -> Array to store faces.
  !!
  type, public, extends(topologicalObjectShelf) :: faceShelf
    private
  contains
    procedure          :: addEdgeToFace
    procedure          :: addElementIdxToFace
    procedure          :: addVertexToFace
    procedure          :: computeFaceIntersection
    procedure          :: distanceSquaredFromFace
    procedure          :: getAllFaceBoundingBoxes
    procedure          :: getAllFaceCentroids
    procedure          :: getFaceArea
    generic            :: getFaceBox => getFaceBox_shortInt, getFaceBox_shortIntArray
    procedure, private :: getFaceBox_shortInt
    procedure, private :: getFaceBox_shortIntArray
    procedure          :: getFaceBoundingBox
    procedure          :: getFaceCentroid
    procedure          :: getFaceEdges
    generic            :: getFaceElementIdxs => getFaceElementIdxs_shortInt, getFaceElementIdxs_shortIntArray
    procedure, private :: getFaceElementIdxs_shortInt
    procedure, private :: getFaceElementIdxs_shortIntArray
    procedure          :: getFaceHasElements
    procedure          :: getFaceIsBoundary
    procedure          :: getFaceNormal
    procedure          :: getFaceTriangleIdxs
    procedure          :: getFaceType
    procedure          :: getFaceVertices
    procedure          :: init
    procedure          :: initFace
    generic            :: intersectsFace => intersectsFace_BoundingBox
    procedure, private :: intersectsFace_BoundingBox
    generic            :: intersectsFaceBoundingBox => intersectsFaceBoundingBox_BoundingBox
    procedure, private :: intersectsFaceBoundingBox_BoundingBox
    procedure          :: splitFace
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
  subroutine addEdgeToFace(self, idx, edge)
    class(faceShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx
    type(edgeBox), intent(in)       :: edge
    type(faceBox)                   :: box

    box = self % getFaceBox(idx)
    call box % ptr % addEdge(edge)

  end subroutine addEdgeToFace

  !! Subroutine 'addElementIdxToFace'
  !!
  !! Basic description:
  !!   Adds the index of an element to a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in]        -> Index of the face in the shelf.
  !!   elementIdx [in] -> Index of the element containing the face.
  !!
  subroutine addElementIdxToFace(self, idx, elementIdx)
    class(faceShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx, elementIdx
    type(faceBox)                   :: box

    box = self % getFaceBox(idx)
    call box % ptr % addElementIdx(elementIdx)

  end subroutine addElementIdxToFace

  !! Subroutine 'addVertexIdxToFace'
  !!
  !! Basic description:
  !!   Adds the index of a vertex to a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in]       -> Index of the face in the shelf.
  !!   vertexIdx [in] -> Index of the vertex in the face.
  !!
  subroutine addVertexToFace(self, idx, vertex)
    class(faceShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx
    type(vertexBox), intent(in)     :: vertex
    type(faceBox)                   :: box

    box = self % getFaceBox(idx)
    call box % ptr % addVertex(vertex)

  end subroutine addVertexToFace

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
  subroutine computeFaceIntersection(self, idx, coords, d)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    type(coord), intent(in)       :: coords
    real(defReal), intent(out)    :: d
    type(faceBox)                 :: box

    box = self % getFaceBox(idx)
    call box % ptr % computeIntersection(coords, d)

  end subroutine computeFaceIntersection

  !! Function 'distanceSquaredFromFace'
  !!
  !!
  function distanceSquaredFromFace(self, idx, r) result(dSquared)
    class(faceShelf), intent(in)            :: self
    integer(shortInt), intent(in)           :: idx
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: dSquared
    type(faceBox)                           :: box

    box = self % getFaceBox(idx)
    dSquared = box % ptr % distanceSquared(r)

  end function distanceSquaredFromFace

  !! Function 'getAllFaceBoundingBoxes'
  !!
  !! Basic description:
  !!   Returns the bounding boxes of all the faces in the shelf.
  !!
  !! Results:
  !!   boundingBoxes -> A defReal array containing the bounding boxes of all the faces in the shelf.
  !!
  function getAllFaceBoundingBoxes(self) result(boundingBoxes)
    class(faceShelf), intent(in)                              :: self
    type(axisAlignedBoundingBox), dimension(self % getSize()) :: boundingBoxes
    integer(shortInt)                                         :: i
    type(faceBox)                                             :: box

    do i = 1, self % getSize()
      box = self % getFaceBox(i)
      boundingBoxes(i) = box % ptr % getBoundingBox()

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
  function getAllFaceCentroids(self) result(centroids)
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
  function getFaceArea(self, idx) result(area)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    real(defReal)                 :: area
    type(faceBox)                 :: box

    box = self % getFaceBox(idx)
    area = box % ptr % getArea()

  end function getFaceArea

  !!
  !!
  !!
  function getFaceBox_shortInt(self, idx) result(box)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    type(faceBox)                 :: box
    type(topologicalObjectBox)    :: objectBox
    character(*), parameter       :: here = 'getFaceBox_shortInt (faceShelf_class.f90)'

    ! First get a pointer to a polymorphic topological object from the shelf.
    objectBox = self % getObjectBox(idx)
    if (.not. associated(objectBox % ptr)) call fatalError(here, 'Invalid pointer for face with index: '//numToChar(idx)//'.')

    select type(ptr => objectBox % ptr)
      type is (face)
        box % ptr => ptr

      class default
        ! Should never happen.
        call fatalError(here, 'Object in faceShelf with idx: '//numToChar(idx)//' is not a face.')

    end select


  end function getFaceBox_shortInt

  !!
  !!
  !!
  function getFaceBox_shortIntArray(self, idxs) result(boxes)
    class(faceShelf), intent(in)                      :: self
    integer(shortInt), dimension(:), intent(in)       :: idxs
    type(faceBox), dimension(size(idxs))              :: boxes
    type(topologicalObjectBox), dimension(size(idxs)) :: objectBoxes
    integer(shortInt)                                 :: i
    character(*), parameter                           :: here = 'getFaceBox_shortIntArray (faceShelf_class.f90)'

    objectBoxes = self % getObjectBox(idxs)
    do i = 1, size(idxs)
      if (.not. associated(objectBoxes(i) % ptr)) &
      call fatalError(here, 'Invalid pointer for face with index: '//numToChar(idxs(i))//'.')

      select type(ptr => objectBoxes(i) % ptr)
        type is (face)
          boxes(i) % ptr => ptr

        class default
          ! Should never happen.
          call fatalError(here, 'Object in faceShelf with idx: '//numToChar(idxs(i))//' is not a face.')

      end select

    end do

  end function getFaceBox_shortIntArray

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
  function getFaceBoundingBox(self, idx) result(boundingBox)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    type(axisAlignedBoundingBox)  :: boundingBox
    type(faceBox)                 :: box

    box = self % getFaceBox(idx)
    boundingBox = box % ptr % getBoundingBox()

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
  function getFaceCentroid(self, idx) result(centroid)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    real(defReal), dimension(3)   :: centroid
    type(faceBox)                 :: box

    box = self % getFaceBox(idx)
    centroid = box % ptr % getCentroid()

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
  function getFaceEdges(self, idx) result(edges)
    class(faceShelf), intent(in)             :: self
    integer(shortInt), intent(in)            :: idx
    type(edgeBox), dimension(:), allocatable :: edges
    type(faceBox)                            :: box

    box = self % getFaceBox(idx)
    edges = box % ptr % getEdges()

  end function getFaceEdges

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
  function getFaceElementIdxs_shortInt(self, idx) result(elementIdxs)
    class(faceShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: elementIdxs
    type(faceBox)                                :: box

    box = self % getFaceBox(idx)
    elementIdxs = box % ptr % getElementIdxs()

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
  function getFaceElementIdxs_shortIntArray(self, idxs) result(elementIdxs)
    class(faceShelf), intent(in)                 :: self
    integer(shortInt), dimension(:), intent(in)  :: idxs
    integer(shortInt), dimension(:), allocatable :: elementIdxs, faceElementIdxs, tempIdxs
    integer(shortInt)                            :: i, idx, j, nIdxs, nTempIdxs
    type(faceBox), dimension(size(idxs))         :: boxes

    ! Compute nIdxs and initialise nTempIdxs = 0
    nIdxs = size(idxs)
    nTempIdxs = 0
    
    ! Do a first pass and count the number of elements sharing each face.
    boxes = self % getFaceBox(idxs)
    do i = 1, nIdxs
      nTempIdxs = nTempIdxs + size(boxes(i) % ptr % getElementIdxs())

    end do
    allocate(tempIdxs(nTempIdxs))

    ! Do a second pass and populate tempIdxs.
    idx = 0
    do i = 1, nIdxs
      faceElementIdxs = boxes(i) % ptr % getElementIdxs()
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
  function getFaceHasElements(self, idx) result(hasElements)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    logical(defBool)              :: hasElements
    type(faceBox)                 :: box

    box = self % getFaceBox(idx)
    hasElements = box % ptr % getHasElements()

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
  function getFaceIsBoundary(self, idx) result(isBoundary)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    logical(defBool)              :: isBoundary
    type(faceBox)                 :: box

    box = self % getFaceBox(idx)
    isBoundary = box % ptr % getIsBoundary()

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
  function getFaceNormal(self, idx) result(normal)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    real(defReal), dimension(3)   :: normal
    type(faceBox)                 :: box

    box = self % getFaceBox(abs(idx))
    normal = box % ptr % getNormal(idx)

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
  function getFaceTriangleIdxs(self, idx) result(triangleIdxs)
    class(faceShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: triangleIdxs
    type(faceBox)                                :: box

    box = self % getFaceBox(idx)
    triangleIdxs = box % ptr % getTriangleIdxs()

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
  function getFaceType(self, idx) result(type)
    class(faceShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    character(:), allocatable     :: type
    type(faceBox)                 :: box

    box = self % getFaceBox(idx)
    type = box % ptr % getType()

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
  function getFaceVertices(self, idx) result(vertices)
    class(faceShelf), intent(in)               :: self
    integer(shortInt), intent(in)              :: idx
    type(vertexBox), dimension(:), allocatable :: vertices
    type(faceBox)                              :: box

    box = self % getFaceBox(idx)
    vertices = box % ptr % getVertices()

  end function getFaceVertices

  !!
  !!
  !!
  subroutine init(self, infos)
    class(faceShelf), intent(inout)               :: self
    type(buildFaceInfo), dimension(:), intent(in) :: infos
    type(faceBox)                                 :: box
    integer(shortInt)                             :: i, nFaces

    ! Allocate shelf then populate it.
    nFaces = size(infos)
    call self % allocateShelf(nFaces)
    do i = 1, nFaces
      call newFaceBox(infos(i), box)
      call self % addObject(box % ptr)

    end do

  end subroutine init

  !! Subroutine 'initFace'
  !!
  !! Basic description:
  !!   Initialises a face in the shelf.
  !!
  !! Arguments:
  !!   idx [in]            -> Index of the face.
  !!   nInternalFaces [in] -> Number of internal faces in the shelf.
  !!
  subroutine initFace(self, info)
    class(faceShelf), intent(inout) :: self
    type(buildFaceInfo), intent(in) :: info
    type(faceBox)                   :: box

    call newFaceBox(info, box)
    call self % addObject(box % ptr)

  end subroutine initFace

  !!
  !!
  !!
  function intersectsFace_BoundingBox(self, idx, boundingBox) result(doesIt)
    class(faceShelf), intent(in)             :: self
    integer(shortInt), intent(in)            :: idx
    type(axisAlignedBoundingBox), intent(in) :: boundingBox
    logical(defBool)                         :: doesIt
    type(faceBox)                            :: box

    box = self % getFaceBox(idx)
    call box % ptr % intersects(boundingBox, doesIt)

  end function intersectsFace_BoundingBox

  !!
  !!
  !!
  function intersectsFaceBoundingBox_BoundingBox(self, idx, boundingBox) result(doesIt)
    class(faceShelf), intent(in)             :: self
    integer(shortInt), intent(in)            :: idx
    type(axisAlignedBoundingBox), intent(in) :: boundingBox
    logical(defBool)                         :: doesIt
    type(faceBox)                            :: box

    box = self % getFaceBox(idx)
    call box % ptr % intersectsBoundingBox(boundingbox, doesIt)

  end function intersectsFaceBoundingBox_BoundingBox

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
  subroutine splitFace(self, idx)
    class(faceShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx
    type(faceBox)                   :: box

    box = self % getFaceBox(idx)
    call box % ptr % split()

  end subroutine splitFace
  
end module faceShelf_class