module edge_class
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use genericProcedures,            only : append, areEqual
  use numPrecision
  use topologicalObject_inter,      only : topologicalObject, topologicalObjectBox, kill_super => kill
  use vertex_class,                 only : vertexBox
  
  implicit none
  private

  !!
  !!
  !!
  type, public :: buildEdgeInfo
    integer(shortInt)             :: idx = 0
    type(vertexBox), dimension(2) :: vertices
  end type buildEdgeInfo

  !!
  !!
  !!
  type, public :: edgeBox
    type(edge), pointer :: ptr => null()
  end type edgeBox
  
  !!
  !! Edge of a mesh linking two vertices.
  !!
  !! Private members:
  !!   idx            -> Index of the edge.
  !!   startVertexIdx -> Index of the first vertex in the edge.
  !!   endVertexIdx   -> Index of the end vertex in the edge.
  !!   edgeToFaces    -> Array that stores edge-to-faces connectivity information.
  !!   edgeToElements -> Array that stores edge-to-elements connectivity information.
  !!
  type, public, extends(topologicalObject)                :: edge
    private
    type(vertexBox), dimension(2)                         :: vertices
    real(defReal), dimension(3)                           :: edgeVector = ZERO, unitEdgeVector = ZERO
    integer(shortInt), dimension(:), allocatable          :: faceIdxs
    type(topologicalObjectBox), dimension(:), allocatable :: elements
    type(axisAlignedBoundingBox)                          :: boundingBox
  contains
    ! Build procedures.
    procedure :: addElement
    procedure :: addFaceIdx
    procedure :: init
    procedure :: kill
    ! Runtime procedures.
    procedure :: distanceSquared
    procedure :: getBoundingBoxBounds
    procedure :: getCentroid
    procedure :: getEdgeVector
    procedure :: getElements
    procedure :: getFaceIdxs
    procedure :: getVertices
    procedure :: intersects_BoundingBox
  end type edge

contains

  !! Subroutine 'addElementIdx'
  !!
  !! Basic description:
  !!   Adds the index of an element sharing the edge.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element.
  !!
  subroutine addElement(self, box)
    class(edge), intent(inout)                            :: self
    type(topologicalObjectBox), intent(in)                :: box
    integer(shortInt)                                     :: nElements
    type(topologicalObjectBox), dimension(:), allocatable :: tempElements

    if (allocated(self % elements)) then
      nElements = size(self % elements)
      allocate(tempElements(nElements + 1))
      tempElements(1:nElements) = self % elements
      tempElements(nElements + 1) = box
      call move_alloc(tempElements, self % elements)

    else
      allocate(self % elements(1))
      self % elements(1) = box

    end if

  end subroutine addElement

  !! Subroutine 'addFaceIdx'
  !!
  !! Basic description:
  !!   Adds the index of a face sharing the edge.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the face.
  !!
  elemental subroutine addFaceIdx(self, idx)
    class(edge), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx

    call append(self % faceIdxs, idx)

  end subroutine addFaceIdx

  !!
  !!
  !!
  function distanceSquared(self, r) result(dSquared)
    class(edge), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: dSquared
    real(defReal), dimension(3)             :: pointVector
    real(defReal)                           :: lSquared, t

    ! First pointVector and the square of the edge length.
    pointVector = r - self % vertices(1) % ptr % getCoordinates()
    lSquared = dot_product(self % edgeVector, self % edgeVector)

    ! Handle the case of a zero-length segment.
    if (areEqual(lSquared, ZERO)) then
        dSquared = dot_product(pointVector, pointVector)
        return
        
    end if

    ! Compute the normalisation parameter t by projecting pointVector onto edgeVector and
    ! snap it to the range [0, 1].
    t = max(ZERO, min(ONE, dot_product(pointVector, self % edgeVector) / lSquared))

    ! Now compute dSquared.
    pointVector = pointVector - self % edgeVector * t
    dSquared = dot_product(pointVector, pointVector)

  end function distanceSquared

  !!
  !!
  !!
  pure function getBoundingBoxBounds(self) result(bounds)
    class(edge), intent(in)        :: self
    real(defReal), dimension(3, 2) :: bounds

    bounds = self % boundingBox % getBounds()

  end function getBoundingBoxBounds

  !!
  !!
  !!
  pure function getCentroid(self) result(centroid)
    class(edge), intent(in)     :: self
    real(defReal), dimension(3) :: centroid

    centroid = HALF * (self % vertices(1) % ptr % getCoordinates() + self % vertices(2) % ptr % getCoordinates())

  end function getCentroid

  !!
  !!
  !!
  pure function getEdgeVector(self) result(edgeVector)
    class(edge), intent(in)     :: self
    real(defReal), dimension(3) :: edgeVector

    edgeVector = self % edgeVector

  end function getEdgeVector

  !! Function 'getElementIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the elements sharing the edge.
  !!
  !! Result:
  !!   elementIdxs -> Indices of the elements sharing the edge.
  !!
  function getElements(self) result(elements)
    class(edge), target, intent(in)                       :: self
    type(topologicalObjectBox), dimension(:), allocatable :: elements

    if (allocated(self % elements)) then
      elements = self % elements

    else
      allocate(elements(0))

    end if

  end function getElements

  !! Function 'getFaceIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the faces sharing the edge.
  !!
  !! Result:
  !!   faceIdxs -> Indices of the faces sharing the edge.
  !!
  pure function getFaceIdxs(self) result(faceIdxs)
    class(edge), intent(in)                             :: self
    integer(shortInt), dimension(size(self % faceIdxs)) :: faceIdxs

    faceIdxs = self % faceIdxs

  end function getFaceIdxs

  !! Function 'getVertexIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in the edge.
  !!
  !! Result:
  !!   vertexIdxs -> Indices of the vertices in the edge.
  !!
  function getVertices(self) result(boxes)
    class(edge), intent(in)       :: self
    type(vertexBox), dimension(2) :: boxes

    boxes = self % vertices

  end function getVertices

  !!
  !!
  !!
  subroutine init(self, info)
    class(edge), intent(inout)      :: self
    type(buildEdgeInfo), intent(in) :: info
    integer(shortInt), dimension(2) :: vertexIdxs
    integer(shortInt)               :: i
    real(defReal), dimension(3, 2)  :: allCoords

    ! Set information from payload.
    call self % setIdx(info % idx)

    ! Sort vertices according to their indices,
    do i = 1, 2
      vertexIdxs(i) = info % vertices(i) % ptr % getIdx()
      allCoords(:, i) = info % vertices(i) % ptr % getCoordinates()

    end do
    self % vertices(1) = info % vertices(minloc(vertexIdxs, 1))
    self % vertices(2) = info % vertices(maxloc(vertexIdxs, 1))
    
    self % edgeVector = self % vertices(2) % ptr % getCoordinates() - self % vertices(1) % ptr % getCoordinates()
    self % unitEdgeVector = self % edgeVector / norm2(self % edgeVector)

    ! Compute bounding box.
    call self % boundingBox % computeBounds(allCoords)

  end subroutine init

  !!
  !!
  !!
  elemental function intersects_BoundingBox(self, boundingBox) result(doesIt)
    class(edge), intent(in)                  :: self
    type(axisAlignedBoundingBox), intent(in) :: boundingBox
    logical(defBool)                         :: doesIt

    ! First check if bounding boxes overlap.
    doesIt = boundingBox % intersects(self % vertices(1) % ptr % getCoordinates(), self % vertices(2) % ptr % getCoordinates())

  end function intersects_BoundingBox

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an unitialised state.
  !!
  elemental subroutine kill(self)
    class(edge), intent(inout) :: self
    integer(shortInt)          :: i

    ! Superclass.
    call kill_super(self)
    
    ! Local.
    self % edgeVector = ZERO
    self % unitEdgeVector = ZERO
    do i = 1, 2
      nullify(self % vertices(i) % ptr)

    end do
    if (allocated(self % faceIdxs)) deallocate(self % faceIdxs)
    
    if (allocated(self % elements)) then
      do i = 1, size(self % elements)
        nullify(self % elements(i) % ptr)

      end do
      deallocate(self % elements)

    end if

  end subroutine kill

end module edge_class