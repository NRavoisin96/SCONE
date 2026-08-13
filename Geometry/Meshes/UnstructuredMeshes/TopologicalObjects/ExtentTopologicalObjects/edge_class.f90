module edge_class
  
  use axisAlignedBoundingBox_class,  only : axisAlignedBoundingBox
  use extentTopologicalObject_inter, only : buildExtentTopologicalObjectPayload, extentTopologicalObject, kill_super => kill
  use genericProcedures,             only : append, areEqual, fatalError, numToChar
  use numPrecision
  use publicObjects,                 only : intersectionTestResult, newIntersectionTestPayload
  use topologicalObject_inter,       only : buildTopologicalObjectPayload, topologicalObjectBox
  use vertex_class,                  only : vertexBox
  
  implicit none
  private

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
  type, public, extends(extentTopologicalObject)          :: edge
    private
    type(vertexBox), dimension(2)                         :: vertices
    real(defReal)                                         :: length
    real(defReal), dimension(3)                           :: unitEdgeVector = ZERO
    type(topologicalObjectBox), dimension(:), allocatable :: sharingElements, sharingFaces
  contains
    ! Build procedures.
    procedure :: addSharingElement
    procedure :: addSharingFace
    procedure :: build
    procedure :: connectComponents
    procedure :: kill
    ! Runtime procedures.
    procedure :: distanceSquared
    procedure :: getEdgeVector
    procedure :: getSharingElements
    procedure :: getSharingFaceIdxs
    procedure :: getSharingFaces
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
  subroutine addSharingElement(self, box)
    class(edge), intent(inout)                            :: self
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

  !! Subroutine 'addFaceIdx'
  !!
  !! Basic description:
  !!   Adds the index of a face sharing the edge.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the face.
  !!
  subroutine addSharingFace(self, box)
    class(edge), intent(inout)                            :: self
    type(topologicalObjectBox), intent(in)                :: box
    integer(shortInt)                                     :: nSharingFaces
    type(topologicalObjectBox), dimension(:), allocatable :: tempSharingFaces

    if (allocated(self % sharingFaces)) then
      nSharingFaces = size(self % sharingFaces)
      allocate(tempSharingFaces(nSharingFaces + 1))
      tempSharingFaces(1:nSharingFaces) = self % sharingFaces
      tempSharingFaces(nSharingFaces + 1) = box
      call move_alloc(tempSharingFaces, self % sharingFaces)

    else
      allocate(self % sharingFaces(1))
      self % sharingFaces(1) = box

    end if

  end subroutine addSharingFace

  !!
  !!
  !!
  pure function distanceSquared(self, r) result(dSquared)
    class(edge), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: dSquared
    real(defReal), dimension(3)             :: edgeVector, pointVector
    real(defReal)                           :: lSquared, t

    ! First pointVector and the square of the edge length.
    edgeVector = self % unitEdgeVector * self % length
    pointVector = r - self % vertices(1) % ptr % getCoordinates()
    lSquared = self % length * self % length

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

  end function distanceSquared

  !!
  !!
  !!
  pure function getEdgeVector(self) result(edgeVector)
    class(edge), intent(in)     :: self
    real(defReal), dimension(3) :: edgeVector

    edgeVector = self % unitEdgeVector * self % length

  end function getEdgeVector

  !! Function 'getElementIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the elements sharing the edge.
  !!
  !! Result:
  !!   elementIdxs -> Indices of the elements sharing the edge.
  !!
  function getSharingElements(self) result(sharingElements)
    class(edge), target, intent(in)                       :: self
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
    class(edge), intent(in)                      :: self
    integer(shortInt)                            :: i, nSharingFaces
    integer(shortInt), dimension(:), allocatable :: sharingFaceIdxs

    if(allocated(self % sharingFaces)) then
      nSharingFaces = size(self % sharingFaces)
      allocate(sharingFaceIdxs(nSharingFaces))

      do i = 1, nSharingFaces
        sharingFaceIdxs(i) = self % sharingFaces(i) % ptr % getIdx()

      end do

    else
      allocate(sharingFaceIdxs(0))

    end if

  end function getSharingFaceIdxs

  !! Function 'getFaceIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the faces sharing the edge.
  !!
  !! Result:
  !!   faceIdxs -> Indices of the faces sharing the edge.
  !!
  function getSharingFaces(self) result(sharingFaces)
    class(edge), target, intent(in)                       :: self
    type(topologicalObjectBox), dimension(:), allocatable :: sharingFaces

    if (allocated(self % sharingFaces)) then
      sharingFaces = self % sharingFaces

    else
      allocate(sharingFaces(0))

    end if

  end function getSharingFaces

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
  subroutine build(self, payload)
    class(edge), intent(inout)                          :: self
    class(buildTopologicalObjectPayload), intent(inout) :: payload
    type(buildExtentTopologicalObjectPayload), pointer  :: payloadPtr
    integer(shortInt), dimension(2)                     :: vertexIdxs
    integer(shortInt)                                   :: i
    real(defReal), dimension(3)                         :: edgeVector
    real(defReal)                                       :: length
    character(*), parameter                             :: here = 'build (edge_class.f90)'

    ! Downcast payload to correct type.
    select type(ptr => payload)
      type is(buildExtentTopologicalObjectPayload)
        payloadPtr => ptr

      class default
        call fatalError(here, 'Invalid payload type.')

    end select

    ! Sort vertices according to their indices.
    if (size(payloadPtr % vertices) /= 2) call fatalError(here, 'Invalid number of vertices.')
    allocate(payloadPtr % allCoords(3, 2))
    do i = 1, 2
      vertexIdxs(i) = payloadPtr % vertices(i) % ptr % getIdx()
      payloadPtr % allCoords(:, i) = payloadPtr % vertices(i) % ptr % getCoordinates()

    end do

    self % vertices(1) = payloadPtr % vertices(minloc(vertexIdxs, 1))
    self % vertices(2) = payloadPtr % vertices(maxloc(vertexIdxs, 1))
    
    edgeVector = self % vertices(2) % ptr % getCoordinates() - self % vertices(1) % ptr % getCoordinates()
    length = norm2(edgeVector)
    self % unitEdgeVector = edgeVector / length
    self % length = length
    payloadPtr % centroid = HALF * sum(payloadPtr % allCoords, 2)

  end subroutine build

  !!
  !!
  !!
  subroutine connectComponents(self)
    class(edge), target, intent(inout) :: self
    type(topologicalObjectBox)         :: box
    integer(shortInt)                  :: i

    box % ptr => self
    do i = 1, 2
      call self % vertices(i) % ptr % addSharingEdge(box)

    end do

  end subroutine connectComponents

  !!
  !!
  !!
  subroutine intersects_BoundingBox(self, boundingBox, doesIt)
    class(edge), intent(in)                  :: self
    type(axisAlignedBoundingBox), intent(in) :: boundingBox
    logical(defBool), intent(out)            :: doesIt
    type(intersectionTestResult)             :: result

    ! First check if bounding boxes overlap.
    call boundingBox % intersects(newIntersectionTestPayload(self % vertices(1) % ptr % getCoordinates(), &
                                                             self % unitEdgeVector, self % length), result)
    doesIt = result % intersects

  end subroutine intersects_BoundingBox

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
    self % length = ZERO
    self % unitEdgeVector = ZERO
    do i = 1, 2
      nullify(self % vertices(i) % ptr)

    end do
    
    if (allocated(self % sharingElements)) then
      do i = 1, size(self % sharingElements)
        nullify(self % sharingElements(i) % ptr)

      end do
      deallocate(self % sharingElements)

    end if

    if (allocated(self % sharingFaces)) then
      do i = 1, size(self % sharingFaces)
        nullify(self % sharingFaces(i) % ptr)

      end do
      deallocate(self % sharingFaces)

    end if

  end subroutine kill

end module edge_class