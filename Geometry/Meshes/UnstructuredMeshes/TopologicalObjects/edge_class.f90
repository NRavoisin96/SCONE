module edge_class
  
  use genericProcedures,       only : append, areEqual
  use numPrecision
  use topologicalObject_inter, only : topologicalObject, kill_super => kill
  use vertex_class,            only : vertexBox
  
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
  type, public, extends(topologicalObject)       :: edge
    private
    type(vertexBox), dimension(2)                :: vertices
    real(defReal), dimension(3)                  :: edgeVector = ZERO, unitEdgeVector = ZERO
    integer(shortInt), dimension(:), allocatable :: faceIdxs, elementIdxs
  contains
    ! Build procedures.
    procedure                                    :: addElementIdx
    procedure                                    :: addFaceIdx
    procedure                                    :: init
    procedure                                    :: kill
    ! Runtime procedures.
    procedure                                    :: distanceSquared
    procedure                                    :: getEdgeVector
    procedure                                    :: getElementIdxs
    procedure                                    :: getFaceIdxs
    procedure                                    :: getVertices
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
  elemental subroutine addElementIdx(self, idx)
    class(edge), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx

    call append(self % elementIdxs, idx, .true.)

  end subroutine addElementIdx

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
    class(edge), intent(in)     :: self
    real(defReal), dimension(3) :: r
    real(defReal)               :: dSquared
    real(defReal), dimension(3) :: pointVector
    real(defReal)               :: lSquared, t

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
  pure function getElementIdxs(self) result(faceIdxs)
    class(edge), intent(in)                                :: self
    integer(shortInt), dimension(size(self % elementIdxs)) :: faceIdxs

    faceIdxs = self % elementIdxs

  end function getElementIdxs

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
  subroutine init(self, idx, vertices)
    class(edge), intent(inout)                :: self
    integer(shortInt), intent(in)             :: idx
    type(vertexBox), dimension(2), intent(in) :: vertices

    call self % setIdx(idx)
    self % vertices = vertices
    
    self % edgeVector = vertices(2) % ptr % getCoordinates() - vertices(1) % ptr % getCoordinates()
    self % unitEdgeVector = self % edgeVector / norm2(self % edgeVector)

  end subroutine init

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
    if (allocated(self % elementIdxs)) deallocate(self % elementIdxs)

  end subroutine kill

end module edge_class