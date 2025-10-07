module edge_class
  
  use numPrecision
  use genericProcedures, only : append
  
  implicit none
  private
  
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
  type, public                                   :: edge
    private
    integer(shortInt)                            :: idx = 0
    integer(shortInt), dimension(2)              :: vertexIdxs = 0
    integer(shortInt), dimension(:), allocatable :: faceIdxs, elementIdxs, elementIdxsArray
    real(defReal), dimension(3)                  :: unitVector = ZERO, localBasis1 = ZERO, localBasis2 = ZERO, &
                                                    vector = ZERO
    real(defReal)                                :: length = ZERO, dotProductOfVector = ZERO
    real(defReal), dimension(:), allocatable     :: anglesArray
    !logical                                      :: isBoundary = .FALSE.

  contains

    ! Build procedures.
    procedure                                    :: addElementIdx
    procedure                                    :: addFaceIdx
    procedure                                    :: kill
    procedure                                    :: setIdx
    procedure                                    :: setVertexIdxs
    procedure                                    :: setUnitVector
    procedure                                    :: setVector
    procedure                                    :: setLocalBasis1
    procedure                                    :: setLocalBasis2
    procedure                                    :: setLength
    procedure                                    :: setDotProductOfVector
    procedure                                    :: setAnglesArray
    procedure                                    :: setElementIdxsArray
    procedure                                    :: isAllocatedAnglesArray
    !procedure                                    :: setIsBoundary
    ! Runtime procedures.
    procedure                                    :: getElementIdxs
    procedure                                    :: getFaceIdxs
    procedure                                    :: getIdx
    procedure                                    :: getVertexIdxs
    procedure                                    :: getUnitVector
    procedure                                    :: getVector
    procedure                                    :: getLocalBasis1
    procedure                                    :: getLocalBasis2
    procedure                                    :: getLength
    procedure                                    :: getDotProductOfVector
    procedure                                    :: getAnglesArray
    procedure                                    :: getElementIdxsArray
    !procedure                                    :: getIsBoundary
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

  !! Function 'getIdx'
  !!
  !! Basic description:
  !!   Returns the index of the edge.
  !!
  !! Result:
  !!   idx -> Index of the edge.
  !!
  elemental function getIdx(self) result(idx)
    class(edge), intent(in) :: self
    integer(shortInt)       :: idx

    idx = self % idx

  end function getIdx

  !! Function 'getVertexIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in the edge.
  !!
  !! Result:
  !!   vertexIdxs -> Indices of the vertices in the edge.
  !!
  pure function getVertexIdxs(self) result(vertexIdxs)
    class(edge), intent(in)         :: self
    integer(shortInt), dimension(2) :: vertexIdxs

    vertexIdxs = self % vertexIdxs

  end function getVertexIdxs

  !!
  !!
  !!
  pure function getUnitVector(self) result(unitVector)
    class(edge), intent(in)       :: self
    real(defReal), dimension(3)   :: unitVector

    unitVector = self % unitVector

  end function getUnitVector

  !!
  !!
  !!
  pure function getVector(self) result(vector)
    class(edge), intent(in)       :: self
    real(defReal), dimension(3)   :: vector

    vector = self % vector

  end function getVector

  !!
  !!
  !!
  pure function getLocalBasis1(self) result(localBasis1)
    class(edge), intent(in)       :: self
    real(defReal), dimension(3)   :: localBasis1

    localBasis1 = self % localBasis1

  end function getLocalBasis1

  !!
  !!
  !!
  pure function getLocalBasis2(self) result(localBasis2)
    class(edge), intent(in)       :: self
    real(defReal), dimension(3)   :: localBasis2

    localBasis2 = self % localBasis2

  end function getLocalBasis2

  !!
  !!
  !!
  elemental function getLength(self) result(length)
    class(edge), intent(in) :: self
    real(defReal)           :: length

    length = self % length

  end function getLength

  !!
  !!
  !!
  elemental function getDotProductOfVector(self) result(dotProductOfVector)
    class(edge), intent(in) :: self
    real(defReal)           :: dotProductOfVector

    dotProductOfVector = self % dotProductOfVector

  end function getDotProductOfVector

  !!
  !!
  !!
  pure function getAnglesArray(self) result(anglesArray)
    class(edge), intent(in)                    :: self
    real(defReal), dimension(:), allocatable   :: anglesArray

    anglesArray = self % anglesArray

  end function getAnglesArray

  !!
  !!
  !!
  pure function getElementIdxsArray(self) result(elementIdxsArray)
    class(edge), intent(in)                        :: self
    integer(shortInt), dimension(:), allocatable   :: elementIdxsArray

    elementIdxsArray = self % elementIdxsArray

  end function getElementIdxsArray

  ! !!
  ! !!
  ! !!
  ! elemental function getIsBoundary(self) result(isBoundary)
  !   class(edge), intent(in)                        :: self
  !   logical                                        :: isBoundary

  !   isBoundary = self % isBoundary

  ! end function getIsBoundary


  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an unitialised state.
  !!
  elemental subroutine kill(self)
    class(edge), intent(inout) :: self

    self % idx = 0
    self % vertexIdxs = 0
    if (allocated(self % faceIdxs)) deallocate(self % faceIdxs)
    if (allocated(self % elementIdxs)) deallocate(self % elementIdxs)

  end subroutine kill

  !! Subroutine 'setIdx'
  !!
  !! Basic description:
  !!   Sets the index of the edge.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the edge.
  !!
  elemental subroutine setIdx(self, idx)
    class(edge), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx

    self % idx = idx

  end subroutine setIdx

  !! Subroutine 'setVertexIdxs'
  !!
  !! Basic description:
  !!   Sets the indices of the vertices in the edge.
  !!
  !! Arguments:
  !!   vertexIdxs [in] -> Array containing the indices of the vertices in the edge.
  !!
  pure subroutine setVertexIdxs(self, vertexIdxs)
    class(edge), intent(inout)                  :: self
    integer(shortInt), dimension(2), intent(in) :: vertexIdxs

    self % vertexIdxs = vertexIdxs

  end subroutine setVertexIdxs

  !!
  !!
  !!
  pure subroutine setUnitVector(self, unitVector)
    class(edge), intent(inout)               :: self
    real(defReal), intent(in), dimension(3)  :: unitVector

    self % unitVector = unitVector

  end subroutine setUnitVector

  !!
  !!
  !!
  pure subroutine setVector(self, vector)
    class(edge), intent(inout)               :: self
    real(defReal), intent(in), dimension(3)  :: vector

    self % vector = vector

  end subroutine setVector

  !!
  !!
  !!
  pure subroutine setLocalBasis1(self, localBasis1)
    class(edge), intent(inout)               :: self
    real(defReal), intent(in), dimension(3)  :: localBasis1

    self % localBasis1 = localBasis1

  end subroutine setlocalBasis1

  !!
  !!
  !!
  pure subroutine setLocalBasis2(self, localBasis2)
    class(edge), intent(inout)               :: self
    real(defReal), intent(in), dimension(3)  :: localBasis2

    self % localBasis2 = localBasis2

  end subroutine setLocalBasis2

  !!
  !!
  !!
  elemental subroutine setLength(self, length)
    class(edge), intent(inout) :: self
    real(defReal), intent(in)  :: length

    self % length = length

  end subroutine setLength

  !!
  !!
  !!
  elemental subroutine setDotProductOfVector(self, DotProductOfVector)
    class(edge), intent(inout) :: self
    real(defReal), intent(in)  :: dotProductOfVector

    self % dotProductOfVector = dotProductOfVector

  end subroutine setDotProductOfVector

  !!
  !!
  !!
  pure subroutine setAnglesArray(self, anglesArray)
    class(edge), intent(inout)               :: self
    real(defReal), intent(in), dimension(:)  :: anglesArray

    self % anglesArray = anglesArray

  end subroutine setAnglesArray

  !!
  !!
  !!
  pure subroutine setElementIdxsArray(self, elementIdxsArray)
    class(edge), intent(inout)                   :: self
    integer(shortInt), intent(in), dimension(:)  :: elementIdxsArray

    self % elementIdxsArray = elementIdxsArray

  end subroutine setElementIdxsArray

  !!
  !!
  !!
  elemental function isAllocatedAnglesArray(self) result(isAllocated)
    class(edge), intent(in)                      :: self
    logical                                      :: isAllocated

    isAllocated = allocated(self % anglesArray)

  end function isAllocatedAnglesArray

  ! !!
  ! !! no need to know if isBoundary == .TRUE. because if this is called, isBoundary == .TRUE.
  ! !! Otherwise, we keep isBoundary == .FALSE. from the initialisation
  ! elemental subroutine setIsBoundary(self)
  !   class(edge), intent(inout)                      :: self

  !   self % isBoundary = .TRUE.

  ! end subroutine setIsBoundary

end module edge_class