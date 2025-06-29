module edgeShelf_class
  
  use edge_class,                   only : edge, edgeBox
  use edgeFactory_func,             only : newEdgeBox
  use genericProcedures,            only : fatalError, numToChar
  use numPrecision
  use topologicalObject_inter,      only : topologicalObjectBox
  use topologicalObjectShelf_inter, only : topologicalObjectShelf
  use vertex_class,                 only : vertexBox
  
  implicit none
  private
  
  type, public, extends(topologicalObjectShelf) :: edgeShelf
    private
  contains
    procedure                             :: addElementIdxToEdge
    procedure                             :: addFaceIdxToEdge
    generic                               :: getEdgeBox => getEdgeBox_shortInt, getEdgeBox_shortIntArray
    procedure, private                    :: getEdgeBox_shortInt
    procedure, private                    :: getEdgeBox_shortIntArray
    procedure                             :: getEdgeElementIdxs
    procedure                             :: getEdgeFaceIdxs
    procedure                             :: initEdge
  end type edgeShelf

contains

  !! Subroutine 'addElementIdxToEdge'
  !!
  !! Basic description:
  !!   Adds the index of an element to an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]        -> Index of the edge in the shelf.
  !!   elementIdx [in] -> Index of the element containing the edge.
  !!
  subroutine addElementIdxToEdge(self, idx, elementIdx)
    class(edgeShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx, elementIdx
    type(edgeBox)                   :: box

    box = self % getEdgeBox(idx)
    call box % ptr % addElementIdx(elementIdx)

  end subroutine addElementIdxToEdge

  !! Subroutine 'addFaceIdxToEdge'
  !!
  !! Basic description:
  !!   Adds the index of a face to an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the edge in the shelf.
  !!   faceIdx [in] -> Index of the face containing the edge.
  !!
  subroutine addFaceIdxToEdge(self, idx, faceIdx)
    class(edgeShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx, faceIdx
    type(edgeBox)                   :: box

    box = self % getEdgeBox(idx)
    call box % ptr % addFaceIdx(faceIdx)

  end subroutine addFaceIdxToEdge

  !!
  !!
  !!
  function getEdgeBox_shortInt(self, idx) result(box)
    class(edgeShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    type(edgeBox)                 :: box
    type(topologicalObjectBox)    :: objectBox
    character(*), parameter       :: here = 'getEdgeBox_shortInt (edgeShelf_class.f90)'

    ! First get a pointer to a polymorphic topological object from the shelf.
    objectBox = self % getObjectBox(idx)
    if (.not. associated(objectBox % ptr)) call fatalError(here, 'Invalid pointer for edge with index: '//numToChar(idx)//'.')

    select type(ptr => objectBox % ptr)
      type is (edge)
        box % ptr => ptr

      class default
        ! Should never happen.
        call fatalError(here, 'Object in edgeShelf with idx: '//numToChar(idx)//' is not an edge.')

    end select

  end function getEdgeBox_shortInt

  !!
  !!
  !!
  function getEdgeBox_shortIntArray(self, idxs) result(boxes)
    class(edgeShelf), intent(in)                      :: self
    integer(shortInt), dimension(:), intent(in)       :: idxs
    type(edgeBox), dimension(size(idxs))              :: boxes
    type(topologicalObjectBox), dimension(size(idxs)) :: objectBoxes
    integer(shortInt)                                 :: i
    character(*), parameter                           :: here = 'getEdgeBox_shortIntArray (edgeShelf_class.f90)'

    objectBoxes = self % getObjectBox(idxs)
    do i = 1, size(idxs)
      if (.not. associated(objectBoxes(i) % ptr)) &
      call fatalError(here, 'Invalid pointer for edge with index: '//numToChar(idxs(i))//'.')

      select type(ptr => objectBoxes(i) % ptr)
        type is (edge)
          boxes(i) % ptr => ptr

        class default
          ! Should never happen.
          call fatalError(here, 'Object in edgeShelf with idx: '//numToChar(idxs(i))//' is not an edge.')

      end select

    end do

  end function getEdgeBox_shortIntArray

  !! Function 'getEdgeElementIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the elements sharing an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]    -> Index of the edge in the shelf.
  !!
  !! Result:
  !!   elementIdxs -> Indices of the elements sharing the edge.
  !!
  function getEdgeElementIdxs(self, idx) result(elementIdxs)
    class(edgeShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: elementIdxs
    type(edgeBox)                                :: box

    box = self % getEdgeBox(idx)
    elementIdxs = box % ptr % getElementIdxs()

  end function getEdgeElementIdxs

  !! Function 'getEdgeFaceIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the faces sharing an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the edge in the shelf.
  !!
  !! Result:
  !!   faceIdxs -> Indices of the faces sharing the edge.
  !!
  function getEdgeFaceIdxs(self, idx) result(faceIdxs)
    class(edgeShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: faceIdxs
    type(edgeBox)                                :: box

    box = self % getEdgeBox(idx)
    faceIdxs = box % ptr % getFaceIdxs()

  end function getEdgeFaceIdxs

  !! Subroutine 'initEdge'
  !!
  !! Basic description:
  !!   Initialises an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]        -> Index of the edge in the shelf.
  !!   vertexIdxs [in] -> Indices of the vertices in the edge.
  !!
  subroutine initEdge(self, idx, vertices)
    class(edgeShelf), intent(inout)           :: self
    integer(shortInt), intent(in)             :: idx
    type(vertexBox), dimension(2), intent(in) :: vertices
    type(edgeBox)                             :: box

    call newEdgeBox(idx, vertices, box)
    call self % addObject(idx, box % ptr)

  end subroutine initEdge
  
end module edgeShelf_class