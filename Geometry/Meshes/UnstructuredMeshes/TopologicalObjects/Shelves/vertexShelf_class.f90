module vertexShelf_class
  
  use numPrecision
  use genericProcedures,            only : append, fatalError, findCommon, numToChar
  use topologicalObject_inter,      only : topologicalObject, topologicalObjectBox
  use topologicalObjectShelf_inter, only : topologicalObjectShelf, kill_super => kill
  use vertex_class,                 only : vertex, vertexBox
  use vertexFactory_func,           only : newVertexBox
  
  implicit none
  private
  
  !!
  !! Storage space for vertices of a given OpenFOAM mesh.
  !!
  !! Public members:
  !!   shelf               -> Array to store vertices.
  !!   offset              -> User-supplied 3-D offset applied to the coordinates 
  !!                          of all the vertices in the shelf.
  !!   extremalCoordinates -> Array of minimum and maximum x-, y- and z-
  !!                          coordinates in the shelf.
  !!
  type, public, extends(topologicalObjectShelf) :: vertexShelf
    private
    real(defReal), dimension(6)                 :: extremalCoordinates = ZERO
  contains
    generic            :: addEdgeIdxToVertex => addEdgeIdxToVertex_shortInt, addEdgeIdxToVertex_shortIntArray
    procedure, private :: addEdgeIdxToVertex_shortInt
    procedure, private :: addEdgeIdxToVertex_shortIntArray
    generic            :: addElementToVertex => addElementToVertex_shortInt, addElementToVertex_shortIntArray
    procedure, private :: addElementToVertex_shortInt
    procedure, private :: addElementToVertex_shortIntArray
    generic            :: addFaceIdxToVertex => addFaceIdxToVertex_shortInt, addFaceIdxToVertex_shortIntArray
    procedure, private :: addFaceIdxToVertex_shortInt
    procedure, private :: addFaceIdxToVertex_shortIntArray
    procedure          :: findCommonEdgeIdx
    procedure          :: findCommonFaceIdx
    procedure          :: getAllCoordinates
    procedure          :: getExtremalCoordinates
    generic            :: getVertexCoordinates => getVertexCoordinates_shortInt, getVertexCoordinates_shortIntArray
    procedure, private :: getVertexCoordinates_shortInt
    procedure, private :: getVertexCoordinates_shortIntArray
    procedure          :: getVertexEdgeIdxs
    procedure          :: getVertexElements
    generic            :: getVertexFaceIdxs => getVertexFaceIdxs_shortInt, getVertexFaceIdxs_shortIntArray
    procedure, private :: getVertexFaceIdxs_shortInt
    procedure, private :: getVertexFaceIdxs_shortIntArray
    generic            :: getVertexBox => getVertexBox_shortInt, getVertexBox_shortIntArray
    procedure, private :: getVertexBox_shortInt
    procedure, private :: getVertexBox_shortIntArray
    procedure          :: init
    procedure          :: initVertex
    procedure          :: kill
    procedure          :: setExtremalCoordinates
  end type 

contains
  !! Subroutine 'addEdgeIdxToVertex'
  !!
  !! Basic description:
  !!   Adds the index of an edge to a vertex in the shelf.
  !!
  !! Arguments:
  !!   vertexIdx [in] -> Index of the vertex in the shelf.
  !!   edgeIdx [in]   -> Index of the edge containing the vertex.
  !!
  subroutine addEdgeIdxToVertex_shortInt(self, vertexIdx, edgeIdx)
    class(vertexShelf), intent(inout) :: self
    integer(shortInt), intent(in)     :: vertexIdx, edgeIdx
    type(vertexBox)                   :: box

    box = self % getVertexBox(vertexIdx)
    call box % ptr % addEdgeIdx(edgeIdx)

  end subroutine addEdgeIdxToVertex_shortInt

  !!
  !!
  !!
  subroutine addEdgeIdxToVertex_shortIntArray(self, vertexIdxs, edgeIdx)
    class(vertexShelf), intent(inout)            :: self
    integer(shortInt), dimension(:), intent(in)  :: vertexIdxs
    integer(shortInt), intent(in)                :: edgeIdx
    type(vertexBox), dimension(size(vertexIdxs)) :: boxes
    integer(shortInt)                            :: i

    boxes = self % getVertexBox(vertexIdxs)
    do i = 1, size(vertexIdxs)
      call boxes(i) % ptr % addEdgeIdx(edgeIdx)

    end do

  end subroutine addEdgeIdxToVertex_shortIntArray

  !! Subroutine 'addElementIdxToVertex'
  !!
  !! Basic description:
  !!   Adds the index of an element to a vertex in the shelf.
  !!
  !! Arguments:
  !!   vertexIdx [in]  -> Index of the vertex in the shelf.
  !!   elementIdx [in] -> Index of the element containing the vertex.
  !!
  subroutine addElementToVertex_shortInt(self, idx, element)
    class(vertexShelf), intent(inout)      :: self
    integer(shortInt), intent(in)          :: idx
    type(topologicalObjectBox), intent(in) :: element
    type(vertexBox)                        :: box

    box = self % getVertexBox(idx)
    call box % ptr % addElement(element)

  end subroutine addElementToVertex_shortInt

  !!
  !!
  !!
  subroutine addElementToVertex_shortIntArray(self, idxs, element)
    class(vertexShelf), intent(inout)           :: self
    integer(shortInt), dimension(:), intent(in) :: idxs
    type(topologicalObjectBox), intent(in)      :: element
    type(vertexBox), dimension(size(idxs))      :: boxes
    integer(shortInt)                           :: i

    boxes = self % getVertexBox(idxs)
    do i = 1, size(idxs)
      call boxes(i) % ptr % addElement(element)

    end do

  end subroutine addElementToVertex_shortIntArray

  !! Subroutine 'addFaceIdxToVertex'
  !!
  !! Basic description:
  !!   Adds the index of a face to a vertex in the shelf.
  !!
  !! Arguments:
  !!   vertexIdx [in] -> Index of the vertex in the shelf.
  !!   faceIdx [in]   -> Index of the face containing the vertex.
  !!
  subroutine addFaceIdxToVertex_shortInt(self, vertexIdx, faceIdx)
    class(vertexShelf), intent(inout) :: self
    integer(shortInt), intent(in)     :: vertexIdx, faceIdx
    type(vertexBox)                   :: box

    box = self % getVertexBox(vertexIdx)
    call box % ptr % addFaceIdx(faceIdx)

  end subroutine addFaceIdxToVertex_shortInt

  !!
  !!
  !!
  subroutine addFaceIdxToVertex_shortIntArray(self, vertexIdxs, faceIdx)
    class(vertexShelf), intent(inout)            :: self
    integer(shortInt), dimension(:), intent(in)  :: vertexIdxs
    integer(shortInt), intent(in)                :: faceIdx
    type(vertexBox), dimension(size(vertexIdxs)) :: boxes
    integer(shortInt)                            :: i

    boxes = self % getVertexBox(vertexIdxs)
    do i = 1, size(vertexIdxs)
      call boxes(i) % ptr % addFaceIdx(faceIdx)

    end do

  end subroutine addFaceIdxToVertex_shortIntArray

  !! Function 'findCommonEdgeIdx'
  !!
  !! Basic description:
  !!   Finds the index of the edge containing two vertices.
  !!
  !! Arguments:
  !!   firstVertexIdx [in]  -> Index of the first vertex in the edge.
  !!   secondVertexIdx [in] -> Index of the second vertex in the edge.
  !!
  !! Result:
  !!   edgeIdx              -> Index of the edge containing the two vertices.
  !!
  function findCommonEdgeIdx(self, firstVertexIdx, secondVertexIdx) result(edgeIdx)
    class(vertexShelf), intent(in)               :: self
    integer(shortInt), intent(in)                :: firstVertexIdx, secondVertexIdx
    type(vertexBox), dimension(2)                :: boxes
    integer(shortInt)                            :: edgeIdx, i
    integer(shortInt), dimension(:), allocatable :: commonIdxs

    ! Initialise edgeIdx = 0 and return immediately if any vertices are not associated with edges.
    edgeIdx = 0
    boxes = self % getVertexBox([firstVertexIdx, secondVertexIdx])
    do i = 1, 2
      if (.not. boxes(i) % ptr % hasEdges()) return

    end do
    
    ! Find common edge indices. Update edgeIdx only if common indices have been found.
    commonIdxs = findCommon(boxes(1) % ptr % getEdgeIdxs(), boxes(2) % ptr % getEdgeIdxs())
    if (size(commonIdxs) > 0) edgeIdx = commonIdxs(1)

  end function findCommonEdgeIdx

  !! Function 'findCommonTriangleIdx'
  !!
  !! Basic description:
  !!   Finds the index of the triangle containing three vertices.
  !!
  !! Arguments:
  !!   vertexIdxs [in] -> Indices of the vertices.
  !!
  !! Result:
  !!   triangleIdx     -> Index of the triangle containing the two vertices.
  !!
  function findCommonFaceIdx(self, vertexIdxs) result(faceIdx)
    class(vertexShelf), intent(in)               :: self
    integer(shortInt), dimension(3), intent(in)  :: vertexIdxs
    type(vertexBox), dimension(3)                :: boxes
    integer(shortInt)                            :: faceIdx, i
    integer(shortInt), dimension(:), allocatable :: commonIdxs

    ! Initialise triangleIdx = 0 and return immediately if any vertices are not associated with triangles.
    faceIdx = 0
    boxes = self % getVertexBox(vertexIdxs)
    do i = 1, 3
      if (.not. boxes(i) % ptr % hasFaces()) return

    end do
    
    ! Find common edge indices. Update edgeIdx only if common indices have been found.
    commonIdxs = boxes(1) % ptr % getFaceIdxs()
    do i = 2, 3
      commonIdxs = findCommon(commonIdxs, boxes(i) % ptr % getFaceIdxs())

    end do
    if (size(commonIdxs) > 0) faceIdx = commonIdxs(1)

  end function findCommonFaceIdx

  !! Function 'getAllCoordinates'
  !!
  !! Basic description:
  !!   Returns the 3-D coordinates of all the vertices in the shelf.
  !!
  !! Result:
  !!   allCoordinates -> Array listing the 3-D coordinates of all the vertices.
  !!
  function getAllCoordinates(self) result(allCoordinates)
    class(vertexShelf), intent(in)                :: self
    real(defReal), dimension(3, self % getSize()) :: allCoordinates
    integer(shortInt)                             :: i
    type(vertexBox)                               :: box

    do i = 1, self % getSize()
      box = self % getVertexBox(i)
      allCoordinates(:, i) = box % ptr % getCoordinates()

    end do

  end function getAllCoordinates
  
  !! Function 'getExtremalCoordinates'
  !!
  !! Basic description:
  !!   Returns the minimum and maximum x-, y- and z- coordinates of the vertices in the shelf.
  !!
  !! Result:
  !!   extremalCoordinates -> Array containing six entries: the first three list the minimum x-, y-
  !!   and z-coordinates, while the last three list the maximum x-, y- and z-coordinates.
  !!
  pure function getExtremalCoordinates(self) result(extremalCoordinates)
    class(vertexShelf), intent(in) :: self
    real(defReal), dimension(6)    :: extremalCoordinates
    
    extremalCoordinates = self % extremalCoordinates

  end function getExtremalCoordinates

  !!
  !!
  !!
  function getVertexBox_shortInt(self, idx) result(box)
    class(vertexShelf), intent(in) :: self
    integer(shortInt), intent(in)  :: idx
    type(vertexBox)                :: box
    type(topologicalObjectBox)     :: objectBox
    character(*), parameter        :: here = 'getVertexBox_shortInt (vertexShelf_class.f90)'

    ! First get a pointer to a polymorphic topological object from the shelf.
    objectBox = self % getObjectBox(idx)
    if (.not. associated(objectBox % ptr)) call fatalError(here, 'Invalid pointer for vertex with index: '//numToChar(idx)//'.')

    select type(ptr => objectBox % ptr)
      type is (vertex)
        box % ptr => ptr

      class default
        ! Should never happen.
        call fatalError(here, 'Object in vertexShelf with idx: '//numToChar(idx)//' is not a vertex.')

    end select

  end function getVertexBox_shortInt

  !!
  !!
  !!
  function getVertexBox_shortIntArray(self, idxs) result(boxes)
    class(vertexShelf), intent(in)                    :: self
    integer(shortInt), dimension(:), intent(in)       :: idxs
    type(vertexBox), dimension(size(idxs))            :: boxes
    type(topologicalObjectBox), dimension(size(idxs)) :: objectBoxes
    integer(shortInt)                                 :: i
    character(*), parameter                           :: here = 'getVertexBox_shortIntArray (vertexShelf_class.f90)'

    objectBoxes = self % getObjectBox(idxs)
    do i = 1, size(idxs)
      if (.not. associated(objectBoxes(i) % ptr)) &
      call fatalError(here, 'Invalid pointer for vertex with index: '//numToChar(idxs(i))//'.')

      select type(ptr => objectBoxes(i) % ptr)
        type is (vertex)
          boxes(i) % ptr => ptr

        class default
          ! Should never happen.
          call fatalError(here, 'Object in vertexShelf with idx: '//numToChar(idxs(i))//' is not a vertex.')

      end select

    end do

  end function getVertexBox_shortIntArray

  !! Function 'getVertexCoordinates_shortInt'
  !!
  !! Basic description:
  !!   Returns the 3-D coordinates of a vertex in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the vertex in the shelf.
  !!
  !! Result:
  !!   coords   -> 3-D coordinates of the vertex.
  !!
  function getVertexCoordinates_shortInt(self, idx) result(coords)
    class(vertexShelf), intent(in) :: self
    integer(shortInt), intent(in)  :: idx
    real(defReal), dimension(3)    :: coords
    type(vertexBox)                :: box

    box = self % getVertexBox(idx)
    coords = box % ptr % getCoordinates()

  end function getVertexCoordinates_shortInt

  !! Function 'getVertexCoordinates_shortInt'
  !!
  !! Basic description:
  !!   Returns the 3-D coordinates of vertices in the shelf.
  !!
  !! Arguments:
  !!   idxs [in] -> Indices of the vertices in the shelf.
  !!
  !! Result:
  !!   coords   -> 3-D coordinates of the vertices.
  !!
  function getVertexCoordinates_shortIntArray(self, idxs) result(coords)
    class(vertexShelf), intent(in)              :: self
    integer(shortInt), dimension(:), intent(in) :: idxs
    real(defReal), dimension(3, size(idxs))     :: coords
    type(vertexBox), dimension(size(idxs))      :: boxes
    integer(shortInt)                           :: i

    boxes = self % getVertexBox(idxs)
    do i = 1, size(idxs)
      coords(:, i) = boxes(i) % ptr % getCoordinates()

    end do

  end function getVertexCoordinates_shortIntArray

  !! Function 'getVertexEdgeIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of all the edges containing a vertex in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the vertex in the shelf.
  !!
  !! Result:
  !!   edgeIdxs -> Indices of all the edges containing the vertex.
  !!
  function getVertexEdgeIdxs(self, idx) result(edgeIdxs)
    class(vertexShelf), intent(in)               :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: edgeIdxs
    type(vertexBox)                              :: box

    box = self % getVertexBox(idx)
    edgeIdxs = box % ptr % getEdgeIdxs()

  end function getVertexEdgeIdxs

  !! Function 'getVertexElementIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of all the elements containing a vertex in the shelf.
  !!
  !! Arguments:
  !!   idx [in]    -> Index of the vertex in the shelf.
  !!
  !! Result:
  !!   elementIdxs -> Indices of all the elements containing the vertex.
  !!
  function getVertexElements(self, idx) result(elements)
    class(vertexShelf), intent(in)                        :: self
    integer(shortInt), intent(in)                         :: idx
    type(topologicalObjectBox), dimension(:), allocatable :: elements
    type(vertexBox)                                       :: box

    box = self % getVertexBox(idx)
    elements = box % ptr % getElements()

  end function getVertexElements

  !! Function 'getVertexFaceIdxs_shortInt'
  !!
  !! Basic description:
  !!   Returns the indices of all the faces containing a vertex in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the vertex in the shelf.
  !!
  !! Result:
  !!   faceIdxs -> Indices of all the faces containing the vertex.
  !!
  function getVertexFaceIdxs_shortInt(self, idx) result(faceIdxs)
    class(vertexShelf), intent(in)               :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: faceIdxs
    type(vertexBox)                              :: box

    box = self % getVertexBox(idx)
    faceIdxs = box % ptr % getFaceIdxs()

  end function getVertexFaceIdxs_shortInt

  !! Function 'getVertexFaceIdxs_shortIntArray'
  !!
  !! Basic description:
  !!   Returns the unique indices of all the faces containing a set of vertices
  !!   in the shelf.
  !!
  !! Arguments:
  !!   idxs [in] -> Indices of the vertices in the shelf.
  !!
  !! Result:
  !!   faceIdxs  -> Unique indices of all the faces containing the vertices.
  !!
  function getVertexFaceIdxs_shortIntArray(self, idxs) result(faceIdxs)
    class(vertexShelf), intent(in)               :: self
    integer(shortInt), dimension(:), intent(in)  :: idxs
    integer(shortInt), dimension(:), allocatable :: faceIdxs
    type(vertexBox), dimension(size(idxs))       :: boxes
    integer(shortInt)                            :: i

    boxes = self % getVertexBox(idxs)
    do i = 1, size(idxs)
      call append(faceIdxs, boxes(i) % ptr % getFaceIdxs(), .true.)

    end do

  end function getVertexFaceIdxs_shortIntArray

  !!
  !!
  !!
  subroutine init(self, coords)
    class(vertexShelf), intent(inout)          :: self
    real(defReal), dimension(:, :), intent(in) :: coords
    integer(shortInt)                          :: i, nVertices
    real(defReal), dimension(6)                :: extremalCoordinates
    type(vertexBox)                            :: box

    ! Allocate shelf then populate.
    nVertices = size(coords, 2)
    call self % allocateShelf(nVertices)
    do i = 1, nVertices
      call newVertexBox(i, coords(:, i), box)
      call self % addObject(box % ptr)

      ! Update extremal coordinates.
      if (i == 1) then
        extremalCoordinates = [coords(:, i), coords(:, i)]

      else
        extremalCoordinates(1:3) = min(extremalCoordinates(1:3), coords(:, i))
        extremalCoordinates(4:6) = max(extremalCoordinates(4:6), coords(:, i))

      end if

    end do

    ! Set extremal coordinates.
    self % extremalCoordinates = extremalCoordinates

  end subroutine init

  !! Subroutine 'initVertex'
  !!
  !! Basic description:
  !!   Initialises a vertex in the shelf.
  !!
  !! Arguments:
  !!   idx [in]    -> Index of the vertex.
  !!   coords [in] -> 3-D coordinates of the vertex.
  !!
  subroutine initVertex(self, idx, coords)
    class(vertexShelf), intent(inout)       :: self
    integer(shortInt), intent(in)           :: idx
    real(defReal), dimension(3), intent(in) :: coords
    type(vertexBox)                         :: box

    call newVertexBox(idx, coords, box)
    call self % addObject(box % ptr)

    ! Good practice.
    nullify(box % ptr)

  end subroutine initVertex
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  subroutine kill(self)
    class(vertexShelf), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % extremalCoordinates = ZERO

  end subroutine kill

  !! Subroutine 'setExtremalCoordinates'
  !!
  !! Basic description:
  !!   Sets the extremal coordinates of the vertices in the mesh.
  !!
  !! Arguments:
  !!   coords [in] -> defReal array of extremal coordinates (x_min, y_min, z_min, x_max, y_max and z_max).
  !!
  pure subroutine setExtremalCoordinates(self, coords)
    class(vertexShelf), intent(inout)       :: self
    real(defReal), dimension(6), intent(in) :: coords

    self % extremalCoordinates = coords

  end subroutine setExtremalCoordinates

end module vertexShelf_class