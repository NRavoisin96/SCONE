module vertex_class
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use numPrecision
  use universalVariables
  use genericProcedures,            only : append, fatalError
  use topologicalObject_inter,      only : buildTopologicalObjectPayload, kill_super => kill, topologicalObject, &
                                           topologicalObjectBox
  
  implicit none
  private

  !!
  !!
  !!
  type, public, extends(buildTopologicalObjectPayload) :: buildVertexPayload
    real(defReal), dimension(3) :: coordinates
  end type buildVertexPayload

  !!
  !!
  !!
  type, public :: vertexBox
    type(vertex), pointer :: ptr => null()
  end type vertexBox
  
  !!
  !! Vertex of a given OpenFOAM mesh.
  !!
  !! Private members:
  !!   idx             -> Index of the vertex.
  !!   coordinates     -> 3-D coordinates of the vertex.
  !!   faceIdxs        -> Array of indices of the faces sharing the vertex.
  !!   elementIdxs     -> Array of indices of the elements sharing the vertex.
  !!   tetrahedronIdxs -> Array of indices of the tetrahedra sharing the vertex.
  !!   triangleIdxs    -> Array of indices of the triangles sharing the vertex.
  !!
  type, public, extends(topologicalObject)                :: vertex
    private
    real(defReal), dimension(3)                           :: coordinates = ZERO
    type(topologicalObjectBox), dimension(:), allocatable :: elements
    integer(shortInt), dimension(:), allocatable          :: faceIdxs, edgeIdxs
  contains
    ! Build procedures.
    procedure :: addFaceIdx
    procedure :: addEdgeIdx
    procedure :: addElement
    procedure :: build
    procedure :: kill
    ! Runtime procedures.
    procedure :: distanceSquared
    procedure :: getBoundingBoxBounds
    procedure :: getCentroid
    procedure :: getCoordinates
    procedure :: getEdgeIdxs
    procedure :: getElements
    procedure :: getFaceIdxs
    procedure :: hasEdges
    procedure :: hasFaces
    procedure :: intersects_BoundingBox
  end type vertex

contains
  
  !! Subroutine 'addFaceIdx'
  !!
  !! Basic description:
  !!   Adds the index of a face sharing the vertex.
  !!
  !! Arguments:
  !!   faceIdx [in] -> Index of the face.
  !!
  elemental subroutine addFaceIdx(self, faceIdx)
    class(vertex), intent(inout)  :: self
    integer(shortInt), intent(in) :: faceIdx
    
    call append(self % faceIdxs, faceIdx)

  end subroutine addFaceIdx

  !! Subroutine 'addEdgeIdx'
  !!
  !! Basic description:
  !!   Adds the index of an edge sharing the vertex.
  !!
  !! Arguments:
  !!   edgeIdx [in] -> Index of the edge.
  !!
  elemental subroutine addEdgeIdx(self, edgeIdx)
    class(vertex), intent(inout)  :: self
    integer(shortInt), intent(in) :: edgeIdx

    call append(self % edgeIdxs, edgeIdx)

  end subroutine addEdgeIdx
  
  !! Subroutine 'addElementIdx'
  !!
  !! Basic description:
  !!   Adds the index of an element sharing the vertex.
  !!
  !! Notes:
  !!   Due to the nature of the mesh importation process, here the index is only added if not
  !!   already present so as to avoid duplicates.
  !!
  !! Arguments:
  !!   elementIdx [in] -> Index of the element.
  !!
  subroutine addElement(self, box)
    class(vertex), intent(inout)                          :: self
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

  !!
  !!
  !!
  subroutine build(self, payload)
    class(vertex), intent(inout)                        :: self
    class(buildTopologicalObjectPayload), intent(inout) :: payload
    type(buildVertexPayload), pointer                   :: payloadPtr
    character(*), parameter                             :: here = 'build (vertex_class.f90)'

    ! Downcast payload to correct type.
    select type(ptr => payload)
      type is(buildVertexPayload)
        payloadPtr => ptr

      class default
        call fatalError(here, 'Invalid payload type.')

    end select

    self % coordinates = payloadPtr % coordinates

  end subroutine build

  !!
  !!
  !!
  function distanceSquared(self, r) result(dSquared)
    class(vertex), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: dSquared
    real(defReal), dimension(3)             :: diff

    diff = self % coordinates - r
    dSquared = dot_product(diff, diff)

  end function distanceSquared

  !!
  !!
  !!
  pure function getBoundingBoxBounds(self) result(bounds)
    class(vertex), intent(in)      :: self
    real(defReal), dimension(3, 2) :: bounds

    bounds = spread(self % coordinates, 2, 2)

  end function getBoundingBoxBounds

  !!
  !!
  !!
  pure function getCentroid(self) result(centroid)
    class(vertex), intent(in)   :: self
    real(defReal), dimension(3) :: centroid

    centroid = self % coordinates

  end function getCentroid
  
  !! Function 'getCoordinates'
  !!
  !! Basic description:
  !!   Returns the 3-D coordinates of the vertex.
  !!
  !! Result:
  !!   coordinates -> Array listing the x-, y- and z-coordinates of the vertex.
  !!
  pure function getCoordinates(self) result(coordinates)
    class(vertex), intent(in)   :: self
    real(defReal), dimension(3) :: coordinates
    
    coordinates = self % coordinates
  end function getCoordinates

  !! Function 'getEdgeIdxs'
  !!
  !! Basic description:
  !!   Returns the edges containing the vertex.
  !!
  !! Result:
  !!   edgeIdxs -> Array listing the indices of the edges containing the vertex.
  !!
  pure function getEdgeIdxs(self) result(edgeIdxs)
    class(vertex), intent(in)                           :: self
    integer(shortInt), dimension(size(self % edgeIdxs)) :: edgeIdxs

    edgeIdxs = self % edgeIdxs

  end function getEdgeIdxs

  !! Function 'getVertexToElements'
  !!
  !! Basic description:
  !!   Returns the elements containing the vertex.
  !!
  !! Result:
  !!   elementIdxs -> Array listing the indices of the elements containing the vertex.
  !!
  function getElements(self) result(elements)
    class(vertex), target, intent(in)                     :: self
    type(topologicalObjectBox), dimension(:), allocatable :: elements
    
    if (allocated(self % elements)) then
      elements = self % elements

    else
      allocate(elements(0))

    end if

  end function getElements
  
  !! Function 'getVertexToFaces'
  !!
  !! Basic description:
  !!   Returns the faces containing the vertex.
  !!
  !! Result:
  !!   faceIdxs -> Array listing the indices of the faces containing the vertex.
  !!
  pure function getFaceIdxs(self) result(faceIdxs)
    class(vertex), intent(in)                           :: self
    integer(shortInt), dimension(size(self % faceIdxs)) :: faceIdxs
    
    faceIdxs = self % faceIdxs

  end function getFaceIdxs

  !! Function 'hasEdges'
  !!
  !! Basic description:
  !!   Returns .true. if edgeIdxs is allocated.
  !!
  elemental function hasEdges(self) result(doesIt)
    class(vertex), intent(in) :: self
    logical(defBool)          :: doesIt

    doesIt = allocated(self % edgeIdxs)

  end function hasEdges

  !! Function 'hasTriangles'
  !!
  !! Basic description:
  !!   Returns .true. if triangleIdxs is allocated.
  !!
  elemental function hasFaces(self) result(doesIt)
    class(vertex), intent(in) :: self
    logical(defBool)          :: doesIt

    doesIt = allocated(self % faceIdxs)

  end function hasFaces

  !!
  !!
  !!
  elemental function intersects_BoundingBox(self, boundingBox) result(doesIt)
    class(vertex), intent(in)                :: self
    type(axisAlignedBoundingBox), intent(in) :: boundingBox
    logical(defBool)                         :: doesIt

    ! Return .true. if vertex is inside the bounds of the bounding box.
    doesIt = boundingBox % contains(self % coordinates)

  end function intersects_BoundingBox
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(vertex), intent(inout) :: self
    integer(shortInt)            :: i

    ! Superclass.
    call kill_super(self)
    
    ! Local.
    self % coordinates = ZERO
    if (allocated(self % faceIdxs)) deallocate(self % faceIdxs)
    if (allocated(self % edgeIdxs)) deallocate(self % edgeIdxs)
    
    if (allocated(self % elements)) then
      do i = 1, size(self % elements)
        nullify(self % elements(i) % ptr)

      end do
      deallocate(self % elements)

    end if
  
  end subroutine kill
  
end module vertex_class