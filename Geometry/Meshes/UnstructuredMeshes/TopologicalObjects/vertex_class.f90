module vertex_class
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use genericProcedures,            only : append, areEqual, fatalError
  use numPrecision
  use publicObjects,                only : intersectionTestPayload, intersectionTestResult
  use topologicalObject_inter,      only : buildTopologicalObjectPayload, kill_super => kill, topologicalObject, &
                                           topologicalObjectBox
  use universalVariables
  use ratint
  
  implicit none
  private

  !!
  !!
  !!
  type, public, extends(buildTopologicalObjectPayload) :: buildVertexPayload
    real(defReal), dimension(3)  :: coordinates = ZERO
    type(ratint_t), dimension(3) :: ratintCoordinates
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
    type(ratint_t), dimension(3) :: ratintCoordinates
    type(topologicalObjectBox), dimension(:), allocatable :: sharingEdges, sharingElements, sharingFaces
  contains
    ! Build procedures.
    procedure :: addSharingEdge
    procedure :: addSharingElement
    procedure :: addSharingFace
    procedure :: build
    procedure :: kill
    ! Runtime procedures.
    procedure :: distanceSquared
    procedure :: getBoundingBoxBounds
    procedure :: getCentroid
    procedure :: getCoordinates
    procedure :: getRatintCoordinates
    procedure :: getSharingEdges
    procedure :: getSharingElements
    procedure :: getSharingFaceIdxs
    procedure :: getSharingFaces
    procedure :: intersects_BoundingBox
    procedure :: intersects_Ray
  end type vertex

contains
  !! Subroutine 'addEdgeIdx'
  !!
  !! Basic description:
  !!   Adds the index of an edge sharing the vertex.
  !!
  !! Arguments:
  !!   edgeIdx [in] -> Index of the edge.
  !!
  subroutine addSharingEdge(self, box)
    class(vertex), intent(inout)                          :: self
    type(topologicalObjectBox), intent(in)                :: box
    integer(shortInt)                                     :: nSharingEdges
    type(topologicalObjectBox), dimension(:), allocatable :: tempSharingEdges

    if (allocated(self % sharingEdges)) then
      nSharingEdges = size(self % sharingEdges)
      allocate(tempSharingEdges(nSharingEdges + 1))
      tempSharingEdges(1:nSharingEdges) = self % sharingEdges
      tempSharingEdges(nSharingEdges + 1) = box
      call move_alloc(tempSharingEdges, self % sharingEdges)

    else
      allocate(self % sharingEdges(1))
      self % sharingEdges(1) = box

    end if

  end subroutine addSharingEdge
  
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
  subroutine addSharingElement(self, box)
    class(vertex), intent(inout)                          :: self
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
  !!   Adds the index of a face sharing the vertex.
  !!
  !! Arguments:
  !!   faceIdx [in] -> Index of the face.
  !!
  subroutine addSharingFace(self, box)
    class(vertex), intent(inout)                          :: self
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
    self % ratintCoordinates = payloadPtr % ratintCoordinates

  end subroutine build

  !!
  !!
  !!
  pure function distanceSquared(self, r) result(dSquared)
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



  pure function getRatintCoordinates(self) result(ratintCoordinates)
    class(vertex), intent(in)   :: self
    type(ratint_t), dimension(3) :: ratintCoordinates
    
    ratintCoordinates = self % ratintCoordinates
  end function getRatintCoordinates

  !! Function 'getEdgeIdxs'
  !!
  !! Basic description:
  !!   Returns the edges containing the vertex.
  !!
  !! Result:
  !!   edgeIdxs -> Array listing the indices of the edges containing the vertex.
  !!
  function getSharingEdges(self) result(sharingEdges)
    class(vertex), intent(in)                             :: self
    type(topologicalObjectBox), dimension(:), allocatable :: sharingEdges

    if (allocated(self % sharingEdges)) then
      sharingEdges = self % sharingEdges

    else
      allocate(sharingEdges(0))

    end if

  end function getSharingEdges

  !! Function 'getVertexToElements'
  !!
  !! Basic description:
  !!   Returns the elements containing the vertex.
  !!
  !! Result:
  !!   elementIdxs -> Array listing the indices of the elements containing the vertex.
  !!
  function getSharingElements(self) result(sharingElements)
    class(vertex), target, intent(in)                     :: self
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
    class(vertex), intent(in)                    :: self
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
  
  !! Function 'getVertexToFaces'
  !!
  !! Basic description:
  !!   Returns the faces containing the vertex.
  !!
  !! Result:
  !!   faceIdxs -> Array listing the indices of the faces containing the vertex.
  !!
  function getSharingFaces(self) result(sharingFaces)
    class(vertex), intent(in)                             :: self
    type(topologicalObjectBox), dimension(:), allocatable :: sharingFaces
    
    if (allocated(self % sharingFaces)) then
      sharingFaces = self % sharingFaces

    else
      allocate(sharingFaces(0))

    end if

  end function getSharingFaces

  !!
  !!
  !!
  subroutine intersects_BoundingBox(self, boundingBox, doesIt)
    class(vertex), intent(in)                :: self
    type(axisAlignedBoundingBox), intent(in) :: boundingBox
    logical(defBool), intent(out)            :: doesIt
    character(*), parameter                  :: here = 'intersects_BoundingBox (vertex_class.f90)'

    call fatalError(here, 'Unsupported procedure.')

  end subroutine intersects_BoundingBox

  !!
  !!
  !!
  subroutine intersects_Ray(self, payload, result)
    class(vertex), intent(in)                    :: self
    class(intersectionTestPayload), intent(in)   :: payload
    class(intersectionTestResult), intent(inout) :: result
    real(defReal), dimension(3)                  :: vector
    real(defReal)                                :: dProjection

    ! Create vector from origin of the ray to the vertex.
    vector = self % coordinates - payload % r
    dProjection = dot_product(vector, payload % u)
    if (dProjection < ZERO .or. payload % dMax < dProjection) return

    if (areEqual(dot_product(vector, vector) - dProjection * dProjection, ZERO)) then
      result % intersects = .true.
      result % d = dProjection

    end if

  end subroutine intersects_Ray
  
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

    if (allocated(self % sharingEdges)) then
      do i = 1, size(self % sharingEdges)
        nullify(self % sharingEdges(i) % ptr)

      end do
      deallocate(self % sharingEdges)

    end if
    
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
  
end module vertex_class