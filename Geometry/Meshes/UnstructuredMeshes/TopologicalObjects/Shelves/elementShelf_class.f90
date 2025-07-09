module elementShelf_class
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use coord_class,                  only : coord
  use edge_class,                   only : edgeBox
  use element_class,                only : buildElementInfo, element, elementBox, inclusionTestResult
  use elementFactory_func,          only : newElementBox
  use face_class,                   only : faceBox, orientatedFaceBox
  use genericProcedures,            only : fatalError, numToChar
  use numPrecision
  use topologicalObject_inter,      only : topologicalObjectBox
  use topologicalObjectShelf_inter, only : topologicalObjectShelf
  use vertex_class,                 only : vertexBox
  
  implicit none
  private
  
  !!
  !! Storage space for elements in a given OpenFOAM mesh.
  !!
  !! Private members:
  !!   shelf -> Array to store elements.
  !!
  type, public, extends(topologicalObjectShelf) :: elementShelf
    private
  contains
    procedure                                   :: addEdgeToElement
    generic                                     :: addElement => addElement_info, addElement_infoArray
    procedure, private                          :: addElement_info
    procedure, private                          :: addElement_infoArray
    procedure                                   :: addFaceToElement
    procedure                                   :: addVertexToElement
    procedure                                   :: computeFaceIntersection
    procedure                                   :: computePotentialFaceIdxs
    procedure                                   :: getElementBoundingBox
    generic                                     :: getElementBox => getElementBox_shortInt, getElementBox_shortIntArray
    procedure, private                          :: getElementBox_shortInt
    procedure, private                          :: getElementBox_shortIntArray
    procedure                                   :: getElementCentroid
    procedure                                   :: getElementEdges
    procedure                                   :: getElementIsActive
    procedure                                   :: getElementIsConvex
    procedure                                   :: getElementLocalId
    procedure                                   :: getElementOrientatedFaces
    procedure                                   :: getElementParentIdx
    procedure                                   :: getElementType
    procedure                                   :: getElementVertices
    procedure                                   :: getElementVolume
    procedure                                   :: init
    procedure                                   :: pushFromElementBoundary
    generic                                     :: setElementLocalId => setElementLocalId_shortInt, setElementLocalId_shortIntArray
    procedure, private                          :: setElementLocalId_shortInt
    procedure, private                          :: setElementLocalId_shortIntArray
    procedure                                   :: isPointInside
  end type elementShelf

contains
  !! Subroutine 'addEdgeIdxToElement'
  !!
  !! Basic description:
  !!   Adds the index of an edge to an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the element in the shelf.
  !!   edgeIdx [in] -> Index of the edge in the element.
  !!
  subroutine addEdgeToElement(self, idx, edge)
    class(elementShelf), intent(inout) :: self
    integer(shortInt), intent(in)      :: idx
    type(edgeBox), intent(in)          :: edge
    type(elementBox)                   :: box

    box = self % getElementBox(idx)
    call box % ptr % addEdge(edge)

  end subroutine addEdgeToElement

  !! Subroutine 'initElement'
  !!
  !! Basic description:
  !!   Initialises an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in]      -> Index of the element in the shelf.
  !!   faces [in]    -> A faceShelf.
  !!   vertices [in] -> A vertexShelf.
  !!
  subroutine addElement_info(self, info)
    class(elementShelf), intent(inout) :: self
    type(buildElementInfo), intent(in) :: info
    type(elementBox)                   :: box

    call newElementBox(info, box)
    call self % addObject(box % ptr)

  end subroutine addElement_info

  !!
  !!
  !!
  subroutine addElement_infoArray(self, infos)
    class(elementShelf), intent(inout)               :: self
    type(buildElementInfo), dimension(:), intent(in) :: infos
    integer(shortInt)                                :: i
    type(elementBox)                                 :: box

    do i = 1, size(infos)
      call newElementBox(infos(i), box)
      call self % addObject(box % ptr)

    end do

  end subroutine addElement_infoArray

  !! Subroutine 'addFaceIdxToElement'
  !!
  !! Basic description:
  !!   Adds the index of a face to an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the element in the shelf.
  !!   faceIdx [in] -> Index of the face in the element.
  !!
  subroutine addFaceToElement(self, idx, orientatedFace)
    class(elementShelf), intent(inout)  :: self
    integer(shortInt), intent(in)       :: idx
    type(orientatedFaceBox), intent(in) :: orientatedFace
    type(elementBox)                    :: box

    box = self % getElementBox(idx)
    call box % ptr % addFace(orientatedFace)

  end subroutine addFaceToElement

  !! Subroutine 'addVertexIdxToElement'
  !!
  !! Basic description:
  !!   Adds the index of a vertex to an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in]       -> Index of the element in the shelf.
  !!   vertexIdx [in] -> Index of the vertex in the element.
  !!
  subroutine addVertexToElement(self, idx, vertex)
    class(elementShelf), intent(inout) :: self
    integer(shortInt), intent(in)      :: idx
    type(vertexBox), intent(in)        :: vertex
    type(elementBox)                   :: box

    box = self % getElementBox(idx)
    call box % ptr % addVertex(vertex)

  end subroutine addVertexToElement

  !! Subroutine 'computeFaceIntersection'
  !!
  !! Basic description:
  !!   Computes the intersection of a line segment with the faces of an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in]                 -> Index of the element in the shelf.
  !!   r [in]                   -> Line segment's origin coordinates.
  !!   rEnd [in]                -> Line segment's end coordinates.
  !!   potentialFaceIdxs [in]   -> Indices of the element faces potentially intersected by the line segment.
  !!   faces [in]               -> A faceShelf.
  !!   intersectedFaceIdx [out] -> Index of the intersected element face.
  !!   lambda [out]             -> Fraction of the line segment to intersection.
  !!
  subroutine computeFaceIntersection(self, idx, r, rEnd, potentialFaceIdxs, intersectedFaceIdx, lambda)
    class(elementShelf), intent(in)             :: self
    integer(shortInt), intent(in)               :: idx
    real(defReal), dimension(3), intent(in)     :: r, rEnd
    integer(shortInt), dimension(:), intent(in) :: potentialFaceIdxs
    integer(shortInt), intent(out)              :: intersectedFaceIdx
    real(defReal), intent(out)                  :: lambda
    type(elementBox)                            :: box

    box = self % getElementBox(idx)
    call box % ptr % computeIntersectedFace(r, rEnd, potentialFaceIdxs, intersectedFaceIdx, lambda)

  end subroutine computeFaceIntersection

  !! Function 'computePotentialFaceIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the potentially intersected faces of an element in the shelf by a line segment.
  !!
  !! Argument:
  !!   idx [in]          -> Index of the element in the shelf.
  !!   rEnd [in]         -> Line segment's end coordinates.
  !!   faces [in]        -> A faceShelf.
  !!
  !! Result:
  !!   potentialFaceIdxs -> Indices of the element's faces potentially intersected by the line segment.
  !!
  function computePotentialFaceIdxs(self, idx, rEnd) result(potentialFaceIdxs)
    class(elementShelf), intent(in)              :: self
    integer(shortInt), intent(in)                :: idx
    real(defReal), dimension(3), intent(in)      :: rEnd
    integer(shortInt), dimension(:), allocatable :: potentialFaceIdxs
    type(elementBox)                             :: box

    box = self % getElementBox(idx)
    potentialFaceIdxs = box % ptr % computePotentialFaces(rEnd)

  end function computePotentialFaceIdxs

  !! Function 'getElementBoundingBox'
  !!
  !! Basic description:
  !!   Returns the bounding box of an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in]    -> Index of the element in the shelf.
  !!
  !! Result:
  !!   boundingBox -> Array containing the bounding box of the element.
  !!
  function getElementBoundingBox(self, idx) result(boundingBox)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    type(axisAlignedBoundingBox)    :: boundingBox
    type(elementBox)                :: box

    box = self % getElementBox(idx)
    boundingBox = box % ptr % getBoundingBox()

  end function getElementBoundingBox

  !!
  !!
  !!
  function getElementBox_shortInt(self, idx) result(box)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    type(elementBox)                :: box
    type(topologicalObjectBox)      :: objectBox
    character(*), parameter         :: here = 'getElementBox_shortInt (elementShelf_class.f90)'

    ! First get a pointer to a polymorphic topological object from the shelf.
    objectBox = self % getObjectBox(idx)
    if (.not. associated(objectBox % ptr)) call fatalError(here, 'Invalid pointer for element with index: '//numToChar(idx)//'.')

    select type(ptr => objectBox % ptr)
      type is (element)
        box % ptr => ptr

      class default
        ! Should never happen.
        call fatalError(here, 'Object in elementShelf with idx: '//numToChar(idx)//' is not an element.')

    end select

  end function getElementBox_shortInt

  !!
  !!
  !!
  function getElementBox_shortIntArray(self, idxs) result(boxes)
    class(elementShelf), intent(in)                   :: self
    integer(shortInt), dimension(:), intent(in)       :: idxs
    type(elementBox), dimension(size(idxs))           :: boxes
    type(topologicalObjectBox), dimension(size(idxs)) :: objectBoxes
    integer(shortInt)                                 :: i
    character(*), parameter                           :: here = 'getElementBox_shortIntArray (elementShelf_class.f90)'

    objectBoxes = self % getObjectBox(idxs)
    do i = 1, size(idxs)
      if (.not. associated(objectBoxes(i) % ptr)) &
      call fatalError(here, 'Invalid pointer for element with index: '//numToChar(idxs(i))//'.')

      select type(ptr => objectBoxes(i) % ptr)
        type is (element)
          boxes(i) % ptr => ptr

        class default
          ! Should never happen.
          call fatalError(here, 'Object in elementShelf with idx: '//numToChar(idxs(i))//' is not an element.')

      end select

    end do

  end function getElementBox_shortIntArray

  !! Function 'getElementCentroid'
  !!
  !! Basic description:
  !!   Returns the centroid of an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element in the shelf.
  !!
  !! Result:
  !!   centroid -> 3-D coordinates of the element's centroid.
  !!
  function getElementCentroid(self, idx) result(centroid)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    real(defReal), dimension(3)     :: centroid
    type(elementBox)                :: box

    box = self % getElementBox(idx)
    centroid = box % ptr % getCentroid()

  end function getElementCentroid

  !! Function 'getElementEdgeIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the edges in an element of the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element in the shelf.
  !!
  !! Result:
  !!   edgeIdxs -> Indices of the edges in the element.
  !!
  function getElementEdges(self, idx) result(edges)
    class(elementShelf), intent(in)          :: self
    integer(shortInt), intent(in)            :: idx
    type(edgeBox), dimension(:), allocatable :: edges
    type(elementBox)                         :: box

    box = self % getElementBox(idx)
    edges = box % ptr % getEdges()

  end function getElementEdges

  !!
  !!
  !!
  function getElementIsActive(self, idx) result(isActive)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    logical(defBool)                :: isActive
    type(elementBox)                :: box

    box = self % getElementBox(idx)
    isActive = box % ptr % getIsActive()

  end function getElementIsActive

  !! Function 'getElementIsConvex'
  !!
  !! Basic description:
  !!   Returns .true. if an element in the shelf is convex.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element in the shelf.
  !!
  !! Result:
  !!   isConvex -> .true. if the element is convex.
  !!
  function getElementIsConvex(self, idx) result(isConvex)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    logical(defBool)                :: isConvex
    type(elementBox)                :: box

    box = self % getElementBox(idx)
    isConvex = box % ptr % getIsConvex()

  end function getElementIsConvex

  !!
  !!
  !!
  function getElementLocalId(self, idx) result(localId)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    integer(shortInt)               :: localId
    type(elementBox)                :: box

    box = self % getElementBox(idx)
    localId = box % ptr % getLocalId()

  end function getElementLocalId

  !! Function 'getElementFaceIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the faces in an element of the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element in the shelf.
  !!
  !! Result:
  !!   faceIdxs -> Indices of the faces in the element.
  !!
  function getElementOrientatedFaces(self, idx) result(orientatedFaces)
    class(elementShelf), intent(in)                    :: self
    integer(shortInt), intent(in)                      :: idx
    type(orientatedFaceBox), dimension(:), allocatable :: orientatedFaces
    type(elementBox)                                   :: box

    box = self % getElementBox(idx)
    orientatedFaces = box % ptr % getOrientatedFaces()

  end function getElementOrientatedFaces

  !! Function 'getElementParentIdx'
  !!
  !! Basic description:
  !!   Returns the index of the parent element of an element of the shelf.
  !!
  !! Arguments:
  !!   idx [in]  -> Index of the element in the shelf.
  !!
  !! Result:
  !!   parentIdx -> Index of the parent element of the element
  !!
  function getElementParentIdx(self, idx) result(parentIdx)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    integer(shortInt)               :: parentIdx
    type(elementBox)                :: box

    box = self % getElementBox(idx)
    parentIdx = box % ptr % getParentIdx()

  end function getElementParentIdx

  !!
  !!
  !!
  function getElementType(self, idx) result(type)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    character(:), allocatable       :: type
    type(elementBox)                :: box

    box = self % getElementBox(idx)
    type = box % ptr % getType()

  end function getElementType

  !! Function 'getElementVertexIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in an element of the shelf.
  !!
  !! Arguments:
  !!   idx [in]   -> Index of the element in the shelf.
  !!
  !! Result:
  !!   vertexIdxs -> Indices of the vertices in the element.
  !!
  function getElementVertices(self, idx) result(vertices)
    class(elementShelf), intent(in)            :: self
    integer(shortInt), intent(in)              :: idx
    type(vertexBox), dimension(:), allocatable :: vertices
    type(elementBox)                           :: box

    box = self % getElementBox(idx)
    vertices = box % ptr % getVertices()

  end function getElementVertices

  !! Function 'getElementVolume'
  !!
  !! Basic description:
  !!   Returns the volume of an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element in the shelf.
  !!
  !! Result:
  !!   volume   -> Volume of the element.
  !!
  function getElementVolume(self, idx) result(volume)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    real(defReal)                   :: volume
    type(elementBox)                :: box

    box = self % getElementBox(idx)
    volume = box % ptr % getVolume()

  end function getElementVolume

  !!
  !!
  !!
  subroutine init(self, infos)
    class(elementShelf), intent(inout)               :: self
    type(buildElementInfo), dimension(:), intent(in) :: infos
    integer(shortInt)                                :: i, nElements
    type(elementBox)                                 :: box

    ! Allocate shelf then populate it.
    nElements = size(infos)
    call self % allocateShelf(nElements)
    do i = 1, nElements
      call newElementBox(infos(i), box)
      call self % addObject(box % ptr)

    end do

  end subroutine init

  !!
  !!
  !!
  subroutine pushFromElementBoundary(self, idx, coords)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    type(coord), intent(inout)      :: coords
    type(elementBox)                :: box

    box = self % getElementBox(idx)
    call box % ptr % pushFromBoundary(coords)

  end subroutine pushFromElementBoundary

  !!
  !!
  !!
  subroutine setElementLocalId_shortInt(self, idx, localId)
    class(elementShelf), intent(inout) :: self
    integer(shortInt), intent(in)      :: idx, localId
    type(elementBox)                   :: box

    box = self % getElementBox(idx)
    call box % ptr % setLocalId(localId)

  end subroutine setElementLocalId_shortInt

  !!
  !!
  !!
  subroutine setElementLocalId_shortIntArray(self, idxs, localId)
    class(elementShelf), intent(inout)          :: self
    integer(shortInt), dimension(:), intent(in) :: idxs
    integer(shortInt), intent(in)               :: localId
    type(elementBox), dimension(size(idxs))     :: boxes
    integer(shortInt)                           :: i

    boxes = self % getElementBox(idxs)
    do i = 1, size(idxs)
      call boxes(i) % ptr % setLocalId(localId)

    end do

  end subroutine setElementLocalId_shortIntArray

  !! Subroutine 'testForInclusion'
  !!
  !! Basic description:
  !!   Tests whether a given element in the shelf contains a point.
  !!
  !! Arguments:
  !!   idx [in]              -> Index of the element in the shelf.
  !!   r [in]                -> 3-D coordinates of the point.
  !!   faces [in]            -> A faceShelf.
  !!   failedFaceIdx [out]   -> Index of the first element's face for which the test fails.
  !!   surfTolFaceIdxs [out] -> Indices of the element's faces on which the point lies.
  !!
  function isPointInside(self, idx, r) result(result)
    class(elementShelf), intent(in)         :: self
    integer(shortInt), intent(in)           :: idx
    real(defReal), dimension(3), intent(in) :: r
    type(inclusionTestResult)               :: result
    type(elementBox)                        :: box

    box = self % getElementBox(idx)
    result = box % ptr % isPointInside(r)

  end function isPointInside
  
end module elementShelf_class