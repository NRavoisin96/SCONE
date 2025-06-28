module objectKDTree_class
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use coord_class,                  only : coord
  use faceShelf_class,              only : faceShelf
  use numPrecision
  use objectNode_class,             only : objectNode
  use universalVariables,           only : INF
  use vertexShelf_class,            only : vertexShelf
  
  implicit none
  private
  
  !!
  !! k-dimensional (kd) tree. Data structure which is used to partition space and perform fast
  !! nearest-neighbour (NN) and mesh entry checks by recursively subdividing a cloud of points along
  !! alternating axis-aligned dimensions.
  !!
  !! Implementation adapted from Matthew B Kennel, University of California, San Diego (2004). 
  !! https://arxiv.org/pdf/physics/0408067.pdf.
  !!
  !! Private members:
  !!   nData    -> Number of data points in the tree.
  !!   data     -> 3-D coordinates of all the data points in the tree. Note that because Fortran
  !!               reads columns first it is preferable to arrange the data such that is it a
  !!               (3 x nData) array.
  !!   dataIdxs -> Internal sorting of the data point indices to enable cache-friendly searches.
  !!   root     -> Root node from which NN searches and mesh entry searches are initialised.
  !!
  type, public                                   :: objectKDTree
    private
    integer(shortInt)                            :: nObjects = 0, nLeaves = 0, nNodes = 0
    real(defReal), dimension(:, :), allocatable  :: data
    integer(shortInt), dimension(:), allocatable :: objectIdxs
    type(objectNode)                             :: root
  contains
    ! Build procedures.
    procedure          :: init
    procedure          :: kill
    ! Runtime procedures.
    procedure          :: findNearestObject
    generic            :: findPotentiallyIntersectedObjects => findPotentiallyIntersectedObjects_BoundingBox
    procedure, private :: findPotentiallyIntersectedObjects_BoundingBox
    procedure          :: getDataNumber
    procedure          :: getLeavesNumber
    procedure          :: getNodesNumber
    procedure          :: getRootBoundingBoxBounds
  end type objectKDTree

contains
  !! Subroutine 'init'
  !!
  !! Basic description:
  !!   Initialises a kd-tree from a supplied set of 3-D coordinates.
  !!
  !! Arguments:
  !!   data [in]                    -> Array of all 3-D coordinates.
  !!   boundingBoxes [in, optional] -> Array containing the bounding box corresponding to each coordinate.
  !!
  subroutine init(self, data, boundingBoxes)
    class(objectKDTree), intent(inout)                     :: self
    real(defReal), dimension(:, :), intent(in)             :: data
    type(axisAlignedBoundingBox), dimension(:), intent(in) :: boundingBoxes
    real(defReal), dimension(:, :), allocatable            :: tempData
    integer(shortInt)                                      :: i, nNodes, nData

    ! Compute the number of data points in the tree.
    nData = size(data, 2)
    self % nObjects = nData
    
    ! Allocate the 'dataIdxs' component of the tree and initialise it.
    allocate(self % objectIdxs(nData))
    do i = 1, nData
      self % objectIdxs(i) = i

    end do
    
    ! Build the tree's root node and all its children nodes.
    call self % root % init(data, self % objectIdxs, 1, nData, self % nNodes, self % nLeaves, boundingBoxes)
    
    ! Rearrange the tree's data for more cache-friendly searches later on.
    allocate(tempData(size(data, 1), nData))
    do i = 1, nData
      tempData(:, i) = data(:, self % objectIdxs(i))

    end do
    self % data = tempData

  end subroutine init
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(objectKDTree), intent(inout) :: self
    
    self % nObjects = 0
    self % nLeaves = 0
    self % nNodes = 0
    if (allocated(self % data)) deallocate(self % data)
    if (allocated(self % objectIdxs)) deallocate(self % objectIdxs)
    
    ! Kill the root node and all its children.
    call self % root % kill()

  end subroutine kill

  !!
  !!
  !!
  function findNearestObject(self, r, vertices, faces) result(idx)
    class(objectKDTree), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r
    type(vertexShelf), intent(in)           :: vertices
    type(faceShelf), intent(in)             :: faces
    integer(shortInt)                       :: idx
    real(defReal)                           :: radiusSquared

    radiusSquared = INF
    idx = 0
    call self % root % search(r, radiusSquared, self % objectIdxs, vertices, faces, idx)
    if (idx > 0) idx = self % objectIdxs(idx)
    
  end function findNearestObject

  !! Subroutine 'findIntersectedFaces'
  !!
  !! Basic description:
  !!   Returns the number and indices of all the faces in the k-d tree which are intersected by a given bounding box.
  !!
  subroutine findPotentiallyIntersectedObjects_BoundingBox(self, boundingBox, idxs)
    class(objectKDTree), intent(in)                           :: self
    type(axisAlignedBoundingBox), intent(in)                  :: boundingBox
    integer(shortInt), dimension(:), allocatable, intent(out) :: idxs

    ! Initialise nIntersectedPrimitives and intersectedPrimitiveIdxs, then descend the tree starting from the root node.
    allocate(idxs(0))
    call self % root % findPotentiallyIntersectedObjects(boundingBox, idxs)
    if (size(idxs) > 0) idxs = self % objectIdxs(idxs)

  end subroutine findPotentiallyIntersectedObjects_BoundingBox

  !! Function 'getDataNumber'
  !!
  !! Basic description:
  !!   Returns the number of data points in the tree.
  !!
  !! Result:
  !!   nData -> Number of data points in the tree.
  !!
  elemental function getDataNumber(self) result(nData)
    class(objectKDTree), intent(in) :: self
    integer(shortInt)               :: nData

    nData = self % nObjects

  end function getDataNumber

  !!
  !!
  !!
  elemental function getLeavesNumber(self) result(nLeaves)
    class(objectKDTree), intent(in) :: self
    integer(shortInt)               :: nLeaves

    nLeaves = self % nLeaves

  end function getLeavesNumber

  !! Function 'getNodesNumber'
  !!
  !! Basic description:
  !!   Returns the number of nodes in the tree.
  !!
  !! Result:
  !!   nNodes -> Number of nodes in the tree.
  !!
  elemental function getNodesNumber(self) result(nNodes)
    class(objectKDTree), intent(in) :: self
    integer(shortInt)               :: nNodes

    nNodes = self % nNodes

  end function getNodesNumber

  !! Function 'getRootBoundingBox'
  !!
  !! Basic description:
  !!   Returns the bounding box of the root node of the tree.
  !!
  !! Result:
  !!   boundingBox -> Bounding box of the root node of the tree.
  !!
  pure function getRootBoundingBoxBounds(self) result(boundingBoxBounds)
    class(objectKDTree), intent(in) :: self
    real(defReal), dimension(6)     :: boundingBoxBounds

    boundingBoxBounds = self % root % getBoundingBoxBounds()

  end function getRootBoundingBoxBounds

end module objectKDTree_class