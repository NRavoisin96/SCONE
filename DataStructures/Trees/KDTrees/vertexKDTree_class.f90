module vertexKDTree_class

  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use numPrecision
  use universalVariables,           only : INF
  use vertexNode_class,             only : vertexNode

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
  type, public                                   :: vertexKDTree
    private
    integer(shortInt)                            :: nData = 0, nLeaves = 0, nNodes = 0
    real(defReal), dimension(:, :), allocatable  :: data
    integer(shortInt), dimension(:), allocatable :: dataIdxs
    type(vertexNode)                             :: root
  contains
    ! Build procedures.
    procedure :: init
    procedure :: kill
    ! Runtime procedures.
    procedure :: findNearestVertex
    procedure :: getDataNumber
    procedure :: getLeavesNumber
    procedure :: getNodesNumber
    procedure :: getRootBoundingBox
  end type vertexKDTree

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
  subroutine init(self, data)
    class(vertexKDTree), intent(inout)                     :: self
    real(defReal), dimension(:, :), intent(in)             :: data
    real(defReal), dimension(:, :), allocatable            :: tempData
    integer(shortInt)                                      :: i, nLeaves, nNodes, nData

    ! Compute the number of data points in the tree.
    nData = size(data, 2)
    self % nData = nData
    
    ! Allocate the 'dataIdxs' component of the tree and initialise it.
    allocate(self % dataIdxs(nData))
    do i = 1, nData
      self % dataIdxs(i) = i

    end do
    
    ! Build the tree's root node and all its children nodes.
    call self % root % init(data, self % dataIdxs, 1, nData, self % nNodes, self % nLeaves)
    
    ! Rearrange the tree's data for more cache-friendly searches later on.
    allocate(tempData(size(data, 1), nData))
    do i = 1, nData
      tempData(:, i) = data(:, self % dataIdxs(i))

    end do
    self % data = tempData

  end subroutine init
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(vertexKDTree), intent(inout) :: self
    
    self % nData = 0
    self % nLeaves = 0
    self % nNodes = 0
    if (allocated(self % data)) deallocate(self % data)
    if (allocated(self % dataIdxs)) deallocate(self % dataIdxs)
    
    ! Kill the root node and all its children.
    call self % root % kill()

  end subroutine kill

  !!
  !!
  !!
  function findNearestVertex(self, r) result(vertexIdx)
    class(vertexKDTree), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r
    integer(shortInt)                       :: idx, vertexIdx
    real(defReal)                           :: radiusSquared

    ! Initialise searchRadius = INF and begin search from root node.
    radiusSquared = INF
    idx = 0
    call self % root % search(self % data, r, radiusSquared, idx)
    
    ! Retrieve the actual vertex index from the 'verticesIdxs' component of the tree.
    if (idx > 0) vertexIdx = self % dataIdxs(idx)
    
  end function findNearestVertex

  !! Function 'getDataNumber'
  !!
  !! Basic description:
  !!   Returns the number of data points in the tree.
  !!
  !! Result:
  !!   nData -> Number of data points in the tree.
  !!
  elemental function getDataNumber(self) result(nData)
    class(vertexKDTree), intent(in) :: self
    integer(shortInt)               :: nData

    nData = self % nData

  end function getDataNumber

  !!
  !!
  !!
  elemental function getLeavesNumber(self) result(nLeaves)
    class(vertexKDTree), intent(in) :: self
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
    class(vertexKDTree), intent(in) :: self
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
  pure function getRootBoundingBox(self) result(boundingBox)
    class(vertexKDTree), intent(in) :: self
    type(axisAlignedBoundingBox)    :: boundingBox

    boundingBox = self % root % getBoundingBox()

  end function getRootBoundingBox

end module vertexKDTree_class