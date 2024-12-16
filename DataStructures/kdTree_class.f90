module kdTree_class
  
  use coord_class,         only : coord
  use faceShelf_class,     only : faceShelf
  use numPrecision
  use universalVariables,  only : INF
  use node_class,          only : node
  use triangleShelf_class, only : triangleShelf
  use vertexShelf_class,   only : vertexShelf
  
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
  !!   nVertices    -> Number of vertices in the tree.
  !!   data         -> 3-D coordinates of all the vertices in the tree. Note that because Fortran
  !!                   reads columns first it is preferable to arrange the data such that is it a
  !!                   (3 x nVertices) array.
  !!   verticesIdxs -> Internal sorting of the vertices indices to enable cache-friendly searches.
  !!   root         -> Root node from which NN searches and mesh entry searches are initialised.
  !!
  type, public :: kdTree
    private
    integer(shortInt)                            :: nVertices = 0
    real(defReal), dimension(:, :), allocatable  :: data
    integer(shortInt), dimension(:), allocatable :: verticesIdxs
    type(node)                                   :: root
  contains
    ! Build procedures.
    procedure :: init
    procedure :: kill
    ! Runtime procedures.
    procedure :: findNearestVertex
    procedure :: findIntersectedFace
    procedure :: findIntersectedTriangle
  end type kdTree

contains
  !! Subroutine 'init'
  !!
  !! Basic description:
  !!   Initialises a kd-tree from a supplied set of 3-D coordinates.
  !!
  !! Arguments:
  !!   data          -> Array of all 3-D coordinates.
  !!   transposeData -> Does the data need to be transposed?
  !!
  subroutine init(self, data, transposeData)
    class(kdTree), intent(inout)                :: self
    real(defReal), dimension(:, :), intent(in)  :: data
    logical(defBool), intent(in)                :: transposeData
    real(defReal), dimension(:, :), allocatable :: tempData
    integer(shortInt)                           :: i, nVertices
    
    ! First check if the data needs to be transposed.
    if (transposeData) then
      self % data = transpose(data)
    else
      self % data = data

    end if

    ! Compute the number of vertices in the tree.
    nVertices = size(self % data, 2)
    ! Allocate the 'verticesIdxs' component of the tree and initialise it.
    allocate(self % verticesIdxs(nVertices))
    do i = 1, nVertices
      self % verticesIdxs(i) = i

    end do
    
    ! Build the tree's root node and all its children nodes.
    call self % root % init(self % data, self % verticesIdxs, 1, nVertices)
    
    ! Rearrange the tree's data for more cache-friendly searches later on.
    allocate(tempData(3, nVertices))
    do i = 1, nVertices
      tempData(:, i) = self % data(:, self % verticesIdxs(i))

    end do
    deallocate(self % data)
    self % data = tempData

  end subroutine init
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(kdTree), intent(inout) :: self
    
    self % nVertices = 0
    if (allocated(self % data)) deallocate(self % data)
    if (allocated(self % verticesIdxs)) deallocate(self % verticesIdxs)
    
    ! Kill the root node and all its children.
    call self % root % kill()

  end subroutine kill
  
  !! Function 'findNearestVertex'
  !!
  !! Basic description:
  !!   Returns the index of the vertex in the tree whose coordinates are closest to the supplied 3-D
  !!   coordinates.
  !!
  !! Result:
  !!   vertexIdx -> Index of the vertex in the tree whose coordinates are closest to the supplied 
  !!                3-D coordinates.
  !!
  pure function findNearestVertex(self, r) result(vertexIdx)
    class(kdTree), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    integer(shortInt)                       :: idx, vertexIdx
    real(defReal)                           :: searchRadius
    
    ! Initialise searchRadius = INF and begin search from root node.
    searchRadius = INF
    call self % root % search(self % data, r, searchRadius, idx)
    
    ! Retrieve the actual vertex index from the 'verticesIdxs' component of the tree.
    vertexIdx = self % verticesIdxs(idx)

  end function findNearestVertex

  !! Subroutine 'findIntersectedFace'
  !!
  !! Basic description:
  !!   Returns the index of the closest face intersected by a line segment.
  !!
  !! Arguments:
  !!   vertices [in]     -> A vertexShelf.
  !!   faces [in]        -> A faceShelf.
  !!   d [inout]         -> Distance from the line segment's origin to the point of intersection.
  !!   coords [inout]    -> Particle's coordinates.
  !!   edgeIdx [inout]   -> Index of the edge intersected by the particle.
  !!   vertexIdx [inout] -> Index of the vertex intersected by the particle.
  !!
  elemental subroutine findIntersectedFace(self, vertices, faces, d, coords, edgeIdx, vertexIdx)
    class(kdTree), intent(in)        :: self
    type(vertexShelf), intent(in)    :: vertices
    type(faceShelf), intent(in)      :: faces
    real(defReal), intent(inout)     :: d
    type(coord), intent(inout)       :: coords
    integer(shortInt), intent(inout) :: edgeIdx, vertexIdx
    
    ! Start searching the root node for potential face intersections.
    call self % root % findIntersectedFace(vertices, faces, coords % r, self % verticesIdxs, d, &
                                           coords, edgeIdx, vertexIdx)

  end subroutine findIntersectedFace
  
  !! Subroutine 'findIntersectedTriangle'
  !!
  !! Basic description:
  !!   Returns the index of the closest triangle intersected by a line segment.
  !!
  !! Arguments:
  !!   vertices [in]  -> A vertexShelf.
  !!   triangles [in] -> A triangleShelf.
  !!   d [inout]      -> Distance from the line segment's origin to the point of intersection.
  !!   coords [inout] -> Particle's coordinates.
  !!   edgeIdx [inout]   -> Index of the edge intersected by the particle.
  !!   vertexIdx [inout] -> Index of the vertex intersected by the particle.
  !!
  elemental subroutine findIntersectedTriangle(self, vertices, triangles, d, coords, edgeIdx, vertexIdx)
    class(kdTree), intent(in)        :: self
    type(vertexShelf), intent(in)    :: vertices
    type(triangleShelf), intent(in)  :: triangles
    real(defReal), intent(inout)     :: d
    type(coord), intent(inout)       :: coords
    integer(shortInt), intent(inout) :: edgeIdx, vertexIdx
    
    ! Start searching the root node for potential triangle intersections.
    call self % root % findIntersectedTriangle(vertices, triangles, coords % r, self % verticesIdxs, d, &
                                               coords, edgeIdx, vertexIdx)

  end subroutine findIntersectedTriangle

end module kdTree_class