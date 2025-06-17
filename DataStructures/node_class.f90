module node_class
  
  use coord_class,         only : coord
  use faceShelf_class,     only : faceShelf
  use genericProcedures,   only : append, areEqual, quickSort, removeDuplicates, swap
  use numPrecision
  use universalVariables,  only : ZERO, HALF, SURF_TOL, INF
  use vertexShelf_class,   only : vertexShelf
  
  implicit none
  private
  
  !!
  !! Node of a kd-tree. Used to recursively subdivide space into smaller and smaller numbers of
  !! vertices. Each node subdivides space along an axis-aligned dimension known as the cut
  !! dimension, which is taken as the dimension of greatest extent among the vertices in the node.
  !! The node is then split into (generally) two children nodes along the cut dimension. The value
  !! used to split the node along the cut dimension is taken as the average of the vertices'
  !! coordinates along the cut dimension. The child node containing vertices whose coordinates along
  !! the cut dimension less than the cut value is termed the left node and vice versa.
  !!
  !! Private members:
  !!   cutDimension   -> Index of the dimension (1 = x, 2 = y, 3 = z) used to split the node into
  !!                     children nodes.
  !!   lowerBound     -> Lowest-index vertex in the node.
  !!   upperBound     -> Greatest-index vertex in the node.
  !!   bucketSize     -> Maximum number of vertices in a terminal node.
  !!   hasLeft        -> Does the node have a left child node?
  !!   hasRight       -> Does the node have a right child node?
  !!   cutValue       -> Value along the cutDimension used to split the node into children nodes.
  !!   cutValue_left  -> CutValue of the node's left child node.
  !!   cutValue_right -> CutValue of the node's right child node.
  !!   left           -> Left child node.
  !!   right          -> Right child node.
  !!   boundingBox    -> Axis-aligned bounding box (AABB) of the node.
  !!
  type, public                  :: node
    private
    integer(shortInt)           :: cutDimension = 0, lowerBound = 0, upperBound = 0, bucketSize = 20, &
                                   idx = 0, parentIdx = 0
    logical(defBool)            :: hasLeft = .false., hasRight = .false.
    real(defReal)               :: cutValue = ZERO, cutValue_left = ZERO, cutValue_right = ZERO
    type(node), allocatable     :: left, right
    real(defReal), dimension(6) :: boundingBox = [INF, INF, INF, -INF, -INF, -INF]
  contains
    ! Build procedures.
    procedure                   :: computeBoundingBox
    procedure                   :: computeCutVertexIdx
    procedure                   :: init
    procedure                   :: kill
    ! Runtime procedures.
    procedure                   :: findIntersectedFace
    procedure                   :: findPotentialElementIdxs
    procedure                   :: getBoundingBox
    procedure                   :: getVertices
    procedure                   :: process
    procedure                   :: search
  end type node
contains
  
  !! Subroutine 'computeBoundingBox'
  !!
  !! Basic description:
  !!   Computes the bounding box of the node along a specified dimension.
  !!
  !! Detailed description:
  !!   First sorts coordinates in the node in increasing order along the specified dimension. The
  !!   bounding box of the node along this dimension is then simply given by the minimum and maximum
  !!   coordinates along the dimension.
  !!
  !! Arguments:
  !!   dimension [in] -> Dimension along which to compute the bounding box.
  !!   coords [in]    -> Coordinates of the vertices in the node along the cut dimension.
  !!
  pure subroutine computeBoundingBox(self, dimension, nVertices, coords)
    class(node), intent(inout)                 :: self
    integer(shortInt), intent(in)              :: dimension, nVertices
    real(defReal), dimension(:, :), intent(in) :: coords
    integer(shortInt)                          :: i, j
    real(defReal), dimension(:), allocatable   :: dimensionCoords
    real(defReal)                              :: lowerCoord, upperCoord, maxCoord, minCoord
    
    ! Loop through all supplied coordinate dimensions.
    do i = 1, size(coords, 1)
      ! Retrieve the minimum coordinate along the supplied dimension and initialise the maximum
      ! and minimum coordinates.
      dimensionCoords = coords(i, :)
      lowerCoord = dimensionCoords(1)
      upperCoord = dimensionCoords(nVertices)
      minCoord = min(lowerCoord, upperCoord)
      maxCoord = max(lowerCoord, upperCoord)
      
      ! Sweep through all vertices in the input data between lowerBound and upperBound.
      do j = 3, nVertices, 2
        ! Retrieve the coordinate of the next two adjacent vertices along the supplied dimension.
        lowerCoord = dimensionCoords(j - 1)
        upperCoord = dimensionCoords(j)

        ! Update minCoord and maxCoord.
        minCoord = min(minCoord, lowerCoord, upperCoord)
        maxCoord = max(maxCoord, lowerCoord, upperCoord)

      end do
      
      ! Set the lower and upper bound for the interval.
      self % boundingBox(dimension) = min(minCoord, self % boundingBox(dimension))
      self % boundingBox(dimension + 3) = max(maxCoord, self % boundingBox(dimension + 3))

    end do

  end subroutine computeBoundingBox

  !! Function 'computeCutVertexIdx'
  !!
  !! Basic description:
  !!   Returns the index of the vertex used to subdivide the node into children nodes along a given
  !!   cut dimension.
  !!
  !! Detailed description:
  !!   First performs a sweep through all supplied vertices and sorts the entries in the 
  !!   verticesIdxs array depending on whether the coordinate of a given vertex along the cut 
  !!   dimension is less than or equal to the average value of all the vertices' coordinates along 
  !!   this dimension. The reason this is done is to later re-organise the order of the data in the 
  !!   kd-tree to enable cache-friendly search optimisation. After the vertices' indices have been 
  !!   sorted, the cut vertex is simply chosen as that whose coordinate along the cut dimension is 
  !!   closest to the average value, while still being inferior to it.
  !!
  !! Arguments:
  !!   coordinates [in]     -> Coordinates of the vertices in the node along the cut dimension.
  !!   verticesIdxs [inout] -> Internal re-ordering of the vertices indices in the tree.
  !!   average [in]         -> Average value of the coordinates along the cut dimension.
  !!   lowerBound [in]      -> Vertex of smallest index in the node.
  !!   upperBound [in]      -> Vertex of greatest index in the node.
  !!
  !! Result:
  !!   cutVertexIdx -> Index of the vertex used to split the node into children nodes.
  !!
  function computeCutVertexIdx(self, coordinates, verticesIdxs, average, lowerBound, upperBound) &
  result(cutVertexIdx)
    class(node), intent(in)                        :: self
    real(defReal), dimension(:), intent(in)        :: coordinates
    integer(shortInt), dimension(:), intent(inout) :: verticesIdxs
    integer(shortInt), intent(in)                  :: lowerBound, upperBound
    real(defReal), intent(in)                      :: average
    integer(shortInt)                              :: counterUp, counterDown, tempIdx, cutVertexIdx
    
    ! Initialise counters to the values of the supplied lower and upper bounds.
    counterUp = lowerBound
    counterDown = upperBound
    ! Sweep through all vertices between the lower and upper bounds.
    do while (counterUp < counterDown)
      ! If the vertex coordinate along the cut dimension is less than or equal to the supplied
      ! average value of the vertices' coordinates along said dimension, then there is no need to
      ! rearrange the corresponding vertex index and we simply update counterUp.
      if (coordinates(verticesIdxs(counterUp)) <= average) then
        ! Update counterUp.
        counterUp = counterUp + 1
      ! Else, swap the current vertex index with that of a vertex higher up and update counterDown.
      else
        ! Swap indices and update counterDown.
        call swap(verticesIdxs(counterUp), verticesIdxs(counterDown))
        counterDown = counterDown - 1

      end if

    end do
    ! If the vertex coordinate along the cut dimension at the current value of counterUp is less 
    ! than or equal to the supplied average value of the vertices' coordinates along said dimension,
    ! then this vertex is chosen as the cut vertex.
    if (coordinates(verticesIdxs(counterUp)) <= average) then
      cutVertexIdx = counterUp
    ! Else, the vertex immediately before is chosen as the cut vertex.
    else
      cutVertexIdx = counterUp - 1

    end if

  end function computeCutVertexIdx
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  pure recursive subroutine kill(self)
    class(node), intent(inout) :: self
    
    ! Recursively kill children nodes.
    if (allocated(self % left)) then
      call self % left % kill()
      deallocate(self % left)

    end if
    
    if (allocated(self % right)) then
      call self % right % kill()
      deallocate(self % right)

    end if
    
    ! Kill local.
    self % boundingBox = [INF, INF, INF, -INF, -INF, -INF]
    self % cutDimension = 0
    self % lowerBound = 0
    self % upperBound = 0
    self % idx = 0
    self % parentIdx = 0
    self % cutValue = ZERO
    self % cutValue_left = ZERO
    self % cutValue_right = ZERO
    self % hasLeft = .false.
    self % hasRight = .false.

  end subroutine kill
  
  !! Subroutine 'init'
  !!
  !! Basic description:
  !!   Recursively initialises the node and its children nodes.
  !!
  !! Detailed description:
  !!   If the node is a leaf, the subroutine simply computes its bounding box and sets its lower 
  !!   and upper bounds. Else, it computes an initial, approximate bounding box for the node from 
  !!   that of its parent node (if present). Then, the cut dimension for the node is computed from 
  !!   the coordinate with the greatest extent amongst the node's vertices; the cut value is then 
  !!   the average coordinate from all of the node's vertices along the cut dimension. From this cut 
  !!   value, the vertex to use to split the current node into children nodes is computed; finally, 
  !!   the present node's bounding box is recomputed exactly from that of its children nodes.
  !!
  !! Arguments:
  !!   coords [in]          -> Coordinates of all the vertices in the tree.
  !!   verticesIdxs [inout] -> Internal re-ordering of the vertices indices in the tree.
  !!   lowerBound [in]      -> Vertex of smallest index in the node.
  !!   upperBound [in]      -> Vertex of greatest index in the node.
  !!   parent [in]          -> Parent node.
  !!
  recursive subroutine init(self, coords, verticesIdxs, lowerBound, upperBound, isBoundingBoxTree, idx, parent)
    class(node), intent(inout)                                :: self
    real(defReal), dimension(:, :), intent(in)                :: coords
    integer(shortInt), dimension(:), intent(inout)            :: verticesIdxs
    integer(shortInt), intent(in)                             :: lowerBound, upperBound
    logical(defBool), intent(in)                              :: isBoundingBoxTree
    integer(shortInt), intent(inout)                          :: idx
    type(node), intent(in), optional                          :: parent
    integer(shortInt)                                         :: i, j, cutDimension, cutVertexIdx, nVertices, middleIdx
    integer(shortInt), dimension(:), allocatable              :: indicesArray
    integer(shortInt), dimension(upperBound - lowerBound + 1) :: nodeVertexIdxs, sortedNodeVertexIdxs
    real(defReal)                                             :: cutValue, mean, variance, maxVariance, difference
    logical(defBool)                                          :: isChild
    real(defReal), dimension(upperBound - lowerBound + 1)     :: sortedCoords
    
    ! Set the node's index, lower and upper bounds, and compute the number of vertices.
    idx = idx + 1
    self % idx = idx
    self % lowerBound = lowerBound
    self % upperBound = upperBound
    nVertices = upperBound - lowerBound + 1
    nodeVertexIdxs = verticesIdxs(lowerBound:upperBound)

    ! Set the node's index and parent index.
    if (present(parent)) self % parentIdx = parent % idx
    
    ! If nVertices <= bucketSize, the node is a terminal node and there is no need to further subdivide. 
    ! Simply compute the node's bounding box along each dimension and return.
    if (nVertices <= self % bucketSize) then
      do i = 1, 3
        indicesArray = [i]
        if (isBoundingBoxTree) indicesArray = [i, i + 3]
        call self % computeBoundingBox(i, nVertices, coords(indicesArray, nodeVertexIdxs))

      end do
      return
    
    end if

    ! Initialise maxVariance = ZERO, determine if the current node is a child node, then loop over
    ! all dimensions.
    maxVariance = ZERO
    isChild = present(parent)
    do i = 1, size(coords, 1)
      ! If the parent node is allocated and the current dimension is not equal to its cut 
      ! dimension, then the bounding box for the child node is set to its parent's bounding box. 

      ! Compute the mean value along the current dimension and initialise variance = ZERO.
      mean = sum(coords(i, nodeVertexIdxs)) / nVertices
      variance = ZERO
      ! Loop over all vertices and update the variance along the current dimension.
      do j = 1, nVertices
        difference = coords(i, nodeVertexIdxs(j)) - mean
        variance = variance + difference * difference

      end do

      ! If variance > maxVariance, update cutDimension and maxVariance.
      if (variance > maxVariance) then
        cutDimension = i
        maxVariance = variance

      end if

    end do
    ! Set the node's cut dimension.
    self % cutDimension = cutDimension
    
    ! Compute the cut value by averaging the coordinates of all the vertices in the node along the
    ! cut dimension, then set the node's cut value.
    sortedCoords = coords(cutDimension, nodeVertexIdxs)
    sortedNodeVertexIdxs = nodeVertexIdxs
    call quickSort(sortedCoords, sortedNodeVertexIdxs)
    if (mod(nVertices, 2) == 0) then
      middleIdx = nVertices / 2
      cutValue = HALF * (sortedCoords(middleIdx) + sortedCoords(middleIdx + 1))

    else
      middleIdx = nVertices / 2 + 1
      cutValue = sortedCoords(middleIdx)

    end if
    self % cutValue = cutValue
    verticesIdxs(lowerBound:upperBound) = sortedNodeVertexIdxs
    cutVertexIdx = lowerBound + middleIdx - 1
    
    ! Build new children nodes. Catch degenerate cases for which there are no vertices on the left
    ! or right, in which case only a single child node is built. Then, recompute the current 
    ! node's bounding box exactly from its children nodes. If one such children is missing, then 
    ! the current node's bounding box is simply that of its only child node, and the cut value is 
    ! set accordingly. If the two children nodes are present, then the current node's bounding box 
    ! is taken as the average of its two children's bounding boxes.
    allocate(self % left)
    allocate(self % right)
    self % hasLeft = .true.
    self % hasRight = .true.
    call self % left % init(coords, verticesIdxs, lowerBound, cutVertexIdx, isBoundingBoxTree, idx, self)
    call self % right % init(coords, verticesIdxs, cutVertexIdx + 1, upperBound, isBoundingBoxTree, idx, self)
    if (.not. isBoundingBoxTree) then
      self % cutValue_right = self % right % boundingBox(cutDimension)
      self % cutValue_left = self % left % boundingBox(cutDimension + 3)

    end if

    ! Update bounding box from children bounding boxes.
    self % boundingBox(1:3) = min(self % left % boundingBox(1:3), self % right % boundingBox(1:3))
    self % boundingBox(4:6) = max(self % left % boundingBox(4:6), self % right % boundingBox(4:6))

  end subroutine init

  !! Subroutine 'findIntersectedFace'
  !!
  !! Basic description:
  !!   Recursively searches a node and its children nodes and returns the index of the first 
  !!   triangle intersected by a line segment.
  !!
  !! Detailed description:
  !!   Starts at the root node. Then, depending on which side of the cut plane the segment's origin
  !!   is, the subroutine descends into either the left or right child node. Repeats the process
  !!   until a terminal node is encountered, at which point all the triangles containing the 
  !!   vertices in the node are tested for an intersection. If an intersection is found at this 
  !!   point we are guaranteed that it is the closest possible intersection and the subroutine ends
  !!   early. If not, the recursion goes up one level and checks if the line segment intersects the
  !!   cut plane of the parent node. If it does, the other child node of the parent node is visited
  !!   until another terminal node is found. The process is then repeated until either an
  !!   intersection is found or the entire tree has been visited.
  !!
  !! Arguments:
  !!   vertices [in]     -> A vertexShelf.
  !!   faces [in]        -> A faceShelf.
  !!   startPos [in]     -> Origin of the line segment.
  !!   vertexIdxs [in]   -> Internal re-ordering of the vertices in the tree.
  !!   d [inout]         -> Distance to the closest intersected triangle.
  !!   coords [inout]    -> Particle's coordinates.
  !!   edgeIdx [inout]   -> Index of the edge intersected by the line segment.
  !!   vertexIdx [inout] -> Index of the vertex intersected by the line segment.
  !!
  pure recursive subroutine findIntersectedFace(self, vertices, faces, startPos, vertexIdxs, d, coords, edgeIdx, vertexIdx)
    class(node), intent(in)                      :: self
    type(vertexShelf), intent(in)                :: vertices
    type(faceShelf), intent(in)                  :: faces
    real(defReal), dimension(3), intent(in)      :: startPos
    integer(shortInt), dimension(:), intent(in)  :: vertexIdxs
    real(defReal), intent(inout)                 :: d
    type(coord), intent(inout)                   :: coords
    integer(shortInt), intent(inout)             :: edgeIdx, vertexIdx
    integer(shortInt)                            :: cutDimension, i, faceIdx, newEdgeIdx, newVertexIdx
    integer(shortInt), dimension(:), allocatable :: nodeVertices, potentialFaces, faceToElements
    logical(defBool)                             :: isIntersecting
    real(defReal)                                :: cutValue, uComponent, startPosComponent, t, update
    real(defReal), dimension(3)                  :: newStartPos

    ! If the node has left or right children, we need to descend deeper into the tree.
    if (self % hasLeft .or. self % hasRight) then
      ! Retrieve the cut dimension and cut value of the current node.
      cutDimension = self % cutDimension
      cutValue = self % cutValue
      uComponent = coords % dir(cutDimension)
      startPosComponent = startPos(cutDimension)

      ! Initialise t = INF and compute t, which is the fraction of the line segment necessary to reach
      ! the node's cut plane.
      t = INF
      if (.not. areEqual(uComponent, ZERO)) t = (cutValue - startPosComponent) / uComponent
      isIntersecting = ZERO <= t .and. t < norm2(coords % rEnd - startPos)

      if (startPosComponent < cutValue .and. self % hasLeft) then
        call self % left % findIntersectedFace(vertices, faces, startPos, vertexIdxs, d, coords, edgeIdx, vertexIdx)
        
        ! If intersection has been found, or if the line segment does not intersect the node's cut plane,
        ! return early.
        if (coords % elementIdx > 0 .or. (.not. isIntersecting)) return
        
        ! Visit the right child node if it exists and look for an intersection.
        newStartPos = startPos + coords % dir * t
        if (self % hasRight) call self % right % findIntersectedFace(vertices, faces, newStartPos, vertexIdxs, d, &
                                                                     coords, edgeIdx, vertexIdx)

      else if (self % hasRight) then
        call self % right % findIntersectedFace(vertices, faces, startPos, vertexIdxs, d, coords, edgeIdx, vertexIdx)
        
        ! If intersection has been found, or if the line segment does not intersect the node's cut plane,
        ! return early.
        if (coords % elementIdx > 0 .or. (.not. isIntersecting)) return
        
        ! Visit the left child node if it exists and look for an intersection.
        newStartPos = startPos + coords % dir * t
        if (self % hasLeft) call self % left % findIntersectedFace(vertices, faces, newStartPos, vertexIdxs, d, &
                                                                   coords, edgeIdx, vertexIdx)

      end if

    end if

    ! If reached here, we are at a terminal node. Retrieve the vertices in the node and use the tree cache
    ! to re-order the vertices.
    nodeVertices = self % getVertices()
    
    ! Retrieve the faces containing the node vertices from mesh connectivity and remove duplicates.
    potentialFaces = vertices % getVertexFaceIdxs(vertexIdxs(nodeVertices))
    
    ! Loop over all potentially intersected faces.
    do i = 1, size(potentialFaces)
      faceIdx = potentialFaces(i)
      
      ! If the current face is not a boundary face we can cycle to the next face. Else test the face for
      ! an intersection.
      if (.not. faces % getFaceIsBoundary(faceIdx)) cycle
      call faces % computeFaceIntersection(faceIdx, coords % r, coords % rEnd, coords % dir, &
                                           vertices, update, newEdgeIdx, newVertexIdx)
      
      ! If the distance to intersection is less than the lowest distance known update d and the
      ! particle's elementIdx.
      if (update < d) then
        d = update
        faceToElements = faces % getFaceElementIdxs(faceIdx)
        coords % elementIdx = faceToElements(1)
        edgeIdx = newEdgeIdx
        vertexIdx = newVertexIdx

      end if

    end do
    
  end subroutine findIntersectedFace

  !!
  !!
  !!
  recursive subroutine findPotentialElementIdxs(self, data, r, potentialElementIdxs)
    class(node), intent(in)                                     :: self
    real(defReal), dimension(:, :), intent(in)                  :: data
    real(defReal), dimension(3), intent(in)                     :: r
    integer(shortInt), dimension(:), allocatable, intent(inout) :: potentialElementIdxs
    integer(shortInt)                                           :: i, cutDimension
    real(defReal)                                               :: cutValue

    ! First check that the coordinates are within the node's bounding box. Return if not.
    if (any(r < self % boundingBox(1:3)) .or. any(self % boundingBox(4:6) < r)) return
    
    ! If the current node is a leaf, loop through all bounding boxes in the leaf and only 
    ! retain those containing the coordinates.
    if (.not. self % hasLeft .and. .not. self % hasRight) then
      do i = self % lowerBound, self % upperBound
        if (any(r < data(1:3, i)) .or. any(data(4:6, i) < r)) cycle
        call append(potentialElementIdxs, i)

      end do

    else
      ! Retrieve the cut dimension and cut value of the current node.
      cutDimension = self % cutDimension
      cutValue = self % cutValue

      ! Check if the cut dimension is one of the minimum coordinates.
      if (cutDimension <= 3) then
        call self % left % findPotentialElementIdxs(data, r, potentialElementIdxs)
        call self % right % findPotentialElementIdxs(data, r, potentialElementIdxs)

      else
        call self % right % findPotentialElementIdxs(data, r, potentialElementIdxs)
        call self % left % findPotentialElementIdxs(data, r, potentialElementIdxs)

      end if

    end if

  end subroutine findPotentialElementIdxs

  !! Function 'getBoundingBox'
  !!
  !! Basic description:
  !!   Returns the bounding box of the node.
  !!
  !! Result:
  !!   boundingBox -> Bounding box of the node.
  !!
  pure function getBoundingBox(self) result(boundingBox)
    class(node), intent(in)     :: self
    real(defReal), dimension(6) :: boundingBox

    boundingBox = self % boundingBox

  end function getBoundingBox
  
  !! Function 'getVertices'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in the node.
  !!
  !! Result:
  !!   vertices -> Indices of the vertices in the node. Note: these indices are internally sorted so
  !!               they will not correspond to the indices of the vertices in the mesh!
  !!
  pure function getVertices(self) result(vertices)
    class(node), intent(in)                                                 :: self
    integer(shortInt), dimension(self % upperBound - self % lowerBound + 1) :: vertices
    integer(shortInt)                                                       :: i
    
    do i = 0, self % upperBound - self % lowerBound
      vertices(i + 1) = self % lowerBound + i

    end do

  end function getVertices
  
  !! Subroutine 'process'
  !!
  !! Basic description:
  !!   Returns the index of the vertex in the terminal node which is closest to the supplied 3-D
  !!   coordinates.
  !!
  !! Arguments:
  !!   treeData [in]    -> 3-D coordinates of all the vertices in the tree.
  !!   coordinates [in] -> Supplied 3-D coordinates.
  !!   ballsize [inout] -> Smallest search radius reached up to this point.
  !!   idx [out]        -> Index of the vertex closest to the supplied 3-D coordinates.
  !!
  subroutine process(self, treeData, coordinates, ballSize, idx)
    class(node), intent(in)                    :: self
    real(defReal), dimension(:, :), intent(in) :: treeData
    real(defReal), dimension(:), intent(in)    :: coordinates
    real(defReal), intent(inout)               :: ballSize
    integer(shortInt), intent(out)             :: idx
    integer(shortInt)                          :: i, j
    real(defReal)                              :: distanceSquared, diff
    
    ! Loop over all vertices in the terminal node.
    mainLoop: do i = self % lowerBound, self % upperBound
      ! Set the square of the distance between the current vertex and the supplied coordinates to 
      ! ZERO.
      distanceSquared = ZERO
      ! Loop over all dimensions.
      do j = 1, 3
        ! Update distanceSquared.
        diff = treeData(j, i) - coordinates(j)
        distanceSquared = distanceSquared + diff * diff
        ! If distanceSquared is greater than the current lowest distance, move on to the next vertex
        ! in the node.
        if (distanceSquared > ballSize) cycle mainLoop

      end do
      ! Set idx to the index of the vertex corresponding to the current lowest distance and ballSize
      ! to said lowest distance.
      idx = i
      ballSize = distanceSquared

    end do mainLoop

  end subroutine process
  
  !! Subroutine 'search'
  !!
  !! Basic description:
  !!   Searches a node and its children node for the vertex which is closest to the supplied 3-D
  !!   coordinates.
  !!
  !! Detailed description:
  !!   First checks whether the current node is a terminal node. If it is, the node is simply
  !!   processed. If not, the subroutine checks which child node of the current node is closest to
  !!   the supplied coordinates. If it exists, the closer node is then searched. In some cases,
  !!   however, the further node may actually contain a vertex which is closer than any vertices in
  !!   the closer node. In this case, the further node is searched only provided that the square of
  !!   the distance between the coordinates and the further node's bounding box is lower than the
  !!   minimum ball size already returned when searching the closer node.
  !!
  !! Arguments:
  !!   treeData [in]    -> 3-D coordinates of all the vertices in the tree.
  !!   coordinates [in] -> Supplied 3-D coordinates.
  !!   ballsize [inout] -> Smallest search radius reached up to this point.
  !!   idx [out]        -> Index of the vertex closest to the supplied 3-D coordinates. 
  !!
  recursive subroutine search(self, treeData, coordinates, ballSize, idx)
    class(node), intent(in)                    :: self
    real(defReal), dimension(:, :), intent(in) :: treeData
    real(defReal), dimension(:), intent(in)    :: coordinates
    real(defReal), intent(inout)               :: ballSize
    integer(shortInt), intent(out)             :: idx
    integer(shortInt)                          :: i, cutDimension
    real(defReal)                              :: coordinate, distanceSquared
    real(defReal), dimension(6)                :: boundingBox
    
    ! If the current node is a leaf simply process it.
    if (.not. self % hasLeft .and. .not. self % hasRight) then
      call self % process(treeData, coordinates, ballSize, idx)

    else
      ! Set the cut dimension to the node's cut dimension.
      cutDimension = self % cutDimension
      ! Set the coordinate to test for as the value of the supplied coordinates along the cut
      ! dimension.
      coordinate = coordinates(cutDimension)
      ! Depending on whether the coordinate is lower or greater than the node's cut value, the
      ! 'closer' node is the current node's left child node and the 'further' node is the current
      ! node's right child node (and vice-versa).
      ! NOTE: Sorry, I know that the below loop is a bit ugly. A more elegant approach would be to
      !       use pointers, however that would require the 'node' structure to have pointer
      !       components which would then prevent pure procedures to be used when calling a kd-tree.
      if (coordinate < self % cutValue) then
        ! Compute the square of the distance between the 'further' node and the coordinate along the
        ! cut dimension.
        distanceSquared = (self % cutValue_right - coordinate) ** 2
        ! If the 'closer' node is allocated then search it.
        if (allocated(self % left)) call self % left % search(treeData, coordinates, ballSize, idx)
        ! In some cases, the 'further' node away may actually contain vertices which are closer to 
        ! the supplied coordinates. However, if it is allocated it only makes sense to search this 
        ! node if distanceSquared is less than the current ballSize.
        if (allocated(self % right) .and. distanceSquared <= ballSize) then
        
          boundingBox = self % boundingBox
          
          ! Loop over all dimensions.
          do i = 1, 3
            ! Add components for the dimensions other than the cut dimension to distanceSquared.
            if (i /= cutDimension) then
            
              if (coordinates(i) > boundingBox(i + 3)) then
              
                distanceSquared = distanceSquared + (coordinates(i) - boundingBox(i + 3)) ** 2
              
              else if (coordinates(i) < boundingBox(i)) then
              
                distanceSquared = distanceSquared + (boundingBox(i) - coordinates(i)) ** 2
              
              end if
              ! Now if at this point the updated distanceSquared is larger than ballSize, there is 
              ! no need to search the node.
              if (distanceSquared > ballSize) return

            end if

          end do
          ! If reached this point search the node.
          call self % right % search(treeData, coordinates, ballSize, idx)

        end if
      ! Same as above except that the 'closer' and 'further' nodes are swapped.
      else
        distanceSquared = (self % cutValue_left - coordinate) ** 2
        if (self % hasRight) then
          call self % right % search(treeData, coordinates, ballSize, idx)

        end if
        if (self % hasLeft .and. distanceSquared <= ballSize) then
          boundingBox = self % boundingBox
          do i = 1, 3
            if (i /= cutDimension) then
              if (coordinates(i) > boundingBox(i + 3)) then
                distanceSquared = distanceSquared + (coordinates(i) - boundingBox(i + 3)) ** 2
              else if (coordinates(i) < boundingBox(i)) then
                distanceSquared = distanceSquared + (boundingBox(i) - coordinates(i)) ** 2

              end if
              if (distanceSquared > ballSize) return

            end if

          end do
          call self % left % search(treeData, coordinates, ballSize, idx)

        end if

      end if

    end if

  end subroutine search
  
end module node_class