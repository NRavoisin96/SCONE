module pyramid_class
  
  use numPrecision
  use genericProcedures,      only : append, findCommon
  use edgeShelf_class,        only : edgeShelf
  use face_class,             only : face
  use faceShelf_class,        only : faceShelf
  use triangle_class,         only : triangle
  use triangleShelf_class,    only : triangleShelf
  use tetrahedronShelf_class, only : tetrahedronShelf
  use universalVariables,     only : ZERO
  use vertex_class,           only : vertex
  use vertexShelf_class,      only : vertexShelf
  
  implicit none
  private
  
  !!
  !! Pyramid of an OpenFOAM mesh. Intermediate product in the decomposition of polyhedral meshes 
  !! into tetrahedral ones (element -> pyramids -> tetrahedra). Consists of a list of vertex and 
  !! triangle indices. Also lists the index of the parent element from which the pyramid
  !! has been created as well as the index of the face forming the pyramid's base.
  !!
  !! Private members:
  !!   idx          -> Index of the pyramid.
  !!   elementIdx   -> Index of the parent element from which the pyramid originates.
  !!   faceIdx      -> Index of the face defining the pyramid's base.
  !!   triangleIdxs -> Array of triangles indices making the pyramid up.
  !!   vertexIdxs   -> Array of vertices indices making the pyramid up.
  !!   volume       -> Volume of the pyramid.
  !!   centroid     -> Vector pointing to the pyramid's centroid.
  !!
  type, public :: pyramid
    private
    integer(shortInt)                            :: idx = 0, elementIdx = 0, faceIdx = 0
    integer(shortInt), dimension(:), allocatable :: triangleIdxs, vertexIdxs
    real(defReal)                                :: volume = ZERO
    real(defReal), dimension(3)                  :: centroid
  contains
    procedure                                    :: addTriangle
    procedure                                    :: computeVolume
    procedure                                    :: computeCentroid
    procedure                                    :: getVolume
    procedure                                    :: getCentroid
    procedure                                    :: getFace
    procedure                                    :: getTriangles
    procedure                                    :: getVertices
    procedure                                    :: kill
    procedure                                    :: setFace
    procedure                                    :: setIdx
    procedure                                    :: setElement
    procedure                                    :: setVertices
    procedure                                    :: split
  end type pyramid

contains
  
  !! Subroutine 'addTriangle'
  !!
  !! Basic description:
  !!   Adds the index of a triangle in the pyramid.
  !!
  !! Arguments:
  !!   triangleIdx [in] -> Index of the triangle.
  !!
  elemental subroutine addTriangle(self, triangleIdx)
    class(pyramid), intent(inout) :: self
    integer(shortInt), intent(in) :: triangleIdx
    
    call append(self % triangleIdxs, triangleIdx)
  end subroutine addTriangle
  
  !! Subroutine 'computeCentroid'
  !!
  !! Basic description:
  !!   Computes the centroid of the pyramid from the centroid of its base and its apex.
  !!
  !! Arguments:
  !!   faces [in]       -> A faceShelf.
  !!   apexCoords [in]  -> 3-D coordinates of the pyramid's apex.
  !!
  pure subroutine computeCentroid(self, faces, apexCoords)
    class(pyramid), intent(inout)           :: self
    type(faceShelf), intent(in)             :: faces
    real(defReal), dimension(3), intent(in) :: apexCoords

    self % centroid = (3.0_defReal * faces % getFaceCentroid(abs(self % faceIdx)) + apexCoords) / 4.0_defReal

  end subroutine computeCentroid
  
  !! Subroutine 'computeVolume'
  !!
  !! Basic description:
  !!   Computes the volume of the pyramid from the normal vector of its base and its height.
  !!
  !! Arguments:
  !!   base [in] -> Face making the pyramid's base.
  !!   apex [in] -> Vector pointing to the pyramid's apex.
  !!
  pure subroutine computeVolume(self, base, apex)
    class(pyramid), intent(inout)           :: self
    type(face), intent(in)                  :: base
    real(defReal), dimension(3), intent(in) :: apex

    ! Compute volume. Note that we need to multiply the normal vector of the base by its area
    ! to un-normalise the normal vector and get the correct volume for the pyramid.
    self % volume = THIRD * abs(dot_product(base % getNormal() * base % getArea(), apex - base % getCentroid()))

  end subroutine computeVolume
  
  !! Function 'getVolume'
  !!
  !! Basic description:
  !!   Returns the volume of the pyramid.
  !!
  !! Result:
  !!   volume -> Volume of the pyramid.
  !!
  elemental function getVolume(self) result(volume)
    class(pyramid), intent(in) :: self
    real(defReal)              :: volume
    
    volume = self % volume
  end function getVolume
  
  !! Function 'getCentroid'
  !!
  !! Basic description:
  !!   Returns the centroid of the pyramid.
  !!
  !! Result:
  !!   centroid -> A vector pointing to the centroid of the pyramid.
  !!
  pure function getCentroid(self) result(centroid)
    class(pyramid), intent(in)  :: self
    real(defReal), dimension(3) :: centroid
    
    centroid = self % centroid
  end function getCentroid
  
  !! Function 'getFace'
  !!
  !! Basic description:
  !!   Returns the index of the face making the pyramid's base.
  !!
  !! Result:
  !!   faceIdx -> Index of the pyramid's base.
  !!
  elemental function getFace(self) result(faceIdx)
    class(pyramid), intent(in) :: self
    integer(shortInt)          :: faceIdx
    
    faceIdx = self % faceIdx
  end function getFace
  
  !! Function 'getTriangles'
  !!
  !! Basic description:
  !!   Returns the indices of the triangles in the pyramid.
  !!
  !! Result:
  !!   triangleIdxs -> Array of indices of the triangles in the pyramid.
  !!
  pure function getTriangles(self) result(triangleIdxs)
    class(pyramid), intent(in)                           :: self
    integer(shortInt), dimension(size(self % triangleIdxs)) :: triangleIdxs
    
    triangleIdxs = self % triangleIdxs
  end function getTriangles
  
  !! Function 'getVertices'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in the pyramid.
  !!
  !! Result:
  !!   vertexIdxs -> Array of indices of the vertices in the pyramid.
  !!
  pure function getVertices(self) result(vertexIdxs)
    class(pyramid), intent(in)                          :: self
    integer(shortInt), dimension(size(self % vertexIdxs)) :: vertexIdxs
    
    vertexIdxs = self % vertexIdxs
  end function getVertices
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(pyramid), intent(inout) :: self
    
    self % idx = 0
    self % faceIdx = 0
    self % elementIdx = 0
    self % volume = ZERO
    self % centroid = ZERO
    if (allocated(self % vertexIdxs)) deallocate(self % vertexIdxs)
    if (allocated(self % triangleIdxs)) deallocate(self % triangleIdxs)
  end subroutine kill
  
  !! Subroutine 'setFace'
  !!
  !! Basic description:
  !!   Sets the index of the face making the pyramid's base.
  !!
  !! Arguments:
  !!   faceIdx [in] -> Index of the pyramid's base.
  !!
  elemental subroutine setFace(self, faceIdx)
    class(pyramid), intent(inout) :: self
    integer(shortInt), intent(in) :: faceIdx
    
    self % faceIdx = faceIdx
  end subroutine setFace
  
  !! Subroutine 'setIdx'
  !!
  !! Basic description:
  !!   Sets the index of the pyramid.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the pyramid.
  !!
  elemental subroutine setIdx(self, idx)
    class(pyramid), intent(inout) :: self
    integer(shortInt), intent(in) :: idx
    
    self % idx = idx
  end subroutine setIdx
  
  !! Subroutine 'setElement'
  !!
  !! Basic description:
  !!   Sets the index of the element from which the pyramid originates.
  !!
  !! Arguments:
  !!   elementIdx [in] -> Index of the pyramid's parent element.
  !!
  elemental subroutine setElement(self, elementIdx)
    class(pyramid), intent(inout) :: self
    integer(shortInt), intent(in) :: elementIdx
    
    self % elementIdx = elementIdx
  end subroutine setElement
  
  !! Subroutine 'setVertices'
  !!
  !! Basic description:
  !!   Sets the indices of the vertices in the pyramid.
  !!
  !! Arguments:
  !!   vertexIdxs [in] -> Array of indices of the vertices in the pyramid.
  !!
  pure subroutine setVertices(self, vertexIdxs)
    class(pyramid), intent(inout)               :: self
    integer(shortInt), dimension(:), intent(in) :: vertexIdxs
    
    self % vertexIdxs = vertexIdxs
  end subroutine setVertices
  
  !! Subroutine 'split'
  !!
  !! Basic description:
  !!   Splits the pyramid into tetrahedra.
  !!
  !! Detailed description:
  !!   Creates a number of tetrahedra equal to the number of triangles forming the pyramid's base.
  !!   One of the four triangles forming a given tetrahedron is one of the base's triangles; two 
  !!   other triangles are given by two of the pyramid's sides which were created during element
  !!   splitting. The remaining triangle is an internal triangle which is generated here by joining
  !!   two of the base's vertices to the pyramid's apex.
  !!
  !! Arguments:
  !!   base [in]                  -> Face making the pyramid's base up.
  !!   vertices [in]              -> A vertexShelf.
  !!   triangles [inout]          -> A triangleShelf.
  !!   tetrahedra [inout]         -> A tetrahedronShelf.
  !!   lastTriangleIdx [inout]    -> Index of the last item in triangleShelf.
  !!   lastTetrahedronIdx [inout] -> Index of the last item in tetrahedronShelf.
  !!
  elemental subroutine split(self, faces, edges, vertices, triangles, tetrahedra, lastTriangleIdx, lastTetrahedronIdx)
    class(pyramid), intent(in)                   :: self
    type(faceShelf), intent(in)                  :: faces
    type(edgeShelf), intent(inout)               :: edges
    type(vertexShelf), intent(inout)             :: vertices
    type(triangleShelf), intent(inout)           :: triangles
    type(tetrahedronShelf), intent(inout)        :: tetrahedra
    integer(shortInt), intent(inout)             :: lastTriangleIdx, lastTetrahedronIdx
    integer(shortInt), dimension(3)              :: triangleVertexIdxs
    integer(shortInt), dimension(4)              :: tetrahedronVertexIdxs
    integer(shortInt), dimension(:), allocatable :: faceVertexIdxs, faceTriangleIdxs, edgeIdxs
    integer(shortInt)                            :: apexIdx, faceIdx, i, j, k, minVertexIdx, nTrianglesFound, &
                                                    minVertexLoc, nVertices, absTriangleIdx, edgeIdx
    real(defReal), dimension(3)                  :: minVertexCoords, apexCoords
    real(defReal), dimension(3, 3)               :: triangleCoordsArray
    real(defReal), dimension(4, 3)               :: tetrahedronCoordsArray
    
    ! Retrieve the index of the pyramid's base and indices of the vertices and triangles in the base.
    faceIdx = self % faceIdx
    faceVertexIdxs = faces % getFaceVertexIdxs(abs(faceIdx))
    faceTriangleIdxs = faces % getFaceTriangleIdxs(abs(faceIdx))
    
    ! Compute the number of vertices in the base and find the index and location of the vertex in the
    ! base with the smallest index.
    nVertices = size(faceVertexIdxs)
    minVertexLoc = minloc(faceVertexIdxs, 1)
    minVertexIdx = minval(faceVertexIdxs)
    
    ! Create vectors pointing to minVertexIdx and to the pyramid's apex (note that due to the earlier
    ! splitting of elements this is always the last vertex in the pyramid's list of vertices).
    minVertexCoords = vertices % getVertexCoordinates(minVertexIdx)
    triangleCoordsArray(1, :) = minVertexCoords
    tetrahedronCoordsArray(1, :) = minVertexCoords
    
    apexIdx = self % vertexIdxs(size(self % vertexIdxs))
    apexCoords = vertices % getVertexCoordinates(apexIdx)
    triangleCoordsArray(3, :) = apexCoords
    tetrahedronCoordsArray(4, :) = apexCoords
    
    ! Pre-set two vertices of the new triangles since they are shared and loop through all new tetrahedra.
    triangleVertexIdxs(1) = minVertexIdx
    triangleVertexIdxs(3) = apexIdx
    do i = 1, nVertices - 2
      ! Update freeTetrahedronIdx and create the new tetrahedron vertex indices.
      lastTetrahedronIdx = lastTetrahedronIdx + 1
      tetrahedronVertexIdxs = [triangles % getTriangleVertexIdxs(faceTriangleIdxs(i)), apexIdx]
      tetrahedronCoordsArray(2, :) = vertices % getVertexCoordinates(tetrahedronVertexIdxs(2))
      tetrahedronCoordsArray(3, :) = vertices % getVertexCoordinates(tetrahedronVertexIdxs(3))

      ! Initialise the new tetrahedron in the shelf.
      call tetrahedra % buildTetrahedron(lastTetrahedronIdx, self % elementIdx, tetrahedronVertexIdxs, &
                                         tetrahedronCoordsArray)

      ! Update mesh connectivity information.
      do j = 1, 4
        call vertices % addTetrahedronIdxToVertex(tetrahedronVertexIdxs(j), lastTetrahedronIdx)

      end do
                                                          
      ! Loop through all supplied triangles to determine which ones are in the current tetrahedron.
      ! Initialise nTrianglesFound = 0.
      nTrianglesFound = 0
      do j = 1, size(self % triangleIdxs)
        ! Cycle to the next triangle if the current triangle has less than three vertices in common with the 
        ! current tetrahedron.
        absTriangleIdx = abs(self % triangleIdxs(j))
        if (size(findCommon(tetrahedronVertexIdxs, triangles % getTriangleVertexIdxs(absTriangleIdx))) < 3) cycle

        ! Update nTrianglesFound if reached here and update mesh connectivity.
        nTrianglesFound = nTrianglesFound + 1

        ! If the current triangle is already associated to a tetrahedron then the current tetrahedron is 
        ! its neighbour.
        if (triangles % getTriangleHasTetrahedra(absTriangleIdx)) then
          call tetrahedra % addTriangleIdxToTetrahedron(lastTetrahedronIdx, -absTriangleIdx)

        else
          call tetrahedra % addTriangleIdxToTetrahedron(lastTetrahedronIdx, absTriangleIdx)

        end if

        call triangles % addTetrahedronIdxToTriangle(absTriangleIdx, lastTetrahedronIdx)
        edgeIdxs = triangles % getTriangleEdgeIdxs(absTriangleIdx)
        do k = 1, 3
          edgeIdx = edgeIdxs(k)
          call edges % addTetrahedronIdxToEdge(edgeIdx, lastTetrahedronIdx)
          call tetrahedra % addEdgeIdxToTetrahedron(lastTetrahedronIdx, edgeIdx)

        end do

        if (nTrianglesFound == 3) exit

      end do

      ! Return if we are at the last iteration of this loop (no need to create a new triangle).
      if (i == nVertices - 2) return
      
      ! Create a new internal triangle. Update freeTriangleIdx, Set the second vertex in the new 
      ! triangle and create a vector pointing to it. Compute the triangle's centre, normal vector 
      ! and area. Normalise normal vector after area computation.
      lastTriangleIdx = lastTriangleIdx + 1
      triangleVertexIdxs(2) = faceVertexIdxs(mod(minVertexLoc + i, nVertices) + 1)
      triangleCoordsArray(2, :) = vertices % getVertexCoordinates(triangleVertexIdxs(2))
      call triangles % buildTriangle(lastTriangleIdx, triangleVertexIdxs, 0, .false., triangleCoordsArray, &
                                     tetrahedra % getTetrahedronCentroid(lastTetrahedronIdx))
      
      ! Update mesh connectivity information. Current tetrahedron is the new triangle's owner,
      ! next tetrahedron is its neighbour.
      call tetrahedra % addTriangleIdxToTetrahedron(lastTetrahedronIdx, lastTriangleIdx)
      call tetrahedra % addTriangleIdxToTetrahedron(lastTetrahedronIdx + 1, -lastTriangleIdx)
      call triangles % addTetrahedronIdxToTriangle(lastTriangleIdx, lastTetrahedronIdx)
      call triangles % addTetrahedronIdxToTriangle(lastTriangleIdx, lastTetrahedronIdx + 1)

      do j = 1, 3
        edgeIdx = vertices % findCommonEdgeIdx(triangleVertexIdxs(j), triangleVertexIdxs(mod(j, 3) + 1))
        call edges % addTetrahedronIdxToEdge(edgeIdx, lastTetrahedronIdx)
        call edges % addTetrahedronIdxToEdge(edgeIdx, lastTetrahedronIdx + 1)
        call edges % addTriangleIdxToEdge(edgeIdx, lastTriangleIdx)
        call tetrahedra % addEdgeIdxToTetrahedron(lastTetrahedronIdx, edgeIdx)
        call tetrahedra % addEdgeIdxToTetrahedron(lastTetrahedronIdx + 1, edgeIdx)
        call triangles % addEdgeIdxToTriangle(lastTriangleIdx, edgeIdx)
        call vertices % addTriangleIdxToVertex(triangleVertexIdxs(j), lastTriangleIdx)

      end do
      
    end do

  end subroutine split

end module pyramid_class