module tetrahedron_class
  
  use numPrecision
  use genericProcedures,   only : append, areEqual, crossProduct, linFind
  use triangle_class,      only : triangle
  use triangleShelf_class, only : triangleShelf
  use universalVariables,  only : INF, ONE, SURF_TOL, ZERO, targetNotFound
  
  implicit none
  private
  
  !!
  !! Tetrahedron of a given OpenFOAM mesh. Results from the decomposition of polyhedral elements
  !! during the mesh importation process. Consists of a list of vertices of the triangles and 
  !! vertices in the tetrahedron. Also lists the index of the parent element from which the 
  !! tetrahedron originates.
  !!
  !! Private members:
  !!   idx          -> Index of the tetrahedron.
  !!   elementIdx   -> Index of the parent element from which the tetrahedron originates.
  !!   triangleIdxs -> Array of triangles indices making the tetrahedron up.
  !!   vertexIdxs   -> Array of vertices indices making the tetrahedron up.
  !!   volume       -> Volume of the tetrahedron.
  !!   centroid     -> Vector pointing to the centroid of the tetrahedron.
  !!
  type, public                      :: tetrahedron
    private
    integer(shortInt)               :: idx = 0, elementIdx = 0
    integer(shortInt), dimension(4) :: triangleIdxs = 0, vertexIdxs = 0
    integer(shortInt), dimension(6) :: edgeIdxs = 0
    real(defReal)                   :: volume = ZERO
    real(defReal), dimension(3)     :: centroid = ZERO
  contains
    procedure                       :: addEdgeIdx
    procedure                       :: addTriangle
    procedure                       :: computeCentroid
    procedure                       :: computeIntersectedTriangle
    procedure                       :: computePotentialTriangles
    procedure                       :: computeVolume
    procedure                       :: getCentroid
    procedure                       :: getEdgeIdxs
    procedure                       :: getElement
    procedure                       :: getIdx
    procedure                       :: getTriangles
    procedure                       :: getVertices
    procedure                       :: getVolume
    procedure                       :: kill
    procedure                       :: setCentroid
    procedure                       :: setElement
    procedure                       :: setIdx
    procedure                       :: setVertices
    procedure                       :: testForInclusion
  end type tetrahedron

contains

  !! Subroutine 'addEdgeIdx'
  !!
  !! Basic description:
  !!   Adds the index of an edge in the tetrahedron.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the edge.
  !!
  elemental subroutine addEdgeIdx(self, idx)
    class(tetrahedron), intent(inout) :: self
    integer(shortInt), intent(in)     :: idx

    if (linFind(self % edgeIdxs, idx) == targetNotFound) self % edgeIdxs(findloc(self % edgeIdxs, 0)) = idx

  end subroutine addEdgeIdx
  
  !! Subroutine 'addTriangle'
  !!
  !! Basic description:
  !!   Adds the index of a triangle in the tetrahedron.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the triangle.
  !!
  elemental subroutine addTriangle(self, idx)
    class(tetrahedron), intent(inout) :: self
    integer(shortInt), intent(in)     :: idx
    
    self % triangleIdxs(findloc(self % triangleIdxs, 0)) = idx

  end subroutine addTriangle
  
  !! Function 'getCentroid'
  !!
  !! Basic description:
  !!   Returns the centroid of the tetrahedron.
  !!
  !! Result:
  !!   centroid -> A vector pointing to the centroid of the tetrahedron.
  !!
  pure function getCentroid(self) result(centroid)
    class(tetrahedron), intent(in) :: self
    real(defReal), dimension(3)    :: centroid
    
    centroid = self % centroid
    
  end function getCentroid
  
  !! Subroutine 'computeCentroid'
  !!
  !! Basic description:
  !!   Computes the centroid of the tetrahedron.
  !!
  !! Arguments:
  !!   firstVertexCoords [in]  -> Vector pointing to the first vertex of the tetrahedron.
  !!   secondVertexCoords [in] -> Vector pointing to the second vertex of the tetrahedron.
  !!   thirdVertexCoords [in]  -> Vector pointing to the third vertex of the tetrahedron.
  !!   apexCoords [in]         -> Vector pointing to the apex of the tetrahedron.
  !!
  pure subroutine computeCentroid(self, firstVertexCoords, secondVertexCoords, &
                                  thirdVertexCoords, apexCoords)
    class(tetrahedron), intent(inout)       :: self
    real(defReal), dimension(3), intent(in) :: firstVertexCoords, secondVertexCoords, thirdVertexCoords, &
                                               apexCoords
    
    self % centroid = (firstVertexCoords + secondVertexCoords + thirdVertexCoords + apexCoords) / 4.0_defReal

  end subroutine computeCentroid
  
  !! Subroutine 'computeIntersectedTriangle'
  !!
  !! Basic description:
  !!   Computes the tetrahedron's triangle which is intersected by a line segment.
  !!
  !! Detailed description:
  !!   See Macpherson, et al. (2009). DOI: 10.1002/cnm.1128.
  !!
  !! Arguments:
  !!   startPos [in]                -> Beginning of the line segment.
  !!   endPos [in]                  -> End of the line segment.
  !!   potentialTriangleIdxs [in]   -> Array of potential triangles intersected by the line segment.
  !!   intersectedTriangleIdx [out] -> Index of the triangle intersected by the line segment.
  !!   lambda [out]                 -> Fraction of the line segment to be traversed before reaching the
  !!                                   intersection point.
  !!   triangles [in]               -> A triangleShelf.
  !!
  !! TODO: Review this. Especially for cases where a particle may end up on an edge or at a vertex.
  pure subroutine computeIntersectedTriangle(self, startPos, endPos, potentialTriangleIdxs, &
                                             intersectedTriangleIdx, lambda, triangles)
    class(tetrahedron), intent(in)              :: self
    real(defReal), dimension(3), intent(in)     :: startPos, endPos
    integer(shortInt), dimension(:), intent(in) :: potentialTriangleIdxs
    integer(shortInt), intent(out)              :: intersectedTriangleIdx
    real(defReal), intent(out)                  :: lambda
    type(triangleShelf), intent(in)             :: triangles
    integer(shortInt)                           :: i, triangleIdx, absTriangleIdx
    real(defReal), dimension(3)                 :: centre, normal
    real(defReal)                               :: triangleLambda
    
    ! Initialise lambda = INF.
    lambda = INF
    
    ! Loop over all potentially intersected triangles.
    do i = 1, size(potentialTriangleIdxs)
      ! Retrieve the current triangle's normal vector and flip it if necessary.
      triangleIdx = potentialTriangleIdxs(i)
      absTriangleIdx = abs(triangleIdx)
      normal = triangles % getTriangleNormal(triangleIdx)
      
      ! Retrieve the current triangle's centre and compute triangleLambda.
      centre = triangles % getTriangleCentroid(absTriangleIdx)
      triangleLambda = dot_product(centre - startPos, normal) / dot_product(endPos - startPos, normal)
      
      ! If triangleLambda < lambda, update lambda and intersectedTriangleIdx.
      if (triangleLambda < lambda) then
        lambda = triangleLambda
        intersectedTriangleIdx = absTriangleIdx

      end if

    end do

  end subroutine computeIntersectedTriangle
  
  !! Subroutine 'computePotentialTriangles'
  !!
  !! Basic description:
  !!   Computes a list of indices of the potentially intersected triangles using the tetrahedron's 
  !!   centroid and the end of a line segment.
  !!
  !! Detailed description:
  !!   See Macpherson, et al. (2009). DOI: 10.1002/cnm.1128.
  !!
  !! Arguments:
  !!   endPos [in]           -> End of the line segment.
  !!   triangles [in]        -> A triangleShelf.
  !!
  !! Result:
  !!   potentialTriangleIdxs -> Array of indices of the potential triangles intersected by 
  !!                            the line segment.
  !!
  pure function computePotentialTriangles(self, endPos, triangles) result(potentialTriangleIdxs)
    class(tetrahedron), intent(in)               :: self
    real(defReal), dimension(3), intent(in)      :: endPos
    type(triangleShelf), intent(in)              :: triangles
    integer(shortInt), dimension(:), allocatable :: potentialTriangleIdxs
    integer(shortInt)                            :: i, triangleIdx, absTriangleIdx
    real(defReal)                                :: lambda
    real(defReal), dimension(3)                  :: centroid, centre, normal
    
    ! Retrieve tetrahedron's centroid then loop over all triangles in the tetrahedron.
    allocate(potentialTriangleIdxs(0))
    centroid = self % centroid
    do i = 1, 4
      ! Retrieve current triangle index and create its absolute index.
      triangleIdx = self % triangleIdxs(i)
      absTriangleIdx = abs(triangleIdx)
      
      ! Retrieve the current triangle's signed normal vector and centre then compute lambda.
      normal = triangles % getTriangleNormal(triangleIdx)
      centre = triangles % getTriangleCentroid(absTriangleIdx)
      lambda = dot_product(centre - centroid, normal) / dot_product(endPos - centroid, normal)
      
      ! If ZERO <= lambda <= ONE, append the current triangle to the list of potentially intersected triangles.
      if (ZERO <= lambda .and. lambda <= ONE) call append(potentialTriangleIdxs, triangleIdx)

    end do

  end function computePotentialTriangles
  
  !! Subroutine 'computeVolume'
  !!
  !! Basic description:
  !!   Computes the volume of the tetrahedron using vector multiplication.
  !!
  !! Arguments:
  !!   A [in] -> A vector pointing to the first vertex of the tetrahedron.
  !!   B [in] -> A vector pointing to the second vertex of the tetrahedron.
  !!   C [in] -> A vector pointing to the third vertex of the tetrahedron.
  !!   D [in] -> A vector pointing to the fourth vertex of the tetrahedron.
  !!
  pure subroutine computeVolume(self, A, B, C, D)
    class(tetrahedron), intent(inout)       :: self
    real(defReal), dimension(3), intent(in) :: A, B, C, D
    
    self % volume = abs(dot_product(crossProduct(B - A, C - A), D - A)) / 6.0_defReal

  end subroutine computeVolume
  
  !! Function 'getElement'
  !!
  !! Basic description:
  !!   Returns the index of the parent element from which the tetrahedron originates.
  !!
  !! Result:
  !!   elementIdx -> Index of the parent element from which the tetrahedron originates.
  !!
  elemental function getElement(self) result(elementIdx)
    class(tetrahedron), intent(in) :: self
    integer(shortInt)              :: elementIdx
    
    elementIdx = self % elementIdx

  end function getElement

  !! Function 'getEdgeIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the edges of the tetrahedron.
  !!
  !! Result:
  !!   edgeIdxs -> Indices of the edges of the tetrahedron.
  !!
  pure function getEdgeIdxs(self) result(edgeIdxs)
    class(tetrahedron), intent(in)  :: self
    integer(shortInt), dimension(6) :: edgeIdxs
    
    edgeIdxs = self % edgeIdxs

  end function getEdgeIdxs
  
  !! Function 'getIdx'
  !!
  !! Basic description:
  !!   Returns the index of the tetrahedron.
  !!
  !! Result:
  !!   idx -> Index of the tetrahedron.
  !!
  elemental function getIdx(self) result(idx)
    class(tetrahedron), intent(in) :: self
    integer(shortInt)              :: idx
    
    idx = self % idx

  end function getIdx
  
  !! Function 'getTriangles'
  !!
  !! Basic description:
  !!   Returns the indices of the triangles in the tetrahedron.
  !!
  !! Result:
  !!   triangleIdxs -> Array of indices of the triangles in the tetrahedron.
  !!
  pure function getTriangles(self) result(triangleIdxs)
    class(tetrahedron), intent(in)  :: self
    integer(shortInt), dimension(4) :: triangleIdxs
    
    triangleIdxs = self % triangleIdxs

  end function getTriangles
  
  !! Function 'getVertices'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in the tetrahedron.
  !!
  !! Result:
  !!   vertexIdxs -> Array of indices of the vertices in the tetrahedron.
  !!
  pure function getVertices(self) result(vertexIdxs)
    class(tetrahedron), intent(in)  :: self
    integer(shortInt), dimension(4) :: vertexIdxs
    
    vertexIdxs = self % vertexIdxs

  end function getVertices
  
  !! Function 'getVolume'
  !!
  !! Basic description:
  !!   Returns the volume of the tetrahedron.
  !!
  !! Result:
  !!   volume -> Volume of the tetrahedron.
  !!
  elemental function getVolume(self) result(volume)
    class(tetrahedron), intent(in) :: self
    real(defReal)                  :: volume
    
    volume = self % volume
  end function getVolume
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(tetrahedron), intent(inout) :: self
    
    self % idx = 0
    self % elementIdx = 0
    self % volume = ZERO
    self % centroid = ZERO
    self % edgeIdxs = 0
    self % triangleIdxs = 0
    self % vertexIdxs = 0

  end subroutine kill
  
  !! Subroutine 'setCentroid'
  !!
  !! Basic description:
  !!   Sets the centroid of the tetrahedron.
  !!
  !! Arguments:
  !!   centroid [in] -> Centroid of the tetrahedron.
  !!
  pure subroutine setCentroid(self, centroid)
    class(tetrahedron), intent(inout)       :: self
    real(defReal), dimension(3), intent(in) :: centroid
    
    self % centroid = centroid
  end subroutine setCentroid
  
  !! Subroutine 'setElement'
  !!
  !! Basic description:
  !!   Sets the index of the parent element from which the tetrahedron originates.
  !!
  !! Arguments:
  !!   elementIdx [in] -> Index of the parent element from which the tetrahedron originates.
  !!
  elemental subroutine setElement(self, elementIdx)
    class(tetrahedron), intent(inout) :: self
    integer(shortInt), intent(in)     :: elementIdx
    
    self % elementIdx = elementIdx
  end subroutine setElement
  
  !! Subroutine 'setIdx'
  !!
  !! Basic description:
  !!   Sets the index of the tetrahedron.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the tetrahedron.
  !!
  elemental subroutine setIdx(self, idx)
    class(tetrahedron), intent(inout) :: self
    integer(shortInt), intent(in)     :: idx
    
    self % idx = idx
  end subroutine setIdx
  
  !! Subroutine 'setVertices'
  !!
  !! Basic description:
  !!   Sets the indices of the vertices making the tetrahedron up.
  !!
  !! Arguments:
  !!   vertexIdxs [in] -> Array of indices of the vertices making the pyramid up.
  !!
  pure subroutine setVertices(self, vertexIdxs)
    class(tetrahedron), intent(inout)           :: self
    integer(shortInt), dimension(4), intent(in) :: vertexIdxs
    
    self % vertexIdxs = vertexIdxs

  end subroutine setVertices
  
  !! Subroutine 'testForInclusion'
  !!
  !! Basic description:
  !!   Tests whether a set of 3-D coordinates is inside the tetrahedron.
  !!
  !! Detailed description:
  !!   First retrieves the indices of the triangles in the tetrahedron. For each face, the subroutine then 
  !!   checks whether the dot product between the triangle's normal vector and a second vector going 
  !!   from the set of 3-D coordinates to the triangle's centroid is positive. If it is, then the 
  !!   two vectors point in the same direction. If this test is successful for all triangles then 
  !!   the coordinates are inside the tetrahedron.
  !!
  !! Arguments:
  !!   triangles [in]            -> A triangleShelf.
  !!   r [in]                    -> A set of 3-D coordinates.
  !!   failedTriangleIdx [out]   -> Index of the first triangle for which the inclusion test fails.
  !!   surfTolTriangleIdxs [out] -> An array listing the indices of the triangles for which the dot 
  !!                                product is below SURF_TOL, meaning that the coordinates are on 
  !!                                the triangle. It is used in the main tracking routine to assign 
  !!                                a tetrahedron to the coordinates in case the coordinates are on
  !!                                one or more triangle(s).
  !!
  pure subroutine testForInclusion(self, triangles, r, failedTriangleIdx, surfTolTriangleIdxs)
    class(tetrahedron), intent(in)                            :: self
    type(triangleShelf), intent(in)                           :: triangles
    real(defReal), dimension(3), intent(in)                   :: r
    integer(shortInt), intent(out)                            :: failedTriangleIdx
    integer(shortInt), dimension(:), allocatable, intent(out) :: surfTolTriangleIdxs
    integer(shortInt)                                         :: i, triangleIdx, absTriangleIdx
    integer(shortInt), dimension(:), allocatable              :: triangleIdxs
    real(defReal), dimension(3)                               :: normal
    real(defReal)                                             :: dotProduct
    
    ! Initialise failedTriangleIdx = 0, retrieve the indices of the triangles in the tetrahedron 
    ! and loop through all triangles..
    failedTriangleIdx = 0
    triangleIdxs = self % getTriangles()
    do i = 1, 4
      ! Create absolute index of the current triangle.
      triangleIdx = triangleIdxs(i)
      absTriangleIdx = abs(triangleIdx)
      
      ! Retrieve the normal vector and flip it if necessary.
      normal = triangles % getTriangleNormal(triangleIdx)
      
      ! Create the test vector going from the coordinates to the current triangle's centre and
      ! compute the dot product between the test vector and the current triangle's normal vector.
      dotProduct = dot_product(triangles % getTriangleCentroid(absTriangleIdx) - r, normal)

      ! If the coordinates lie on the triangle, update surfTolTriangleIdxs.
      if (areEqual(dotProduct, ZERO)) then
        call append(surfTolTriangleIdxs, absTriangleIdx)
        cycle
      
      end if

      ! If dotProduct < ZERO update failedTriangleIdx and return.
      if (dotProduct < ZERO) then
        failedTriangleIdx = absTriangleIdx
        return

      end if

    end do
    
  end subroutine testForInclusion
  
end module tetrahedron_class