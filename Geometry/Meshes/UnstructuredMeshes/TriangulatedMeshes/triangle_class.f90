module triangle_class
  
  use numPrecision
  use genericProcedures,  only : append, areEqual, crossProduct
  use universalVariables, only : ONE, ZERO, THIRD, INF
  
  implicit none
  private
  
  !! Triangle of an OpenFOAM mesh.
  !!
  !! Private members:
  !!   idx             -> Index of the triangle.
  !!   face            -> Index of the parent face from which the triangle originates.
  !!   vertexIdxs      -> Array of indices of the vertices in the triangle.
  !!   tetrahedronIdxs -> Array of indices of the tetrahedra sharing the triangle.
  !!   isBoundary      -> .true. if the triangle is at mesh boundary.
  !!   area            -> Area of the triangle.
  !!   AB              -> First edge vector of the triangle.
  !!   AC              -> Second edge vector of the triangle.
  !!   normal          -> Normal vector of the triangle.
  !!   centre          -> Vector pointing to the centre of the triangle.
  !!
  type, public :: triangle
    private
    integer(shortInt)                            :: idx = 0, faceIdx = 0
    integer(shortInt), dimension(3)              :: vertexIdxs = 0
    integer(shortInt), dimension(:), allocatable :: tetrahedronIdxs
    logical(defBool)                             :: isBoundary = .false.
    real(defReal)                                :: area = ZERO
    real(defReal), dimension(3)                  :: AB = ZERO, AC = ZERO, normal = ZERO, centre = ZERO
  contains
    procedure                                    :: addTetrahedronIdx
    procedure                                    :: computeCentre
    procedure                                    :: computeIntersection
    procedure                                    :: computeNormal
    procedure                                    :: computeArea
    procedure                                    :: hasTetrahedra
    procedure                                    :: getIsBoundary
    procedure                                    :: negateNormal
    procedure                                    :: normaliseNormal
    procedure                                    :: getCentre
    procedure                                    :: getArea
    procedure                                    :: getIdx
    procedure                                    :: getNormal
    procedure                                    :: getTetrahedra
    procedure                                    :: getVertices
    procedure                                    :: kill
    procedure                                    :: setArea
    procedure                                    :: setIsBoundary
    procedure                                    :: setCentre
    procedure                                    :: setFace
    procedure                                    :: setIdx
    procedure                                    :: setNormal
    procedure                                    :: setVertices
  end type triangle

contains
  
  !! Subroutine 'addTetrahedronIdx'
  !!
  !! Basic description:
  !!   Adds the index of a tetrahedron sharing the triangle.
  !!
  !! Arguments:
  !!   tetrahedronIdx [in] -> Index of the tetrahedron.
  !!
  elemental subroutine addTetrahedronIdx(self, tetrahedronIdx)
    class(triangle), intent(inout) :: self
    integer(shortInt), intent(in)  :: tetrahedronIdx
    
    call append(self % tetrahedronIdxs, tetrahedronIdx)
  end subroutine addTetrahedronIdx
  
  !! Subroutine 'computeCentre'
  !!
  !! Basic description:
  !!   Computes the centre of the triangle by taking the arithmetic average of its vertices.
  !!
  !! Arguments:
  !!   A [in] -> Vector pointing to the first vertex of the triangle.
  !!   B [in] -> Vector pointing to the second vertex of the triangle.
  !!   C [in] -> Vector pointing to the third vertex of the triangle.
  !!
  pure subroutine computeCentre(self, A, B, C)
    class(triangle), intent(inout)          :: self
    real(defReal), dimension(3), intent(in) :: A, B, C  
    
    self % centre = THIRD * (A + B + C)
  end subroutine computeCentre
  
  !! Subroutine 'computeIntersection'
  !!
  !! Basic description:
  !!   Checks whether a line segment intersects the triangle and if so, computes the distance from
  !!   the segment's origin to the point of intersection.
  !!
  !! Detailed description:
  !!   See http://geomalgorithms.com/a06-_intersect-2.html.
  !!
  !! Arguments:
  !!   startPos [in]          -> 3-D coordinates of the line segment's origin.
  !!   endPos [in]            -> 3-D coordinates of the line segment's end.
  !!   firstVertexCoords [in] -> 3-D coordinates of the triangle's first vertex.
  !!   isIntersecting [out]   -> .true. if the line segment intersects the triangle.
  !!   dist [out]             -> Distance from the line segment's origin to the point of intersection.
  !!
  pure subroutine computeIntersection(self, startPos, endPos, firstVertexCoords, isIntersecting, d)
    class(triangle), intent(in)             :: self
    real(defReal), dimension(3), intent(in) :: startPos, endPos, firstVertexCoords 
    logical(defBool), intent(out)           :: isIntersecting
    real(defReal), intent(out)              :: d
    real(defReal), dimension(3)             :: normal, diff, intersectionCoords, inPlaneVector, AB, AC
    real(defReal)                           :: denominator, inverseDenominator, s, t, AB_dot_AC, AB_dot_AB, &
                                               AC_dot_AC, inPlane_dot_AB, inPlane_dot_AC
    
    ! Initialise isIntersecting = .false. and d = INF.
    isIntersecting = .false.
    d = INF
    
    ! Retrieve the triangle's normal vector and compute s to check if the line segment intersects
    ! the triangle's plane. Return early if not (s < ZERO or s > ONE).
    normal = self % normal
    diff = endPos - startPos
    denominator = dot_product(normal, diff)
    if (areEqual(denominator, ZERO)) return
    s = dot_product(normal, firstVertexCoords - startPos) / denominator
    if (s < ZERO .or. s > ONE) return
      
    ! Compute the coordinates of the intersection point and create a vector in the plane of the 
    ! triangle going from the first vertex to the point of intersection.
    diff = s * diff
    intersectionCoords = startPos + diff
    inPlaneVector = intersectionCoords - firstVertexCoords
    
    ! Retrieve triangle's edge vectors and pre-compute dot products between the different vectors.
    AB = self % AB
    AC = self % AC
    AB_dot_AC = dot_product(AB, AC)
    AB_dot_AB = dot_product(AB, AB)
    AC_dot_AC = dot_product(AC, AC)
    inPlane_dot_AB = dot_product(inPlaneVector, AB)
    inPlane_dot_AC = dot_product(inPlaneVector, AC)
    
    ! Pre-compute the denominator and compute the values of s and t, which are the fraction of the
    ! intersection point's projection along each of the two edges sharing the first vertex.
    inverseDenominator = ONE / (AB_dot_AC * AB_dot_AC - AB_dot_AB * AC_dot_AC)
    s = (AB_dot_AC * inPlane_dot_AC - AC_dot_AC * inPlane_dot_AB) * inverseDenominator
    t = (AB_dot_AC * inPlane_dot_AB - AB_dot_AB * inPlane_dot_AC) * inverseDenominator

    ! If the intersection point's projection along one of the edges sharing the first vertex is
    ! outside said edge we can return early.
    if (s < ZERO .or. t < ZERO .or. s + t > ONE) return
    
    ! If reached here, update isIntersecting = .true. and compute distance from segment origin to
    ! intersection point.
    isIntersecting = .true.
    d = norm2(diff)

  end subroutine computeIntersection
  
  !! Subroutine 'computeNormal'
  !!
  !! Basic description:
  !!   Computes the normal vector of the triangle by performing the cross-product of its edge 
  !!   vectors (the two vectors contain the common point A). Also sets the triangle's edge vectors
  !!   in the process.
  !!
  !! Arguments:
  !!   A [in] -> Vector pointing to the first vertex of the triangle.
  !!   B [in] -> Vector pointing to the second vertex of the triangle.
  !!   C [in] -> Vector pointing to the third vertex of the triangle.
  !!
  pure subroutine computeNormal(self, A, B, C)
    class(triangle), intent(inout)          :: self
    real(defReal), dimension(3), intent(in) :: A, B, C
    real(defReal), dimension(3)             :: AB, AC
    
    ! Compute the two edge vectors of the triangle, set its edge vectors and
    ! compute its normal vector.
    AB = B - A
    AC = C - A
    self % normal = crossProduct(AB, AC)
  end subroutine computeNormal
  
  !! Subroutine 'computeArea'
  !!
  !! Basic description:
  !!   Computes the area of the triangle from the norm of its normal vector.
  !!
  elemental subroutine computeArea(self)
    class(triangle), intent(inout) :: self
    
    self % area = HALF * norm2(self % normal)
  end subroutine computeArea
  
  !! Function 'hasTetrahedra'
  !!
  !! Basic description:
  !!   Returns .true. if the tetrahedronIdxs component of the triangle is allocated.
  !!
  elemental function hasTetrahedra(self) result(doesIt)
    class(triangle), intent(in) :: self
    logical(defBool)            :: doesIt
    
    doesIt = allocated(self % tetrahedronIdxs)
  end function hasTetrahedra
  
  !! Function 'isBoundary'
  !!
  !! Basic description:
  !!   Returns .true. if the triangle is a boundary triangle.
  !!
  elemental function getIsBoundary(self) result(isBoundary)
    class(triangle), intent(in) :: self
    logical(defBool)            :: isBoundary
    
    isBoundary = self % isBoundary
  end function getIsBoundary
  
  !! Subroutine 'negateNormal'
  !!
  !! Basic description:
  !!   Negates the components of the triangle's normal vector.
  !!
  elemental subroutine negateNormal(self)
    class(triangle), intent(inout) :: self
    
    self % normal = -self % normal
  end subroutine negateNormal
  
  !! Subroutine 'normalisedNormal'
  !!
  !! Basic description:
  !!   Normalises the triangle's normal vector by dividing it by its norm.
  !!
  elemental subroutine normaliseNormal(self)
    class(triangle), intent(inout) :: self
    
    self % normal = self % normal / norm2(self % normal)
  end subroutine normaliseNormal
  
  !! Function 'getCentre'
  !!
  !! Basic description:
  !!   Returns the centre vector of the triangle.
  !!
  !! Result:
  !!   centre -> Vector pointing to the centre of the triangle.
  !!
  pure function getCentre(self) result(centre)
    class(triangle), intent(in) :: self
    real(defReal), dimension(3) :: centre
    
    centre = self % centre
  end function getCentre
  
  !! Function 'getArea'
  !!
  !! Basic description:
  !!   Returns the area of the triangle.
  !!
  !! Result:
  !!   area -> Area of the triangle.
  !!
  elemental function getArea(self) result(area)
    class(triangle), intent(in) :: self
    real(defReal)               :: area
    
    area = self % area
  end function getArea
  
  !! Function 'getNormal'
  !!
  !! Basic description:
  !!   Returns the normal vector of the triangle. Negates (flips) the normal vector if idx < 0.
  !!
  !! Arguments:
  !!   idx [in] [optional] -> Used during particle tracking. Index of the element containing the
  !!                          triangle. If positive then the element is the owner of the triangle
  !!                          and the normal points in the correct direction. Else the element is
  !!                          neighbour of the triangle and the normal's direction needs to be 
  !!                          flipped.
  !!
  !! Result:
  !!   normal              -> Normal vector of the triangle.
  !!
  pure function getNormal(self, idx) result(normal)
    class(triangle), intent(in)             :: self
    integer(shortInt), intent(in), optional :: idx
    real(defReal), dimension(3)             :: normal
    
    normal = self % normal
    
    ! Check if the normal vector needs to be flipped.
    if (.not. present(idx)) return
    if (idx < 0) normal = -normal
    
  end function getNormal
  
  !! Function 'getIdx'
  !!
  !! Basic description:
  !!   Returns the index of the triangle.
  !!
  !! Result:
  !!   idx -> Index of the triangle.
  !!
  elemental function getIdx(self) result(idx)
    class(triangle), intent(in) :: self
    integer(shortInt)           :: idx
    
    idx = self % idx
  end function getIdx
  
  !! Function 'getTetrahedra'
  !!
  !! Basic description:
  !!   Returns the indices of the tetrahedra sharing the triangle.
  !!
  !! Result:
  !!   tetrahedronIdxs -> Indices of the tetrahedra sharing the triangle.
  !!
  pure function getTetrahedra(self) result(tetrahedronIdxs)
    class(triangle), intent(in)                                :: self
    integer(shortInt), dimension(size(self % tetrahedronIdxs)) :: tetrahedronIdxs
    
    tetrahedronIdxs = self % tetrahedronIdxs
  end function getTetrahedra
  
  !! Function 'getVertices'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in the triangle.
  !!
  !! Result:
  !!   vertexIdxs -> Indices of the vertices in the triangle.
  !!
  pure function getVertices(self) result(vertexIdxs)
    class(triangle), intent(in)     :: self
    integer(shortInt), dimension(3) :: vertexIdxs
    
    vertexIdxs = self % vertexIdxs
  end function getVertices
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(triangle), intent(inout) :: self
    
    self % idx = 0
    self % faceIdx = 0
    self % vertexIdxs = 0
    self % isBoundary = .false.
    self % area = ZERO
    self % AB = ZERO
    self % AC = ZERO
    self % normal = ZERO
    self % centre = ZERO
    if (allocated(self % tetrahedronIdxs)) deallocate(self % tetrahedronIdxs)
  end subroutine kill
  
  !! Subroutine 'setArea'
  !!
  !! Basic description:
  !!   Sets the area of the triangle.
  !!
  !! Arguments:
  !!   area [in] -> Area of the triangle.
  !!
  elemental subroutine setArea(self, area)
    class(triangle), intent(inout) :: self
    real(defReal), intent(in)      :: area
    
    self % area = area
  end subroutine setArea
  
  !! Subroutine 'setIsBoundary'
  !!
  !! Basic description:
  !!   Sets isBoundary to .true.
  !!
  elemental subroutine setIsBoundary(self)
    class(triangle), intent(inout) :: self
    
    self % isBoundary = .true.
  end subroutine setIsBoundary
  
  !! Subroutine 'setCentre'
  !!
  !! Basic description:
  !!   Sets the centre of the triangle.
  !!
  !! Arguments:
  !!   centre [in] -> Vector pointing to the centre of the triangle.
  !!
  pure subroutine setCentre(self, centre)
    class(triangle), intent(inout)          :: self
    real(defReal), dimension(3), intent(in) :: centre
    
    self % centre = centre
  end subroutine setCentre
  
  !! Subroutine 'setFace'
  !!
  !! Basic description:
  !!   Sets the index of the face from which the triangle originates.
  !!
  !! Arguments:
  !!   faceIdx [in] -> Index of the face from which the triangle originates.
  !!
  elemental subroutine setFace(self, faceIdx)
    class(triangle), intent(inout) :: self
    integer(shortInt), intent(in)  :: faceIdx
    
    self % faceIdx = faceIdx
  end subroutine setFace
  
  !! Subroutine 'setIdx'
  !!
  !! Basic description:
  !!   Sets the index of the triangle.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the triangle.
  !!
  elemental subroutine setIdx(self, idx)
    class(triangle), intent(inout) :: self
    integer(shortInt), intent(in)  :: idx
    
    self % idx = idx
  end subroutine setIdx
  
  !! Subroutine 'setNormal'
  !!
  !! Basic description:
  !!   Sets the normal vector of the triangle.
  !!
  !! Arguments:
  !!   normal [in] -> Normal vector of the triangle.
  !!
  pure subroutine setNormal(self, normal)
    class(triangle), intent(inout)          :: self
    real(defReal), dimension(3), intent(in) :: normal
    
    self % normal = normal
  end subroutine setNormal
  
  !! Subroutine 'setVertices'
  !!
  !! Basic description:
  !!   Sets the indices of the vertices in the triangle.
  !!
  !! Arguments:
  !!   vertexIdxs [in] -> Array of indices of the vertices in the triangle.
  !!
  pure subroutine setVertices(self, vertexIdxs)
    class(triangle), intent(inout)              :: self
    integer(shortInt), dimension(3), intent(in) :: vertexIdxs
    
    self % vertexIdxs = vertexIdxs
  end subroutine setVertices
end module triangle_class