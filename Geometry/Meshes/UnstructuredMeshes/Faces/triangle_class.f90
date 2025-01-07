module triangle_class
  
  use edgeShelf_class,    only : edgeShelf
  use numPrecision
  use face_inter,         only : face, faceBox, kill_super => kill
  use genericProcedures,  only : areEqual, computeTriangleArea, computeTriangleCentre, computeTriangleNormal
  use universalVariables, only : ONE, ZERO, THIRD, INF
  use vertexShelf_class,  only : vertexShelf
  
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
  type, public, extends(face) :: triangle
    private
  contains
    procedure                 :: computeComponents
    procedure                 :: createTriangle
    procedure                 :: kill
    procedure                 :: testForInclusion
  end type triangle

contains

  !!
  !!
  !!
  pure subroutine computeComponents(self, vertexIdxs, vertices, centroid, normal, area)
    class(triangle), intent(inout)              :: self
    integer(shortInt), dimension(:), intent(in) :: vertexIdxs
    type(vertexShelf), intent(in)               :: vertices
    real(defReal), dimension(3), intent(out)    :: centroid, normal
    real(defReal), intent(out)                  :: area
    integer(shortInt)                           :: i
    real(defReal), dimension(3, 3)              :: array

    ! Construct array.
    do i = 1, 3
      array(i, :) = vertices % getVertexCoordinates(vertexIdxs(i))

    end do

    ! Compute centroid, normal vector and area.
    centroid = computeTriangleCentre(array)
    normal = computeTriangleNormal(array)
    area = computeTriangleArea(normal)
    normal = normal / norm2(normal)

  end subroutine computeComponents

  !!
  !!
  !!
  pure subroutine createTriangle(self, lastNewFaceIdx, edgeIdxs, newVertices, newTriangle, vertexIdxs)
    class(triangle), intent(in)                    :: self
    integer(shortInt), intent(in)                  :: lastNewFaceIdx
    integer(shortInt), dimension(3), intent(in)    :: edgeIdxs
    type(vertexShelf), intent(in)                  :: newVertices
    type(faceBox), intent(inout)                   :: newTriangle
    integer(shortInt), dimension(3), intent(inout) :: vertexIdxs

    ! Allocate new triangle and simply copy everything.
    allocate(triangle :: newTriangle % item)
    call newTriangle % item % init(lastNewFaceIdx, self % getIdx(), self % getIsBoundary(), self % getArea(), &
                                   self % getCentroid(), self % getNormal(), self % getAB(), self % getAC(), &
                                   vertexIdxs, 'Triangle', edgeIdxs)

  end subroutine createTriangle
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(triangle), intent(inout) :: self
    
    ! Face.
    call kill_super(self)

  end subroutine kill

  !! Subroutine 'testForInclusion'
  !!
  !! Basic description:
  !!   Checks whether a point is inside the triangle.
  !!
  !! Detailed description:
  !!   See http://geomalgorithms.com/a06-_intersect-2.html.
  !!
  !! Arguments:
  !!   startPos [in]          -> 3-D coordinates of the line segment's origin.
  !!   endPos [in]            -> 3-D coordinates of the line segment's end.
  !!   u [in]                 -> Particle's direction.
  !!   firstVertexCoords [in] -> 3-D coordinates of the triangle's first vertex.
  !!   dist [out]             -> Distance from the line segment's origin to the point of intersection.
  !!
  pure subroutine testForInclusion(self, vertices, intersectionCoords, diff, d, edgeIdx, vertexIdx)
    class(triangle), intent(in)             :: self
    type(vertexShelf), intent(in)           :: vertices
    real(defReal), dimension(3), intent(in) :: intersectionCoords, diff
    real(defReal), intent(inout)            :: d
    integer(shortInt), intent(inout)        :: edgeIdx, vertexIdx
    integer(shortInt), dimension(3)         :: vertexIdxs
    real(defReal), dimension(3)             :: firstVertexCoords, diffCoords, AB, AC
    real(defReal)                           :: AB_dot_AC, AB_dot_AB, AC_dot_AC, diffCoords_dot_AB, diffCoords_dot_AC, &
                                               inverseDenominator, s, t

    ! Retrieve the indices of the vertices in the triangle and compute the coordinates of the first vertex.
    vertexIdxs = self % getVertexIdxs()
    firstVertexCoords = vertices % getVertexCoordinates(vertexIdxs(1))

    ! Compute the difference between the intersection coordinates and the first vertex. If diffCoords = ZERO
    ! then the intersection point is on the first vertex and we can return early.
    diffCoords = intersectionCoords - firstVertexCoords
    if (areEqual(diffCoords, ZERO)) then
      vertexIdx = vertexIdxs(1)
      d = norm2(diff)
      return

    end if
    
    ! Retrieve the triangle's edge vectors and pre-compute dot products between the different vectors.
    AB = self % getAB()
    AC = self % getAC()
    AB_dot_AC = dot_product(AB, AC)
    AB_dot_AB = dot_product(AB, AB)
    AC_dot_AC = dot_product(AC, AC)
    diffCoords_dot_AB = dot_product(diffCoords, AB)
    diffCoords_dot_AC = dot_product(diffCoords, AC)
    
    ! Pre-compute the denominator and compute the values of s and t, which are the fraction of the
    ! intersection point's projection along each of the two edges sharing the first vertex.
    inverseDenominator = ONE / (AB_dot_AC * AB_dot_AC - AB_dot_AB * AC_dot_AC)
    s = (AB_dot_AC * diffCoords_dot_AC - AC_dot_AC * diffCoords_dot_AB) * inverseDenominator
    t = (AB_dot_AC * diffCoords_dot_AB - AB_dot_AB * diffCoords_dot_AC) * inverseDenominator

    ! If the intersection point's projection along one of the edges sharing the first vertex is
    ! outside said edge we can return early.
    if (s < ZERO .or. t < ZERO .or. s + t > ONE) return

    ! If reached here, the intersection point is inside the triangle. Update d.
    d = norm2(diff)

    ! Check if intersection point lies on a triangle vertex or edge and update vertexIdx and edgeIdx
    ! accordingly.
    if (areEqual(s, ONE)) then
      vertexIdx = vertexIdxs(2)

    elseif (areEqual(t, ONE)) then
      vertexIdx = vertexIdxs(3)

    elseif (areEqual(s, ZERO)) then
      edgeIdx = vertices % findCommonEdgeIdx(vertexIdxs(1), vertexIdxs(3))

    elseif (areEqual(t, ZERO)) then
      edgeIdx = vertices % findCommonEdgeIdx(vertexIdxs(1), vertexIdxs(2))

    elseif (areEqual(s + t, ONE)) then
      edgeIdx = vertices % findCommonEdgeIdx(vertexIdxs(2), vertexIdxs(3))

    end if

  end subroutine testForInclusion

end module triangle_class