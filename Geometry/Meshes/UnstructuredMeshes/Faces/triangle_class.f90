module triangle_class
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use edgeShelf_class,              only : edgeShelf
  use numPrecision
  use face_inter,                   only : face, faceBox, kill_super => kill
  use genericProcedures,            only : areEqual, computeTriangleArea, computeTriangleCentre, computeTriangleNormal
  use universalVariables,           only : ONE, ZERO, THIRD, INF
  use vertexShelf_class,            only : vertexShelf
  
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
  end type triangle

contains

  !!
  !!
  !!
  subroutine computeComponents(self, vertexIdxs, vertices, centroid, normal, area)
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
  pure subroutine createTriangle(self, lastNewFaceIdx, edgeIdxs, newVertices, newTriangle, vertexIdxs, boundingBox)
    class(triangle), intent(in)                    :: self
    integer(shortInt), intent(in)                  :: lastNewFaceIdx
    integer(shortInt), dimension(3), intent(in)    :: edgeIdxs
    type(vertexShelf), intent(in)                  :: newVertices
    type(faceBox), intent(inout)                   :: newTriangle
    integer(shortInt), dimension(3), intent(inout) :: vertexIdxs
    type(axisAlignedBoundingBox), intent(in)       :: boundingBox

    ! Allocate new triangle and simply copy everything.
    allocate(triangle :: newTriangle % item)
    call newTriangle % item % init(lastNewFaceIdx, self % getIdx(), self % getIsBoundary(), self % getArea(), &
                                   self % getCentroid(), self % getNormal(), self % getAB(), self % getAC(), &
                                   vertexIdxs, 'Triangle', boundingBox, edgeIdxs)

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

end module triangle_class