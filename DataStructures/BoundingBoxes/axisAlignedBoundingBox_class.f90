module axisAlignedBoundingBox_class

  use coord_class,        only : coord
  use genericProcedures,  only : anyAreEqual, areEqual
  use numPrecision
  use universalVariables, only : INF, ZERO

  implicit none
  private

  !!
  !!
  !!
  type, public :: axisAlignedBoundingBox
    private
    real(defReal), dimension(6) :: bounds = ZERO
    real(defReal), dimension(3) :: centre = ZERO, halfwidths = ZERO
  contains
    ! Build procedure.
    procedure          :: init
    procedure          :: kill
    ! Runtime procedures.
    generic            :: computeBounds => computeBoundsFromCoords, computeBoundsFromBoundingBoxes
    procedure, private :: computeBoundsFromBoundingBoxes
    procedure, private :: computeBoundsFromCoords
    generic            :: contains => containsCoords, containsBoundingBox
    procedure, private :: containsBoundingBox
    procedure, private :: containsCoords
    generic            :: distanceSquared => distanceSquared_Vertex
    procedure, private :: distanceSquared_Vertex
    procedure          :: getBounds
    procedure          :: getCentre
    procedure          :: getHalfwidths
    generic            :: intersects => intersectsRay, intersectsBoundingBox
    procedure, private :: intersectsBoundingBox
    procedure, private :: intersectsRay
    procedure          :: pushFromBoundary
  end type axisAlignedBoundingBox

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
  !!   nVertices [in]   -> Number of vertices in the node.
  !!   boundingBox [in] -> Bounding boxes of each coordinate in the node.
  !!
  pure subroutine computeBoundsFromBoundingBoxes(self, primitives)
    class(axisAlignedBoundingBox), intent(inout)           :: self
    type(axisAlignedBoundingBox), dimension(:), intent(in) :: primitives
    integer(shortInt)                                      :: i, nPrimitives
    real(defReal), dimension(3)                            :: minCoords, maxCoords
    real(defReal), dimension(6)                            :: bounds

    ! Check against zero-sized arrays.
    nPrimitives = size(primitives)
    if (nPrimitives == 0) return

    bounds = primitives(1) % getBounds()
    minCoords = bounds(1:3)
    maxCoords = bounds(4:6)

    do i = 2, nPrimitives
      bounds = primitives(i) % getBounds()
      minCoords = min(minCoords, bounds(1:3))
      maxCoords = max(maxCoords, bounds(4:6))

    end do
    call self % init([minCoords, maxCoords])

  end subroutine computeBoundsFromBoundingBoxes

  !!
  !!
  !!
  pure subroutine computeBoundsFromCoords(self, primitives)
    class(axisAlignedBoundingBox), intent(inout) :: self
    real(defReal), dimension(:, :), intent(in)   :: primitives

    ! Check against zero-sized arrays.
    if (size(primitives, 2) == 0) return
    call self % init([minval(primitives, dim = 2), maxval(primitives, dim = 2)])

  end subroutine computeBoundsFromCoords

  !!
  !!
  !!
  elemental function containsBoundingBox(self, boundingBox) result(doesIt)
    class(axisAlignedBoundingBox), intent(in) :: self
    type(axisAlignedBoundingBox), intent(in)  :: boundingBox
    logical(defBool)                          :: doesIt

    doesIt = all(self % bounds(1:3) <= boundingBox % bounds(1:3)) .and. all(boundingBox % bounds(4:6) <= self % bounds(4:6))

  end function containsBoundingBox

  !!
  !!
  !!
  pure function containsCoords(self, r) result(doesIt)
    class(axisAlignedBoundingBox), intent(in) :: self
    real(defReal), dimension(3), intent(in)   :: r
    logical(defBool)                          :: doesIt

    doesIt = all(self % bounds(1:3) <= r(1:3)) .and. all(r(1:3) <= self % bounds(4:6))

  end function containsCoords

  !!
  !!
  !!
  pure function distanceSquared_Vertex(self, r) result(dSquared)
    class(axisAlignedBoundingBox), intent(in) :: self
    real(defReal), dimension(3), intent(in)   :: r
    real(defReal)                             :: dSquared
    real(defReal), dimension(3)               :: diff
    integer(shortInt)                         :: i

    ! Initialise diff = ZERO then loop over all dimensions.
    diff = ZERO
    do i = 1, 3
      if (r(i) < self % bounds(i)) then
        diff(i) = self % bounds(i) - r(i)

      elseif (r(i) > self % bounds(i + 3)) then
        diff(i) = r(i) - self % bounds(i + 3)

      end if

    end do
    dSquared = dot_product(diff, diff)

  end function distanceSquared_Vertex

  !!
  !!
  !!
  pure function getBounds(self) result(bounds)
    class(axisAlignedBoundingBox), intent(in) :: self
    real(defReal), dimension(6)               :: bounds

    bounds = self % bounds

  end function getBounds

  !!
  !!
  !!
  pure function getCentre(self) result(centre)
    class(axisAlignedBoundingBox), intent(in) :: self
    real(defReal), dimension(3)               :: centre

    centre = self % centre

  end function getCentre
  
  !!
  !!
  !!
  pure function getHalfwidths(self) result(halfwidths)
    class(axisAlignedBoundingBox), intent(in) :: self
    real(defReal), dimension(3)               :: halfwidths

    halfwidths = self % halfwidths

  end function getHalfwidths 

  !!
  !!
  !!
  pure subroutine init(self, bounds)
    class(axisAlignedBoundingBox), intent(inout) :: self
    real(defReal), dimension(6), intent(in)      :: bounds

    self % bounds = bounds
    self % centre = HALF * (bounds(1:3) + bounds(4:6))
    self % halfwidths = HALF * (bounds(4:6) - bounds(1:3))

  end subroutine init

  !!
  !!
  !!
  elemental function intersectsBoundingBox(self, boundingBox) result(doesIt)
    class(axisAlignedBoundingBox), intent(in) :: self
    type(axisAlignedBoundingBox), intent(in)  :: boundingBox
    logical(defBool)                          :: doesIt

    doesIt = all(self % bounds(1:3) <= boundingBox % bounds(4:6)) .and. all(boundingBox % bounds(1:3) <= self % bounds(4:6))

  end function intersectsBoundingBox

  !!
  !!
  !!
  pure function intersectsRay(self, startCoords, endCoords) result(doesIt)
    class(axisAlignedBoundingBox), intent(in) :: self
    real(defReal), dimension(3), intent(in)   :: startCoords, endCoords
    logical(defBool)                          :: doesIt

    doesIt = .true.

  end function intersectsRay

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(axisAlignedBoundingBox), intent(inout) :: self

    self % bounds = ZERO
    self % centre = ZERO
    self % halfwidths = ZERO

  end subroutine kill

  !!
  !!
  !! Note:
  !!   Assumes that the coordinates are already inside the bounding box to begin with.
  !!
  subroutine pushFromBoundary(self, coords, inside)
    class(axisAlignedBoundingBox), intent(in) :: self
    type(coord), intent(inout)                :: coords
    logical(defBool), intent(out)             :: inside
    real(defReal), dimension(3)               :: r, u, nudgeDirection
    logical(defBool)                          :: isOnBoundary
    integer(shortInt)                         :: i

    ! Initialise hasEscaped = .false. and check if coordinates are on the boundary.
    inside = .true.
    r = coords % getPositionToNudge()
    isOnBoundary = anyAreEqual(self % bounds(1:3), r) .or. anyAreEqual(self % bounds(4:6), r)

    ! If the coordinates are not on the boundary simply return.
    if (.not. isOnBoundary) return

    ! Nudge the coordinates until they are not on any boundaries anymore.
    u = coords % getDirection()
    do while (isOnBoundary)
      nudgeDirection = ZERO
      do i = 1, 3
        if (areEqual(self % bounds(i), r(i))) then
          if (areEqual(u(i), ZERO)) nudgeDirection(i) = ONE

        elseif (areEqual(self % bounds(i + 3), r(i))) then
          if (areEqual(u(i), ZERO)) nudgeDirection(i) = -ONE

        end if

      end do

      if (any(nudgeDirection /= ZERO)) then
        call coords % nudgePosition(nudgeDirection / norm2(nudgeDirection))

      else
        call coords % nudgePosition()

      end if

      r = coords % getPositionToNudge()
      isOnBoundary = anyAreEqual(self % bounds(1:3), r) .or. anyAreEqual(self % bounds(4:6), r)

    end do

    inside = self % contains(coords % getPositionToNudge())

  end subroutine pushFromBoundary

end module axisAlignedBoundingBox_class