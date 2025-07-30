module axisAlignedBoundingBox_class

  use genericProcedures,  only : anyAreEqual, areEqual, swap
  use numPrecision
  use publicObjects,      only : intersectionTestPayload, intersectionTestResult, resetIntersectionTestResult
  use universalVariables, only : INF, NUDGE, ONE, SURF_TOL, ZERO

  implicit none
  private

  !!
  !!
  !!
  type, public :: axisAlignedBoundingBox
    private
    real(defReal), dimension(3, 2) :: bounds = ZERO
    real(defReal), dimension(3)    :: centre = ZERO, halfwidths = ZERO
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
    generic            :: intersects => intersects_BoundingBox, intersects_Ray
    procedure, private :: intersects_BoundingBox
    procedure, private :: intersects_Ray
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
  pure subroutine computeBoundsFromBoundingBoxes(self, boundingBoxes)
    class(axisAlignedBoundingBox), intent(inout)           :: self
    type(axisAlignedBoundingBox), dimension(:), intent(in) :: boundingBoxes
    integer(shortInt)                                      :: i, nBoundingBoxes
    real(defReal), dimension(3)                            :: minCoords, maxCoords
    real(defReal), dimension(3, 2)                         :: bounds

    ! Check against zero-sized arrays.
    nBoundingBoxes = size(boundingBoxes)
    if (nBoundingBoxes == 0) return

    bounds = boundingBoxes(1) % getBounds()
    minCoords = bounds(:, 1)
    maxCoords = bounds(:, 2)

    do i = 2, nBoundingBoxes
      bounds = boundingBoxes(i) % getBounds()
      minCoords = min(minCoords, bounds(:, 1))
      maxCoords = max(maxCoords, bounds(:, 2))

    end do
    call self % init(reshape([minCoords, maxCoords], shape=[3, 2]))

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

    doesIt = all(self % bounds(:, 1) <= boundingBox % bounds(:, 1)) .and. all(boundingBox % bounds(:, 2) <= self % bounds(:, 2))

  end function containsBoundingBox

  !!
  !!
  !!
  pure function containsCoords(self, r, u) result(doesIt)
    class(axisAlignedBoundingBox), intent(in) :: self
    real(defReal), dimension(3), intent(in)   :: r, u
    integer(shortInt)                         :: i
    logical(defBool)                          :: doesIt

    do i = 1, 3
      if (areEqual(self % bounds(i, 1), r(i))) then
        doesIt = ZERO < u(i)

      elseif (areEqual(self % bounds(i, 2), r(i))) then
        doesIt = u(i) < ZERO

      else
        doesIt = self % bounds(i, 1) < r(i) .and. r(i) < self % bounds(i, 2)

      end if
      if (.not. doesIt) return

    end do

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
      if (r(i) < self % bounds(i, 1)) then
        diff(i) = self % bounds(i, 1) - r(i)

      elseif (r(i) > self % bounds(i, 2)) then
        diff(i) = r(i) - self % bounds(i, 2)

      end if

    end do
    dSquared = dot_product(diff, diff)

  end function distanceSquared_Vertex

  !!
  !!
  !!
  pure function getBounds(self) result(bounds)
    class(axisAlignedBoundingBox), intent(in) :: self
    real(defReal), dimension(3, 2)            :: bounds

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
    real(defReal), dimension(3, 2), intent(in)   :: bounds

    self % bounds = bounds
    self % centre = HALF * (bounds(:, 1) + bounds(:, 2))
    self % halfwidths = HALF * (bounds(:, 2) - bounds(:, 1))

  end subroutine init

  !!
  !!
  !!
  elemental subroutine intersects_BoundingBox(self, boundingBox, doesIt)
    class(axisAlignedBoundingBox), intent(in) :: self
    type(axisAlignedBoundingBox), intent(in)  :: boundingBox
    logical(defBool), intent(out)             :: doesIt

    doesIt = all(self % bounds(:, 1) <= boundingBox % bounds(:, 2)) .and. all(boundingBox % bounds(:, 1) <= self % bounds(:, 2))

  end subroutine intersects_BoundingBox

  !!
  !!
  !!
  subroutine intersects_Ray(self, payload, result)
    class(axisAlignedBoundingBox), intent(in)    :: self
    class(intersectionTestPayload), intent(in)   :: payload
    class(intersectionTestResult), intent(inout) :: result
    real(defReal)                                :: d, inverseU, tFar, tNear, t1, t2
    integer(shortInt)                            :: i

    ! Initialise d = INF then loop over all halfwidths.
    tNear = -INF
    tFar = INF
    call resetIntersectionTestResult(result)

    do i = 1, 3
      if (areEqual(payload % u(i), ZERO)) then
        if (payload % r(i) < self % bounds(i, 1) .or. self % bounds(i, 2) < payload % r(i)) return

      else
        inverseU = ONE / payload % u(i)
        t1 = (self % bounds(i, 1) - payload % r(i)) * inverseU
        t2 = (self % bounds(i, 2) - payload % r(i)) * inverseU
        if (t2 < t1) call swap(t1, t2)

        tNear = max(tNear, t1)
        tFar = min(tFar, t2)
        if (tFar < tNear .or. tFar < SURF_TOL) return

      end if

    end do

    d = merge(tFar, tNear, tNear < SURF_TOL)
    result % intersects = .true.
    result % d = d

  end subroutine intersects_Ray

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
  subroutine pushFromBoundary(self, u, r, inside)
    class(axisAlignedBoundingBox), intent(in)  :: self
    real(defReal), dimension(3), intent(in)    :: u
    real(defReal), dimension(3), intent(inout) :: r
    logical(defBool), intent(out)              :: inside
    real(defReal), dimension(3)                :: nudgeDirection
    logical(defBool)                           :: isOnBoundary
    integer(shortInt)                          :: i

    ! Initialise hasEscaped = .false. and check if coordinates are on the boundary.
    inside = .true.
    isOnBoundary = anyAreEqual(self % bounds(:, 1), r) .or. anyAreEqual(self % bounds(:, 2), r)

    ! If the coordinates are not on the boundary simply return.
    if (.not. isOnBoundary) return

    ! Nudge the coordinates until they are not on any boundaries anymore.
    do while (isOnBoundary)
      nudgeDirection = ZERO
      do i = 1, 3
        if (areEqual(self % bounds(i, 1), r(i))) then
          if (areEqual(u(i), ZERO)) nudgeDirection(i) = ONE

        elseif (areEqual(self % bounds(i, 2), r(i))) then
          if (areEqual(u(i), ZERO)) nudgeDirection(i) = -ONE

        end if

      end do

      if (any(nudgeDirection /= ZERO)) then
        nudgeDirection = nudgeDirection / norm2(nudgeDirection)

      else
        nudgeDirection = u

      end if
      r = r + nudgeDirection * NUDGE
      isOnBoundary = anyAreEqual(self % bounds(:, 1), r) .or. anyAreEqual(self % bounds(:, 2), r)

    end do

    inside = self % contains(r, u)

  end subroutine pushFromBoundary

end module axisAlignedBoundingBox_class