module compoundSurface_inter

  use numPrecision
  use genericProcedures,  only : fatalError, isEqual, numToChar
  use universalVariables, only : INF
  use surface_inter,      only : surface, kill_super => kill

  implicit none
  private

  !!
  !! Abstract interface to group all surfaces containing halfwidths. These surfaces are called
  !! compound surfaces. A compound surface contains one (eg, truncated quadratic surfaces), two
  !! (eg, square cylinders) or three (eg, boxes) halfwidths.
  !!
  !! Private members:
  !!   halfwidths -> Halfwidth(s) of the surface. It is an allocatable array as different compound
  !!                 surfaces can have different numbers of halfwidths.
  !!
  !! Interface:
  !!   closestDist         -> Returns the minimum and maximum distances to the intersected surface 
  !!                          halwidth.
  !!   isHalfspacePositive -> Returns .true. if a particle is in the positive halfspace with respect
  !!                          to the surface's halfwidth(s). Note that it only checks halfspace sign
  !!                          with respect to the components of a surface having halfwidths.
  !!   evaluateHalfwidths  -> Evaluates the expression F(r) = c with respect to the halfwidths of a
  !!                          surface.
  !!   killHalfwidths      -> Returns the halfwidths component to an unitialised state.
  !!   getHalfwidths       -> Returns the halfwidth(s) of the surface.
  !!   setHalfwidths       -> Sets the halfwidth(s) of the surface.
  !!
  type, public, abstract, extends(surface)       :: compoundSurface
    private
    real(defReal), dimension(:), allocatable     :: halfwidths
  contains
    procedure, non_overridable                   :: chooseDistance
    procedure, non_overridable                   :: closestDist
    procedure, non_overridable                   :: isHalfspacePositive
    procedure, non_overridable                   :: evaluateHalfwidths
    procedure                                    :: killHalfwidths
    procedure                                    :: getHalfwidths
    procedure                                    :: setHalfwidths
  end type compoundSurface

contains

  pure subroutine chooseDistance(self, dMin, dMax, surfTolCondition, cannotIntersect, d)
    class(compoundSurface), intent(in)      :: self
    real(defReal), intent(in)               :: dMin, dMax
    logical(defBool), intent(in)            :: surfTolCondition, cannotIntersect
    real(defReal), intent(inout)            :: d

    ! If there is no intersection we can return.
    if (cannotIntersect) return
    
    ! Update d = dMin and retrieve surface tolerance.
    d = dMin
    
    ! Check if particle is on surface or already inside the box and update d = tMax if yes.
    if ((surfTolCondition .and. abs(dMax) >= abs(dMin)) .or. dMin <= ZERO) d = dMax
    
    ! Cap distance to INF if particle is within surfTol to another surface (eg, particle
    ! previously entered the box very close to a vertex) or if d > INF.
    if (d < self % getSurfTol() .or. d > INF) d = INF

  end subroutine chooseDistance

  !! Function 'closestDist'
  !!
  !! Basic description:
  !!   Computes the distance to the nearest intersection between a ray of origin r and direction
  !!   u and a set of halfwidths.
  !!
  !! Detailed description:
  !!   Intersection test is based on the slab method. In essence, an axis-aligned bounding box can
  !!   be represented by two triples (l_x, l_y, l_z) and (h_x, h_y, h_z) which are
  !!   the low and high bounds of the box along each dimension. Each pair of parallel planes (eg,
  !!   x_min and x_max) defines a slab.
  !!
  !!   A point p(d) = [p_x(d), p_y(d), p_z(d)] along the ray can be represented in parametric form
  !!   as:
  !!
  !!   p(d) = r + d * u
  !!
  !!   Assuming that none of the direction components are ZERO, then we can solve for t as:
  !!
  !!   d = (p - r) / u
  !!
  !!   Therefore, the intersection of the ray with the two slabs orthogonal to the i-th dimension
  !!   are given by:
  !!
  !!   d_l_i = (l_i - r_i) / u
  !!   d_h_i = (h_i - r_i) / u
  !!
  !!   The minimum and maximum extrema of the segment within the i-th slab are then:
  !!
  !!   d_min_i = min(d_l_i, d_h_i)
  !!   d_max_i = max(d_l_i, d_h_i)
  !!
  !!   And the intersection of all these segments therefore is:
  !!
  !!   d_min = max(d_min_x, d_min_y, d_min_z)
  !!   d_max = min(d_max_x, d_max_y, d_max_z)
  !!
  !! Arguments:
  !!   offsetCoords [in]     -> Origin-centred coordinates components of the ray along the set of 
  !!                            halfwidths.
  !!   u [in]                -> Direction of the ray.
  !!   dMin [out]            -> Minimum distance of intersection to the halfwidth.
  !!   dMax [out]            -> Maximum distance of intersection to the halfwidth.
  !!   cannotIntersect [out] -> .true. if the ray cannot intersect the halfwidth(s).
  !!
  !! Notes:
  !!   Special handling is required when one or more component(s) of the ray's direction is (are) ZERO.
  !!   In this case, the ray is parallel to one or more pair(s) of orthogonal planes (eg, the ray only
  !!   moves in the y-direction).
  !!
  pure subroutine closestDist(self, offsetCoords, uComponents, dMin, dMax, cannotIntersect)
    class(compoundSurface), intent(in)                            :: self
    real(defReal), dimension(size(self % halfwidths)), intent(in) :: offsetCoords, uComponents
    real(defReal), intent(inout)                                  :: dMin, dMax
    logical(defBool), intent(out)                                 :: cannotIntersect
    integer(shortInt)                                             :: i
    real(defReal)                                                 :: halfwidth, offsetCoord, uComponent, &
                                                                     inverseU

    ! Initialise cannotIntersect = .true.
    cannotIntersect = .true.

    do i = 1, size(self % halfwidths)
      halfwidth = self % halfwidths(i)
      offsetCoord = offsetCoords(i)
      uComponent = uComponents(i)

      ! Perform early check to see if the particle is outside the box and moving 
      ! away from it along the current dimension. If yes the particle cannot
      ! intersect the box and we can return early.
      if ((offsetCoord <= -halfwidth .and. uComponent <= ZERO) .or. &
      (offsetCoord >= halfwidth .and. uComponent >= ZERO)) return

      ! Cycle to next dimension if the current direction component is ZERO, since
      ! in this case there can be no intersection along the current dimension.
      if (isEqual(uComponent, ZERO)) cycle

      ! Pre-compute inverse direction and update dMin and dMax.
      inverseU = 1 / uComponent
      dMin = max(dMin, (-sign(halfwidth, uComponent) - offsetCoord) * inverseU)
      dMax = min(dMax, (sign(halfwidth, uComponent) - offsetCoord) * inverseU)

    end do

    if (dMax <= dMin * (ONE + 10.0_defReal * epsilon(ONE))) return

    ! If reached here, ray intersects the halfwidth(s) so update cannotIntersect = .false.
    cannotIntersect = .false.

  end subroutine closestDist

  !! Function 'isHalfspacePositive'
  !!
  !! Basic description:
  !!   Returns .true. if particle lies in the positive halfspace with respect to a set of halfwidths.
  !!
  !! Detailed description:
  !!   Performs a check along all dimensions along which the surface contains a halfwidth (we need to do
  !!   this to properly account for particles at corners). If the component of the direction vector along
  !!   a given dimension is ZERO, halfspace allocation is decided based on particle's location. If not, it
  !!   is determined by the sign of the projection of the particle's direction component on the halfwidth's
  !!   normal vector.
  !!
  !! Arguments:
  !!   offsetCoords [in] -> Origin-centred coordinate components of the particle along the dimensions 
  !!                        containing halfwidths.
  !!   uComponents [in]  -> Components of the particle's direction vector along the dimensions containing
  !!                        halfwidths.
  !!
  !! Result:
  !!   isIt              -> .true. if particle is in positive halfspace.
  !!
  pure function isHalfspacePositive(self, offsetCoords, uComponents) result(isIt)
    class(compoundSurface), intent(in)                            :: self
    real(defReal), dimension(size(self % halfwidths)), intent(in) :: offsetCoords, uComponents
    logical(defBool)                                              :: isIt
    real(defReal)                                                 :: surfTol, halfwidth, offsetCoord, &
                                                                     dist, uComponent
    integer(shortInt)                                             :: i

    ! Initialise isIt = .true., retrieve surface tolerance then loop over all compound dimensions (note, 
    ! we need to do this to properly account for particles close to corners).
    isIt = .true.
    surfTol = self % getSurfTol()
    do i = 1, size(self % halfwidths)
      ! Retrieve offset coordinates and halfwidth components for the current dimension. Compute
      ! distance to halfwidth along current dimension and cycle to the next dimension if
      ! abs(dist) >= surface tolerance.
      halfwidth = self % halfwidths(i)
      offsetCoord = offsetCoords(i)
      dist = abs(offsetCoord) - halfwidth
      if (abs(dist) >= surfTol) cycle
      
      ! Retrieve direction component for the current dimension. If ray is parallel to current 
      ! dimension update halfspace based on position.
      uComponent = uComponents(i)
      if (isEqual(uComponent, ZERO)) then
        ! If dist >= ZERO the particle is in the positive halfspace and we can return early. If not
        ! move onto the next dimension.
        if (dist >= ZERO) return
        cycle

      end if

      ! If reached here, halfspace is given by the sign of the projection of uComponent onto
      ! halfwidth's normal. Return early if projection > ZERO.
      if (uComponent * sign(ONE, offsetCoord) > ZERO) return

    end do

    ! If reached here, the particle is in the negative halfspace. Update isIt = .false.
    isIt = .false.

  end function isHalfspacePositive

  !! Function 'evaluateHalfwidths'
  !!
  !! Basic description:
  !!   Evaluates the expression F(r) = c for all dimensions containing halfwidths.
  !!
  !! Arguments:
  !!   offsetCoords -> Origin-centred coordinate components of the particle along the dimensions 
  !!                   containing halfwidths.
  !! Result:
  !!   c            -> Evaluated expression F(r).
  !!
  pure function evaluateHalfwidths(self, offsetCoords) result(c)
    class(compoundSurface), intent(in)                            :: self
    real(defReal), dimension(size(self % halfwidths)), intent(in) :: offsetCoords
    real(defReal)                                                 :: c

    c = maxval(abs(offsetCoords) - self % halfwidths)

  end function evaluateHalfwidths

  !! Subroutine 'killHalfwidths'
  !!
  !! Basic description:
  !!   Returns halfwidths to an unitialised state.
  !!
  elemental subroutine killHalfwidths(self)
    class(compoundSurface), intent(inout) :: self

    if (allocated(self % halfwidths)) deallocate(self % halfwidths)

  end subroutine killHalfwidths

  !! Function 'getHalfwidths'
  !!
  !! Basic description:
  !!   Returns the halfwidths of the surface.
  !!
  !! Result:
  !!   halfwidths -> Halfwidths of the surface.
  !!
  pure function getHalfwidths(self) result(halfwidths)
    class(compoundSurface), intent(in)       :: self
    real(defReal), dimension(:), allocatable :: halfwidths

    halfwidths = self % halfwidths

  end function getHalfwidths

  !! Subroutine 'setHalfwidths'
  !!
  !! Basic description:
  !!   Sets the halfwidths of the surface.
  !!
  !! Arguments:
  !!   halfwidths [in] -> Array of halfwidths.
  !!   nRequired [in]  -> Number of required halfwidths for the surface.
  !!
  !! Errors:
  !!   - fatalError if the number of supplied halfwidths does not match the number of required
  !!     halfwidths.
  !!   - fatalError if any of the supplied halfwidths are negative.
  !!
  subroutine setHalfwidths(self, halfwidths, nRequired)
    class(compoundSurface), intent(inout)   :: self
    real(defReal), dimension(:), intent(in) :: halfwidths
    integer(shortInt), intent(in)           :: nRequired
    integer(shortInt)                       :: nHalfwidths
    character(*), parameter                 :: here = 'setHalfwidths (compoundSurface_inter.f90)'

    nHalfwidths = size(halfwidths)
    if (nHalfwidths /= nRequired) call fatalError(here, &
    'Halfwidths must have size: '//numToChar(nRequired)//'. Has: '//numToChar(nHalfwidths)//'.')
    if (any(halfwidths <= ZERO)) call fatalError(here, 'Halfwidths cannot contain negative entries.')

    self % halfwidths = halfwidths

  end subroutine setHalfwidths
  
end module compoundSurface_inter