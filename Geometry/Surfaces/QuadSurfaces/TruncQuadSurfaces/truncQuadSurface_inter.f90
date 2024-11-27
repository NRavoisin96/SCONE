module truncQuadSurface_inter

  use numPrecision
  use universalVariables, only : INF, MISS_TOL, VACUUM_BC, REFLECTIVE_BC, PERIODIC_BC
  use genericProcedures,  only : isEqual, numToChar
  use quadSurface_inter,  only : quadSurface

  implicit none
  private

  !!
  !! Abstract interface to group all truncated quadratic surfaces.
  !!
  type, public, abstract, extends(quadSurface) :: truncQuadSurface
    private
    real(defReal)                              :: halfwidth = ZERO, halfwidthOrigin = ZERO
    integer(shortInt), dimension(2)            :: BCs = VACUUM_BC
  contains
    procedure, non_overridable                 :: distanceHalfwidth
    procedure, non_overridable                 :: isHalfspacePositive
    procedure, non_overridable                 :: evaluateHalfwidth
    procedure, non_overridable                 :: setHalfwidthBCs
    procedure, non_overridable                 :: killHalfwidth
    procedure, non_overridable                 :: setHalfwidth
    procedure, non_overridable                 :: explicitHalfwidthBCs
    procedure, non_overridable                 :: transformHalfwidthBCs
  end type truncQuadSurface

contains

  !! Subroutine 'distanceHalfwidth'
  !!
  !! Basic description:
  !!   Computes the distance to intersection to the quadratic surface given by the solutions to the
  !!   quadratic equation a * d^2 + 2 * k * d + c = 0.
  !!
  !! Arguments:
  !!   a [in]                         -> Coefficient of d^2 in a * d^2 + 2 * k * d + c = 0.
  !!   k [in]                         -> Coefficient of d in a * d^2 + 2 * k * d + c = 0.
  !!   c [in]                         -> Coefficient c in a * d^2 + 2 * k * d + c = 0.
  !!   surfTolCondition [in]          -> .true. if particle is within surface tolerance of the 
  !!                                     quadratic surface.
  !!   d [inout]                      -> Distance to intersection with the quadratic surface.
  !!   noIntersection [out]           -> .true. if there is no intersection with the quadratic surface.
  !!   extraCondition [in] [optional] -> Extra logical condition for which there is no intersection with
  !!                                     the quadratic surface when .true. (e.g., a = ZERO for a 
  !!                                     cylinder).
  !!
  elemental subroutine distanceHalfwidth(self, rComponent, uComponent, dMin, dMax, noIntersection)
    class(truncQuadSurface), intent(in) :: self
    real(defReal), intent(in)           :: rComponent, uComponent
    real(defReal), intent(inout)        :: dMin, dMax
    logical(defBool), intent(out)       :: noIntersection
    real(defReal)                       :: halfwidth, offsetCoord, inverseU
    
    ! Initialise noIntersection = .true. then compute halfwidth intersections.
    noIntersection = .true.
    halfwidth = self % halfwidth
    offsetCoord = rComponent - self % halfwidthOrigin

    ! Perform early check to see if the particle's direction is parallel to the slab, or if the
    ! particle is outside the slab and moving away from it along the specified dimension. In any
    ! case the particle cannot intersect the slab and we can return early.
    if (isEqual(uComponent, ZERO) .or. (offsetCoord <= -halfwidth .and. uComponent <= ZERO) .or. &
    (offsetCoord >= halfwidth .and. uComponent >= ZERO)) return

    ! Pre-compute inverse direction component and update dMin and dMax.
    inverseU = 1 / uComponent
    dMin = max(dMin, (-sign(halfwidth, uComponent) - offsetCoord) * inverseU)
    dMax = min(dMax, (sign(halfwidth, uComponent) - offsetCoord) * inverseU)

    if (dMax <= dMin * MISS_TOL) return

    ! If reached here, ray intersects the halfwidth so update noIntersection = .false.
    noIntersection = .false.

  end subroutine distanceHalfwidth

  !! Function 'isHalfspacePositive'
  !!
  !! Basic description:
  !!   Returns .true. if particle lies in the positive halfspace with respect to the surface's halfwidth.
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
  elemental function isHalfspacePositive(self, rComponent, uComponent) result(isIt)
    class(truncQuadSurface), intent(in) :: self
    real(defReal), intent(in)           :: rComponent, uComponent
    logical(defBool)                    :: isIt
    real(defReal)                       :: offsetCoord, dist

    ! Initialise isIt = .true., compute offset coordinate and distance to halfwidth. Return
    ! early if abs(dist) >= surface tolerance.
    isIt = .true.
    offsetCoord = rComponent - self % halfwidthOrigin
    dist = abs(offsetCoord) - self % halfwidth
    if (abs(dist) >= self % getSurfTol()) return
    
    ! If ray is parallel to current dimension update halfspace based on position.
    if (isEqual(uComponent, ZERO)) then
      isIt = dist >= ZERO
      return

    end if

    ! If reached here, halfspace is given by the sign of the projection of uComponent onto
    ! halfwidth's normal.
    isIt = uComponent * sign(ONE, offsetCoord) > ZERO

  end function isHalfspacePositive

  !! Function 'evaluateHalfwidth'
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
  elemental function evaluateHalfwidth(self, rComponent) result(c)
    class(truncQuadSurface), intent(in) :: self
    real(defReal), intent(in)           :: rComponent
    real(defReal)                       :: c

    c = abs(rComponent - self % halfwidthOrigin) - self % halfwidth

  end function evaluateHalfwidth

  !! Subroutine 'setHalfwidthBCs'
  !!
  !! Basic description:
  !!   Sets boundary conditions of the compound surface according to BCs. If BCs is absent, initialises
  !!   the boundary conditions to VACUUM_BC.
  !!
  !! Arguments:
  !!   nRequired [in]      -> Number of required boundary conditions for the compound surface.
  !!   BCs [in] [optional] -> Integer array containing the boundary condition flags. If absent boundary
  !!                          conditions are set to VACUUM_BC by default.
  !!
  !! Errors:
  !!   - fatalError if number of required boundary conditions does not match number of supplied boundary
  !!     conditions.
  !!   - fatalError if any boundary condition flags are unrecognised.
  !!   - fatalError if periodic boundary conditions are not applied on opposite surfaces.
  !!
  subroutine setHalfwidthBCs(self, BCs)
    class(truncQuadSurface), intent(inout)      :: self
    integer(shortInt), dimension(:), intent(in) :: BCs
    integer(shortInt)                           :: nBCs, i
    character(*), parameter                     :: here = 'setHalfwidthBCs (truncQuadSurface_inter.f90)'

    nBCs = size(BCs)
    if(nBCs < 2) call fatalError(Here, 'Boundary conditions must have size 2. Has: '//numToChar(nBCs)//'.')

    ! Verify BC flags.
    do i = 1, 2
      select case(BCs(i))
        case (VACUUM_BC, REFLECTIVE_BC, PERIODIC_BC)
          ! Do nothing, pass
        case default
          call fatalError(here, 'Unrecognised BC: '//numToChar(BCs(i))//' in position: '//numToChar(i)//'.')

      end select

    end do

    ! Verify periodic BCs.
    if ((BCs(1) == PERIODIC_BC) .neqv. (BCs(2) == PERIODIC_BC)) &
    call fatalError(here, 'Periodic BCs need to be applied on opposite surfaces.')

    ! Load BC flags.
    self % BCs = BCs(1:2)

  end subroutine setHalfwidthBCs

  !! Subroutine 'killHalfwidth'
  !!
  !! Basic description:
  !!   Returns compound surface to an unitialised state.
  !!
  elemental subroutine killHalfwidth(self)
    class(truncQuadSurface), intent(inout) :: self

    self % halfwidth = ZERO
    self % halfwidthOrigin = ZERO
    self % BCs = VACUUM_BC

  end subroutine killHalfwidth

  !! Subroutine 'setHalfwidth'
  !!
  !! Basic description:
  !!   Sets the halfwidth of the surface.
  !!
  !! Arguments:
  !!   halfwidths [in] -> Array of halfwidths.
  !!
  !! Errors:
  !!   - fatalError if the number of supplied halfwidths does not match the number of required
  !!     halfwidths.
  !!   - fatalError if any of the supplied halfwidths are negative.
  !!
  subroutine setHalfwidth(self, halfwidth, halfwidthOrigin)
    class(truncQuadSurface), intent(inout) :: self
    real(defReal), intent(in)              :: halfwidth, halfwidthOrigin
    character(*), parameter                 :: here = 'setHalfwidth (truncQuadSurface_inter.f90)'

    if (halfwidth <= ZERO) call fatalError(here, 'Halfwidth must be +ve.')
    self % halfwidth = halfwidth
    self % halfwidthOrigin = halfwidthOrigin

  end subroutine setHalfwidth

  !! Subroutine 'explicitHalfwidthBCs'
  !!
  !! Basic description:
  !!   Explicitly applies boundary conditions on particle position and direction for each
  !!   halfwidth of the compound surface along specified dimensions.
  !!
  !! Arguments:
  !!   dimensions [in]       -> Dimensions along which to apply boundary conditions.
  !!   originComponents [in] -> Components of the compound surface's origin along the
  !!                            specified dimensions.
  !!   r [inout]             -> Particle position.
  !!   u [inout]             -> Particle direction.
  !!
  pure subroutine explicitHalfwidthBCs(self, dim, r, u)
    class(truncQuadSurface), intent(in)        :: self
    integer(shortInt), intent(in)              :: dim
    real(defReal), dimension(3), intent(inout) :: r, u
    integer(shortInt)                          :: bc
    real(defReal)                              :: offsetCoord, halfwidth

    ! Retrieve halfwidth and offset coordinate with respect to box centre
    ! for the current dimension.
    halfwidth = self % halfwidth
    offsetCoord = r(dim) - self % halfwidthOrigin

    ! If particle is well inside the domain return early.
    if (abs(offsetCoord) <= halfwidth - self % getSurfTol()) return

    ! Retrieve appropriate BC based on whether offsetCoord < ZERO and apply it.
    if (offsetCoord < ZERO) then
      bc = self % BCs(1)

    else
      bc = self % BCs(2)

    end if
    if (bc == REFLECTIVE_BC) then
      u(dim) = -u(dim)
      return

    end if
    if (bc == PERIODIC_BC) r(dim) = r(dim) - TWO * sign(halfwidth, offsetCoord)

  end subroutine explicitHalfwidthBCs

  !! Subroutine 'transformHalfwidthBCs'
  !!
  !! Basic description:
  !!   Transforms particle position and direction using boundary conditions.
  !!
  !! Arguments:
  !!   dimensions [in]       -> Dimensions along which to apply boundary conditions.
  !!   originComponents [in] -> Components of the compound surface's origin along the
  !!                            specified dimensions.
  !!   r [inout]             -> Particle position.
  !!   u [inout]             -> Particle direction.
  !!
  pure subroutine transformHalfwidthBCs(self, dim, r, u)
    class(truncQuadSurface), intent(in)        :: self
    integer(shortInt), intent(in)              :: dim
    real(defReal), dimension(3), intent(inout) :: r, u
    integer(shortInt)                          :: nTransforms, i, bc
    real(defReal)                              :: offsetCoord, halfwidthOrigin, halfwidth

    ! Calculate number of transforms in each direction. Reduce halfwidths by surface
    ! tolerance to catch particles at the surface.
    halfwidth = self % halfwidth
    halfwidthOrigin = self % halfwidthOrigin
    nTransforms = ceiling(abs(r(dim) - halfwidthOrigin) / (halfwidth - self % getSurfTol())) / 2

    ! Apply number of transforms.
    do i = 1, nTransforms
      ! Offset coordinates with respect to halfwidth origin.
      offsetCoord = r(dim) - halfwidthOrigin

      ! Retrieve appropriate BC based on whether offsetCoord < ZERO and apply it.
      if (offsetCoord < ZERO) then
        bc = self % BCs(1)

      else
        bc = self % BCs(2)

      end if
      if (bc == REFLECTIVE_BC) then
        ! Perform reflection based on distance to the closest plane to offsetCoord.
        r(dim) = -r(dim) + TWO * (sign(halfwidth, offsetCoord) + halfwidthOrigin)
        u(dim) = -u(dim)
        cycle

      end if

      ! If periodic BC perform translation.
      if (bc == PERIODIC_BC) r(dim) = r(dim) - TWO * sign(halfwidth, offsetCoord)

    end do

  end subroutine transformHalfwidthBCs
  
end module truncQuadSurface_inter