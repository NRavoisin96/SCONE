module cone_class

  use numPrecision
  use universalVariables, only : SURF_TOL, INF, X_AXIS, Y_AXIS, Z_AXIS
  use genericProcedures,  only : fatalError, numToChar, isEqual, computeQuadraticSolutions
  use dictionary_class,   only : dictionary
  use quadSurface_inter,  only : quadSurface
  use surface_inter,      only : kill_super => kill

  implicit none
  private

  !!
  !! Axis aligned cone, truncated at hMin and hMax
  !!
  !! Surface description:
  !! F(r) = max{ (ri - vi)^2 + (rj - vj)^2 - t^2 * (rk - vk)^2;
  !!             (rk - vk) - hMax;
  !!            -(rk - vk) + hMin; }
  !!
  !! Where i,j,k are x,y & z axis, and the cone is aligned along the k axis; t is
  !! the tangent of the cone opening angle (defined as the angle between the axis
  !! and the cone surface). The point V(vx, vy, vz) is the cone vertex.
  !!
  !! The cone is truncated by two faces given by their coordinates along the axis
  !! of the cone (hMin and hMax) with the vertex being 0. The sign of hMin and hMax
  !! determines the orientation of the cone.
  !! NOTE: Both entries must have the same sign (single truncated cone), and
  !! hMin < hMax must be true.
  !!
  !! The opening angle of the cone must be provided in degrees, and be in the
  !! range 0 - 90 degrees (extremes excluded)
  !!
  !! Three different types are available
  !!   xCone -> aligned with X-axis
  !!   yCone -> aligned with Y-axis
  !!   zCone -> aligned with Z-axis
  !!
  !! Surface tolerance varies with the position along the cone:
  !!   -SURF_TOL for the cone bases.
  !!   -2 * R(rk) * SURF_TOL for the conic surface, where R(rk) = (rk - vk) * tan(theta)
  !!
  !! Sample dictionary input:
  !!   cone { type xCone; // could be yCone or zCone as well
  !!          id 3;
  !!          vertex (0.0 0.0 0.0);
  !!          angle 45;
  !!          hMin 0.0;
  !!          hMax 10.0; }
  !!
  !! Private Members:
  !!   axis       -> Index of an alignment axis in {X_AXIS, Y_AXIS, Z_AXIS}
  !!   planes     -> Indices of axis in planes of cone {X_AXIS, Y_AXIS, Z_AXIS}\{axis}
  !!   vertex     -> Location of the vertex of the cone
  !!   tan        -> Tangent of the opening angle
  !!   tanSquared -> Square of the tangent of the opening angle
  !!   hMin       -> Cone lower boundary
  !!   hMax       -> Cone upper boundary
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(quadSurface) :: cone
    private
    integer(shortInt)               :: axis = 0
    integer(shortInt), dimension(2) :: planes = 0
    real(defReal), dimension(3)     :: vertex = ZERO
    real(defReal)                   :: tan = ZERO, tanSquared = ZERO, inverseCos = ZERO, hMin = ZERO, hMax = ZERO
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: distance
    procedure :: entersPositiveHalfspace
    procedure :: evaluate
    procedure :: isWithinSurfTol
  end type cone

contains

  !!
  !! Initialise cone from a dictionary
  !!
  !! See surface_inter for more details
  !!
  subroutine init(self, dict)
    class(cone), intent(inout)               :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id, N, axis
    integer(shortInt), dimension(2)          :: planes
    real(defReal), dimension(:), allocatable :: vertex
    character(nameLen)                       :: type
    real(defReal)                            :: angle, tangent, hMin, hMax, maxRadius
    real(defReal), dimension(6)              :: boundingBox
    character(100), parameter                :: Here = 'init (cone_class.f90)'

    ! Load id.
    call dict % get(id, 'id')
    call self % setId(id)

    ! Load vertex.
    call dict % get(vertex, 'vertex')
    N = size(vertex)
    if (N /= 3) call fatalError(Here,'Vertex must have size 3. Has: '//numToChar(N)//'.')
    self % vertex = vertex

    ! Load angle.
    call dict % get(angle, 'angle')
    if (angle <= ZERO .or. angle >= 90.0_defReal) call fatalError(Here, &
    'Cone opening angle must be in the range 0-90 degrees (extremes excluded). Is: '//numToChar(angle)//'.')
    
    ! Convert angle to radians and set tan, tanSquared and inverseCos.
    angle = angle * PI / 180.0_defReal
    tangent = tan(angle)
    self % tan = tangent
    self % tanSquared = tangent * tangent
    self % inverseCos = ONE / cos(angle)

    ! Load hMin and hMax.
    call dict % get(hMin, 'hMin')
    call dict % get(hMax, 'hMax')
    if (hMin >= hMax) call fatalError(Here, 'hMin must be strictly smaller than hMax.')
    if (hMin * hMax < ZERO) call fatalError(Here, 'hMin and hMax have different signs.')
    self % hMin = hMin
    self % hMax = hMax

    ! Load type. Set axis and planes accordingly.
    call dict % get(type, 'type')
    select case(type)
      case('xCone')
        axis = X_AXIS
        planes = [Y_AXIS, Z_AXIS]

      case('yCone')
        axis = Y_AXIS
        planes = [X_AXIS, Z_AXIS]

      case('zCone')
        axis = Z_AXIS
        planes = [X_AXIS, Y_AXIS]

      case default
        call fatalError(Here, 'Unknown type of cone: '//type)

    end select
    self % axis = axis
    self % planes = planes

    ! Set bounding box. Set along cone axis first.
    boundingBox(axis) = vertex(axis) + hMin
    boundingBox(axis + 3) = vertex(axis) + hMax

    ! Compute maximum cone radius and set bounding box along the planes of the bases.
    maxRadius = max(abs(hMin), abs(hMax)) * tangent
    boundingBox(planes) = vertex(planes) - maxRadius
    boundingBox(planes + 3) = vertex(planes) + maxRadius
    call self % setBoundingBox(boundingBox)
    call self % setType(type)

  end subroutine init

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(cone), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: offsetCoords
    real(defReal), dimension(2)             :: offsetCoordsPlanes
    real(defReal)                           :: offsetCoordAxis, radius, c, cMin, cMax

    ! Offset coordinates with respect to cone's vertex and retrieve offset coordinates 
    ! components along the cone's axis and planes.
    offsetCoords = r - self % vertex
    offsetCoordAxis = offsetCoords(self % axis)
    offsetCoordsPlanes = offsetCoords(self % planes)
    
    ! Compute radius along the axis of the cone and evaluate c along the cone's axis.
    radius = self % tan * offsetCoordAxis
    if (isEqual(radius, ZERO)) then
      c = dot_product(offsetCoordsPlanes, offsetCoordsPlanes)

    else
      c = HALF * (dot_product(offsetCoordsPlanes, offsetCoordsPlanes) / radius - radius)

    end if

    ! Evaluate c for the two bases of the cone and return the overall maximum.
    cMin = -offsetCoordAxis + self % hMin
    cMax = offsetCoordAxis - self % hMax
    c = max(c, cMin, cMax)

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! See surface_inter for details
  !!
  !! Solves quadratic intersection equation
  !!   a*d^2 + 2b*d + c = 0
  !!   c = cone surface expression
  !!   b = (r1 - v1)u1 + (r2 - v2)u2 - t^2(r3 - v3)u3
  !!   a = u1^2 + u2^2 - t^2 * u3^2
  !!
  pure function distance(self, r, u) result(d)
    class(cone), intent(in)                  :: self
    real(defReal), dimension(3), intent(in)  :: r, u
    real(defReal), dimension(3)              :: offsetCoords
    real(defReal), dimension(2)              :: offsetCoordsPlanes, uPlanes
    real(defReal)                            :: d, offsetCoordAxis, uAxis, uAxisScaled, radius, &
                                                a, b, c, delta, hMin, hMax, low, high, inverseUAxis, &
                                                dMin, dMax, dBase, bound, surfTol
    real(defReal), dimension(:), allocatable :: solutions
    integer(shortInt), dimension(2)          :: planes
    integer(shortInt)                        :: axis, nSolutions, i
    real(defReal), parameter :: FP_MISS_TOL = ONE + 10.0_defReal * epsilon(ONE)

    ! Initialise d = INF. Offset coordinates with respect to cone vertex and retrieve offset 
    ! coordinates and direction components along the cone's planes and axis.
    d = INF
    offsetCoords = r - self % vertex
    planes = self % planes
    axis = self % axis
    offsetCoordsPlanes = offsetCoords(planes)
    offsetCoordAxis = offsetCoords(axis)
    uPlanes = u(planes)
    uAxis = u(axis)
    hMin = self % hMin
    hMax = self % hMax

    ! Check if particle is outside the bases of the cone and moves parallel to them, since we can
    ! return early in this case. Note, here check for perfect equality since if the particle is
    ! very close to the surface and does not move exactly parallel to it, it may still intersect.
    if ((offsetCoordAxis < hMin .or. offsetCoordAxis > hMax) .and. uAxis == ZERO) return

    ! Pre-compute uAxisScaled = uAxis / cos(angle) and cone radius, then compute a, b and c.
    uAxisScaled = uAxis * self % inverseCos
    radius = self % tan * offsetCoordAxis
    a = ONE - uAxisScaled * uAxisScaled
    b = dot_product(offsetCoordsPlanes , uPlanes) - self % tanSquared * offsetCoordAxis * uAxis
    c = dot_product(offsetCoordsPlanes, offsetCoordsPlanes) - radius * radius

    ! Compute delta (technically, delta / 4).
    delta = b * b - a * c

    ! If delta < ZERO, the solutions are complex. There is no intersection with the conic surface
    ! and there can be no intersection with either base of the cone so we can return early.
    if (delta < ZERO) return

    ! Initialise intersections with the planes of the cone's bases dMin = ZERO and dMax = INF.
    dMin = ZERO
    dMax = INF

    ! Only update dMin and dMax if uAxis /= ZERO.
    if (.not. isEqual(uAxis, ZERO)) then
      ! Pre-compute 1 / uAxis and update dMin and dMax.
      inverseUAxis = ONE / uAxis
      low = (hMin - offsetCoordAxis) * inverseUAxis
      high = (hMax - offsetCoordAxis) * inverseUAxis
      dMin = min(low, high)
      dMax = max(low, high)

    end if

    ! If dMin <= ZERO and dMax <= ZERO, the ray is moving away from both planes of the cone's bases
    ! and we can return early.
    if (dMin <= ZERO .and. dMax <= ZERO) return

    ! Compute distance to intersection of closest plane of the cone's bases.
    dBase = dMin
    if (dMin <= ZERO) dBase = dMax

    ! Compute distances to intersection with the cone surface. If there are no solutions then the ray 
    ! is parallel to the cone opening and fully contained in its surface. If there is only one solution 
    ! then the ray is parallel to the cone opening. If this solution is negative and the particle is 
    ! inside the cone (c < ZERO) then the ray can only intersect with one of the cone bases and we can 
    ! return early.
    solutions = computeQuadraticSolutions(a, b, c, delta)
    nSolutions = size(solutions)
    if (nSolutions == 0 .or. (nSolutions == 1 .and. solutions(1) <= ZERO .and. c < ZERO)) then
      d = dBase
      return

    end if

    ! Pre-compute bound and loop over all solutions retrieved.
    bound = sign(INF, uAxis * hMin)
    do i = 1, nSolutions
      ! Retrieve current solution and cap it to signed INF if the intersection is in the wrong hemicone.
      if ((offsetCoordAxis + uAxis * solutions(i)) * hMin < ZERO) solutions(i) = bound

    end do

    ! Update dMin and dMax.
    dMin = max(minval(solutions), dMin)
    dMax = min(maxval(solutions), dMax)

    ! If dMin > dMax, there is no intersection and we can return.
    if (dMin > dMax) return

    ! If reached here, update d = dMin.
    d = dMin

    ! Check if particle is on surface or already inside the cone and update d = tMax if yes.
    surfTol = self % getSurfTol()
    if ((self % evaluate(r) < surfTol .and. abs(dMax) >= abs(dMin)) .or. dMin <= ZERO) d = dMax

    ! Cap distance to INF if particle is within surfTol to another surface (eg, particle
    ! previously entered the cone very close to the intersection between a base and the surface) 
    ! or if d > INF.
    if (d < surfTol .or. d > INF) d = INF

  end function distance

  !!
  !! Returns TRUE if particle is going into +ve halfspace
  !!
  !! See surface_inter for details
  !!
  pure function entersPositiveHalfspace(self, r, u) result(isHalfspacePositive)
    class(cone), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r, u
    logical(defBool)                        :: isHalfspacePositive
    real(defReal), dimension(3)             :: offsetCoords, normal
    integer(shortInt), dimension(2)         :: planes
    integer(shortInt)                       :: axis
    real(defReal), dimension(2)             :: offsetCoordsPlanes, uPlanes
    real(defReal)                           :: offsetCoordAxis, uAxis, surfTol, radius, &
                                               cMin, cMax, cCone, proj

    ! Retrieve offset coordinates with respect to cone vertex and surface tolerance.
    offsetCoords = r - self % vertex
    axis = self % axis
    offsetCoordAxis = offsetCoords(axis)
    uAxis = u(axis)
    surfTol = self % getSurfTol()
    
    ! Check cone bases first and return as soon as halfspace is positive.
    cMin = self % hMin - offsetCoordAxis
    if (abs(cMin) < surfTol) then
      if (isEqual(uAxis, ZERO)) then
        isHalfspacePositive = cMin >= ZERO

      else
        isHalfspacePositive = uAxis < ZERO

      end if

    end if
    if (isHalfspacePositive) return

    cMax = offsetCoordAxis - self % hMax
    if (abs(cMax) < surfTol) then
      if (isEqual(uAxis, ZERO)) then
        isHalfspacePositive = cMax >= ZERO

      else
        isHalfspacePositive = uAxis > ZERO

      end if

    end if
    if (isHalfspacePositive) return

    ! If reached here check cone surface.
    planes = self % planes
    offsetCoordsPlanes = offsetCoords(planes)
    uPlanes = u(planes)
    radius = self % tan * offsetCoordAxis

    ! First check if particle location along the cone's axis is exactly (to the bit) on the cone vertex.
    if (radius == ZERO) then
      ! Compute cCone. If particle is exactly on the cone vertex the cone normal is undefined. Use particle
      ! direction to determine halfspace in this case and return.
      cCone = dot_product(offsetCoordsPlanes, offsetCoordsPlanes)
      if (cCone == ZERO) isHalfspacePositive = dot_product(uPlanes, uPlanes) - self % tanSquared * uAxis * uAxis >= ZERO
      return

    end if

    ! If particle is not exactly on the cone vertex compute cone normal and project direction onto it. Use
    ! particle location to determine halfspace if projection is ZERO.
    cCone = HALF * (dot_product(offsetCoordsPlanes, offsetCoordsPlanes) / radius - radius)
    if (abs(cCone) < surfTol) then
      normal(axis) = -self % tanSquared * offsetCoordAxis
      normal(planes) = offsetCoordsPlanes
      normal = normal / norm2(normal)
      proj = dot_product(normal, u)

      if (isEqual(proj, ZERO)) then
        isHalfspacePositive = cCone >= ZERO
        
      else
        isHalfspacePositive = proj > ZERO

      end if

    end if

  end function entersPositiveHalfspace

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(cone), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % axis = 0
    self % planes = 0
    self % vertex = ZERO
    self % tan = ZERO
    self % tanSquared = ZERO
    self % inverseCos = ZERO
    self % hMin = ZERO
    self % hMax = ZERO

  end subroutine kill

  !! Function 'isWithinSurfTol'
  !!
  !! Basic description:
  !!   Returns .true. if particle is within surface tolerance of the cone.
  !!
  !! Arguments:
  !!   r [in] -> Location of the particle.
  !!
  !! Result:
  !!   isIt   -> .true. if particle is within surface tolerance of the cone.
  !!
  pure function isWithinSurfTol(self, r) result(isIt)
    class(cone), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r
    logical(defBool)                        :: isIt
    real(defReal), dimension(3)             :: offsetCoords
    integer(shortInt), dimension(2)         :: planes
    integer(shortInt)                       :: axis
    real(defReal)                           :: offsetCoordAxis, radius, c, cMin, cMax
    real(defReal), dimension(2)             :: offsetCoordsPlanes


    offsetCoords = r - self % vertex
    planes = self % planes
    axis = self % axis
    offsetCoordAxis = offsetCoords(axis)
    offsetCoordsPlanes = offsetCoords(planes)
    
    ! If particle's offset coordinate component along the cone's axis is exactly (to the bit)
    ! ZERO, then particle is within surface tolerance if it is exactly on the cone's vertex.
    radius = self % tan * offsetCoordAxis
    cMin = -offsetCoordAxis + self % hMin
    cMax = offsetCoordAxis - self % hMax
    if (radius == ZERO) then
      isIt = dot_product(offsetCoordsPlanes, offsetCoordsPlanes) == ZERO .and. max(cMin, cMax) < self % getSurfTol()
      return

    end if
    c = HALF * (dot_product(offsetCoordsPlanes, offsetCoordsPlanes) / radius - radius)
    isIt = max(c, cMin, cMax) < self % getSurfTol()

  end function isWithinSurfTol

end module cone_class