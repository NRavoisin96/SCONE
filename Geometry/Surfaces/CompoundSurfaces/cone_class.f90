module cone_class

  use numPrecision
  use universalVariables,    only : SURF_TOL, HALF, INF, X_AXIS, Y_AXIS, Z_AXIS
  use genericProcedures,     only : fatalError, numToChar, isEqual, computeQuadraticSolutions
  use dictionary_class,      only : dictionary
  use compoundSurface_inter, only : compoundSurface
  use surface_inter,         only : kill_super => kill

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
  type, public, extends(compoundSurface) :: cone
    private
    integer(shortInt)                    :: axis = 0
    integer(shortInt), dimension(2)      :: planes = 0
    real(defReal), dimension(3)          :: vertex = ZERO
    real(defReal)                        :: tan = ZERO, tanSquared = ZERO, inverseCos = ZERO, hMin = ZERO, hMax = ZERO
  contains
    ! Superclass procedures
    procedure                            :: init
    procedure                            :: kill
    procedure                            :: distance
    procedure                            :: entersPositiveHalfspace
    procedure                            :: evaluate
    procedure                            :: isWithinSurfTol
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

    ! Compute origin along the cone axial dimension and halfwidth.
    call self % setOrigin([vertex(axis) + HALF * (hMin + hMax)], 1)
    call self % setHalfwidths([HALF * (hMax - hMin)], 1)

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
    real(defReal)                           :: offsetCoordAxis, radius, c

    ! Offset coordinates with respect to cone's vertex and retrieve offset coordinates 
    ! components along the cone's axis and planes.
    offsetCoords = r - self % vertex
    offsetCoordAxis = offsetCoords(self % axis)
    offsetCoordsPlanes = offsetCoords(self % planes)
    
    ! Compute radius along the axis of the cone and evaluate c along the cone's axis.
    radius = self % tan * offsetCoordAxis
    if (radius == ZERO) then
      c = dot_product(offsetCoordsPlanes, offsetCoordsPlanes)

    else
      c = HALF * (dot_product(offsetCoordsPlanes, offsetCoordsPlanes) / radius - radius)

    end if

    ! Evaluate c for the two bases of the cone and return the overall maximum.
    c = max(c, self % evaluateCompound([r(self % axis) - self % getOrigin()]))

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
                                                a, k, c, delta, hMin, dMin, dMax, bound, tangent
    real(defReal), dimension(:), allocatable :: solutions
    integer(shortInt), dimension(2)          :: planes
    integer(shortInt)                        :: axis, nSolutions, i
    logical(defBool)                         :: cannotIntersect

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

    ! Pre-compute uAxisScaled = uAxis / cos(angle) and cone radius, then compute a, k and c.
    uAxisScaled = uAxis * self % inverseCos
    tangent = self % tan
    radius = tangent * offsetCoordAxis
    a = ONE - uAxisScaled * uAxisScaled
    k = dot_product(offsetCoordsPlanes , uPlanes) - tangent * radius * uAxis
    c = dot_product(offsetCoordsPlanes, offsetCoordsPlanes) - radius * radius

    ! Compute delta (technically, delta / 4). Return early if delta < ZERO.
    delta = k * k - a * c
    if (delta < ZERO) return

    ! Compute distances to intersection with the cone surface and initialise dMin = ZERO and
    ! dMax = INF.
    solutions = computeQuadraticSolutions(a, k, c, delta)
    nSolutions = size(solutions)
    dMin = ZERO
    dMax = INF

    ! If there are no solutions then the ray is parallel to the cone opening and fully contained in its 
    ! surface. If there is only one solution then the ray is parallel to the cone opening. If this solution 
    ! is negative and the particle is inside the cone (c < ZERO) then the ray can only intersect with one of 
    ! the cone bases so compute this intersection and return.
    if (nSolutions == 0 .or. (nSolutions == 1 .and. solutions(1) <= ZERO .and. c < ZERO)) then
      call self % distancesCompound([r(axis) - self % getOrigin()], [uAxis], dMin, dMax, cannotIntersect)
      call self % chooseDistance(dMin, dMax, self % isWithinSurfTol(r), cannotIntersect, d)
      return

    end if

    ! Pre-compute bound and loop over all solutions retrieved.
    bound = sign(INF, uAxis * hMin)
    do i = 1, nSolutions
      ! Retrieve current solution and cap it to signed INF if the intersection is in the wrong hemicone.
      if ((offsetCoordAxis + uAxis * solutions(i)) * hMin < ZERO) solutions(i) = bound

    end do

    ! Update dMin and dMax, compute distances to intersection with the cone's halfwidth and choose correct
    ! distance.
    dMin = max(minval(solutions), dMin)
    dMax = min(maxval(solutions), dMax)
    call self % distancesCompound([r(axis) - self % getOrigin()], [uAxis], dMin, dMax, cannotIntersect)
    call self % chooseDistance(dMin, dMax, self % isWithinSurfTol(r), cannotIntersect, d)

  end function distance

  !!
  !! Returns TRUE if particle is going into +ve halfspace
  !!
  !! See surface_inter for details
  !!
  pure function entersPositiveHalfspace(self, r, u) result(positiveHalfspace)
    class(cone), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r, u
    logical(defBool)                        :: positiveHalfspace
    real(defReal), dimension(3)             :: offsetCoords, normal
    integer(shortInt), dimension(2)         :: planes
    integer(shortInt)                       :: axis
    real(defReal), dimension(2)             :: offsetCoordsPlanes, uPlanes
    real(defReal)                           :: offsetCoordAxis, uAxis, radius, c, proj

    ! Retrieve offset coordinates with respect to cone vertex and surface tolerance.
    offsetCoords = r - self % vertex
    axis = self % axis
    offsetCoordAxis = offsetCoords(axis)
    uAxis = u(axis)
    
    ! Check cone halfwidth first and return if halfspace is positive.
    positiveHalfspace = self % isHalfspacePositive([r(axis) - self % getOrigin()], [uAxis])
    if (positiveHalfspace) return

    ! Check cone surface.
    planes = self % planes
    offsetCoordsPlanes = offsetCoords(planes)
    uPlanes = u(planes)
    radius = self % tan * offsetCoordAxis

    ! First check if particle location along the cone's axis is exactly (to the bit) on the cone vertex.
    if (radius == ZERO) then
      ! Compute c = F(r). If particle is exactly on the cone vertex the cone normal is undefined. Use particle
      ! direction to determine halfspace in this case and return.
      c = dot_product(offsetCoordsPlanes, offsetCoordsPlanes)
      if (c == ZERO) positiveHalfspace = dot_product(uPlanes, uPlanes) - self % tanSquared * uAxis * uAxis >= ZERO
      return

    end if

    ! If particle is not exactly on the cone vertex compute cone normal and project direction onto it. Use
    ! particle location to determine halfspace if projection is ZERO.
    c = HALF * (dot_product(offsetCoordsPlanes, offsetCoordsPlanes) / radius - radius)
    if (abs(c) < self % getSurfTol()) then
      normal(axis) = -self % tanSquared * offsetCoordAxis
      normal(planes) = offsetCoordsPlanes
      normal = normal / norm2(normal)
      proj = dot_product(normal, u)

      if (isEqual(proj, ZERO)) then
        positiveHalfspace = c >= ZERO
        
      else
        positiveHalfspace = proj > ZERO

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

    ! Compound surface.
    call self % killCompound()

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
    real(defReal)                           :: offsetCoordAxis, radius, c, cHalfwidth
    real(defReal), dimension(2)             :: offsetCoordsPlanes

    offsetCoords = r - self % vertex
    planes = self % planes
    axis = self % axis
    offsetCoordAxis = offsetCoords(axis)
    offsetCoordsPlanes = offsetCoords(planes)
    
    ! If particle's offset coordinate component along the cone's axis is exactly (to the bit)
    ! ZERO, then particle is within surface tolerance if it is exactly on the cone's vertex.
    radius = self % tan * offsetCoordAxis
    cHalfwidth = self % evaluateCompound([r(axis) - self % getOrigin()])
    if (radius == ZERO) then
      isIt = dot_product(offsetCoordsPlanes, offsetCoordsPlanes) == ZERO .and. cHalfwidth < self % getSurfTol()
      return

    end if
    c = HALF * (dot_product(offsetCoordsPlanes, offsetCoordsPlanes) / radius - radius)
    isIt = max(c, cHalfwidth) < self % getSurfTol()

  end function isWithinSurfTol

end module cone_class