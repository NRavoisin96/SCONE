module cylinder_class

  use numPrecision
  use universalVariables, only : SURF_TOL, INF, X_AXIS, Y_AXIS, Z_AXIS
  use genericProcedures,  only : fatalError, numToChar, isEqual
  use dictionary_class,   only : dictionary
  use quadSurface_inter,  only : quadSurface
  use surface_inter,      only : kill_super => kill

  implicit none
  private

  !!
  !! Axis aligned cylinder
  !!
  !! F(r) = (ri - oi)^2 + (rj - oj)^2 - R^2 = 0
  !!
  !! Where i,j (i /= j) can be any of x,y & z axis
  !!
  !! Three diffrent types are avaliable
  !!   xCylinder -> aligned with X-axis
  !!   yCylinder -> aligned with Y-axis
  !!   zCylinder -> aligned with Z-axis
  !!
  !! Surface tolerance: 2 * R * SURF_TOL
  !!
  !! Sample dictionary input:
  !!   cyl { type xCylinder; // could be yCylinder or zCylinder as well
  !!         id 3;
  !!         origin (0.0 0.0 0.0);
  !!         radius 7.34; }
  !!
  !! Private Members:
  !!   axis          -> Index of an alignment axis in {X_AXIS, Y_AXIS, Z_AXIS}
  !!   planes        -> Indexes of axis in planes of cylinder {X_AXIS, Y_AXIS, Z_AXIS}\{axis}
  !!   origin        -> Location of the middle of the cylinder.
  !!   radius        -> Radius of the cylinder.
  !!   radiusSquared -> Radius of the cylinder squared.
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(quadSurface) :: cylinder
    private
    integer(shortInt)               :: axis = 0
    integer(shortInt), dimension(2) :: planes = 0
    real(defReal)                   :: radius = ZERO, radiusSquared = ZERO
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: evaluate
    procedure :: distance
    procedure :: entersPositiveHalfspace
    procedure :: kill

    ! Local procedures.
    procedure :: build

  end type cylinder

contains
  !!
  !! Build cylinder from components
  !!
  !! Args:
  !!   id [in] -> Surface ID
  !!   type [in] -> Cylinder type {'xCylinder', 'yCylinder' or 'zCylinder'}
  !!   origin [in] -> Cylinder origin
  !!   radius [in] -> Cylinder radius
  !!
  !! Errors:
  !!   fatalError if id or radius are -ve
  !!
  subroutine build(self, id, type, origin, radius)
    class(cylinder), intent(inout)          :: self
    integer(shortInt), intent(in)           :: id
    character(*), intent(in)                :: type
    real(defReal), dimension(3), intent(in) :: origin
    real(defReal), intent(in)               :: radius
    integer(shortInt)                       :: axis
    integer(shortInt), dimension(2)         :: planes
    real(defReal), dimension(6)             :: boundingBox
    character(100), parameter               :: Here = 'build (cylinder_class.f90)'

    ! Check values.
    if (id < 1) call fatalError(Here, 'Invalid surface id provided. Id must be > 1.')
    if (radius <= ZERO) call fatalError(Here, 'Radius of cylinder must be +ve. Is: '//numToChar(radius)//'.')

    ! Select type of cylinder
    select case(type)
      case('xCylinder')
        axis = X_AXIS
        planes = [Y_AXIS, Z_AXIS]

      case('yCylinder')
        axis = Y_AXIS
        planes = [X_AXIS, Z_AXIS]

      case('zCylinder')
        axis = Z_AXIS
        planes = [X_AXIS, Y_AXIS]

      case default
        call fatalError(Here, 'Unknown type of cylinder: '//type//'.')

    end select

    ! Load data.
    self % axis = axis
    self % planes = planes
    self % radius = radius
    self % radiusSquared = radius * radius
    call self % setOrigin(origin(planes), 2)

    ! Set id, surface tolerance and bounding box.
    call self % setId(id)
    call self % setSurfTol(TWO * radius * SURF_TOL)
    boundingBox(planes) = origin(planes) - radius
    boundingBox(planes + 3) = origin(planes) + radius
    boundingBox(axis) = -INF
    boundingBox(axis + 3) = INF
    call self % setBoundingBox(boundingBox)
    call self % setType(type)

  end subroutine build

  !!
  !! Initialise cylinder from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!   fatalError if radius or id < 0.
  !!
  subroutine init(self, dict)
    class(cylinder), intent(inout)           :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id
    character(nameLen)                       :: type
    real(defReal), dimension(:), allocatable :: origin
    real(defReal)                            :: radius
    character(100), parameter                :: Here = 'init (cylinder_class.f90)'

    ! Get from dictionary
    call dict % get(id, 'id')
    call dict % get(radius, 'radius')
    call dict % get(origin, 'origin')
    call dict % get(type, 'type')

    ! Build cylinder.
    call self % build(id, type, origin, radius)

  end subroutine init

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(cylinder), intent(in)             :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c
    real(defReal), dimension(2)             :: offsetCoords

    ! Offset coordinates with respect to cylinder's origin and compute c.
    offsetCoords = r(self % planes) - self % getOrigin()
    c = dot_product(offsetCoords, offsetCoords) - self % radiusSquared

  end function evaluate

  !! Function 'distance'
  !!
  !! Basic description:
  !!   Returns the distance to the surface of the cylinder.
  !!
  !! Detailed description:
  !!   Consider an infinite cylinder parallel to the z-axis (the reasoning is identical for cylinders
  !!   parallel to the other two axes); a point lying on the surface of this cylinder satisfies:
  !!
  !!   (x - x0)² + (y - y0)² - R² = 0
  !!
  !!   Now consider a point p(d) = [p_x(d), p_y(d), p_z(d)] along a ray of direction u = [u_x, u_y, u_z].
  !!   Its parametric equation can be written as:
  !!
  !!   p(d) = r + d * u
  !!
  !!   Substituting this into the equation for the cylinder, we get:
  !!
  !!   (r_x + d * u_x - x0)² + (r_y + d * u_y - y0)² - R² = 0
  !!=> (d * u_x + (r_x - x0))² + (d * u_y + (r_y - y0))² - R² = 0
  !!=> (u_x² + u_y²) * d² + 2 * ((r_x - x0) * u_x + (r_y - y0) * u_y) * d + (r_x - x0)² + (r_y - y0)² - R² = 0
  !!
  !!   Now, define:
  !!
  !!   a = u_x² + u_y²
  !!     = 1 - u_z² (since norm(u) = 1)
  !!   k = (r_x - x0) * u_x + (r_y - y0) * u_y
  !!   c = (r_x - x0)² + (r_y - y0)² - R² (note: this is the value returned by the 'evaluate' function)
  !!
  !!   Then the equation becomes:
  !!
  !!   ad² + 2kd + c = 0
  !!
  !!   which is quadratic in d. The determinant (delta) is given by:
  !!
  !!   delta = 4k² - 4ac
  !!         = 4(k² - ac)
  !!
  !!   And the two solutions for d are therefore:
  !!
  !!   d_1 = (-2k + 2sqrt(k² - ac)) / (2a)
  !!       = (-k + sqrt(k² - ac)) / a
  !!   d_2 = (-2k - 2sqrt(k² - ac)) / (2a)
  !!       = -(k + sqrt(k² - ac)) / a
  !!
  !! Arguments:
  !!   r [in] -> Position of the particle.
  !!   u [in] -> Direction of the particle.
  !!
  pure function distance(self, r, u) result(d)
    class(cylinder), intent(in)             :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal)                           :: d, uAxis, a, c, k, delta
    integer(shortInt), dimension(2)         :: planes

    ! Initialise d = INF, retrieve planes and axial direction component and calculate a, k and c.
    d = INF
    planes = self % planes
    uAxis = u(self % axis)
    a = ONE - uAxis * uAxis
    k = dot_product(r(planes) - self % getOrigin(), u(planes))
    c = self % evaluate(r)

    ! Compute delta (technically, delta / 4).
    delta = k * k - a * c

    ! If delta < ZERO, the solutions are complex. If a = ZERO, the ray is parallel to the cylinder's
    ! axis. In any case there is no intersection and we can return early.
    if (delta < ZERO .or. isEqual(a, ZERO)) return

    ! Check if particle is within surface tolerance of the cylinder.
    if (abs(c) < self % getSurfTol()) then
      ! Update d only if k < ZERO (k >= ZERO corresponds to the particle moving away from the cylinder). 
      ! Choose maximum distance and return.
      if (k < ZERO) d = (-k + sqrt(delta)) / a
      return

    end if

    ! If reached here, update d depending on the sign of c and cap distance at infinity.
    d = -(k + sign(sqrt(delta), c)) / a
    if (d <= ZERO .or. d > INF) d = INF

  end function distance

  !!
  !! Returns .true. if particle is going into +ve halfspace.
  !!
  !! See surface_inter for details
  !!
  pure function entersPositiveHalfspace(self, r, u) result(isHalfspacePositive)
    class(cylinder), intent(in)             :: self
    real(defReal), dimension(3), intent(in) :: r, u
    integer(shortInt), dimension(2)         :: planes
    logical(defBool)                        :: isHalfspacePositive

    planes = self % planes
    isHalfspacePositive = dot_product(r(planes) - self % getOrigin() , u(planes)) >= ZERO

  end function entersPositiveHalfspace

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(cylinder), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % axis = 0
    self % planes = 0
    self % radius = ZERO
    self % radiusSquared = ZERO

  end subroutine kill

end module cylinder_class
