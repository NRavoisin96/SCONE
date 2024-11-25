module sphere_class

  use numPrecision
  use universalVariables, only : INF, SURF_TOL
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use surface_inter,      only : kill_super => kill
  use quadSurface_inter,  only : quadSurface

  implicit none
  private

  !!
  !! Sphere surface
  !!
  !! F(r) = (r1-x0)^2 + (r2-y0)^2 + (r3-z0)^2  - R^2
  !!
  !! Surface tolerance: 2 * R * SURF_TOL
  !!
  !! Sample dictionary input:
  !!   sph { type sphere;
  !!         id 17;
  !!         origin (-1.0 -2.0 0.0);
  !!         radius 5.0;
  !!       }
  !!
  !! Private Members:
  !!   r      -> Sphere radius
  !!   r_sq   -> Square of radius r (r^2)
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(quadSurface) :: sphere
    private
    real(defReal) :: radius = ZERO, radiusSquared = ZERO
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: evaluate
    procedure :: distance
    procedure :: entersPositiveHalfspace
    procedure :: kill
    procedure :: getRadius
    procedure :: getRadiusSquared
  end type sphere

contains

  !!
  !! Initialise sphere from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!   fatalError if radius or id < 0.
  !!
  subroutine init(self, dict)
    class(sphere), intent(inout)             :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id
    real(defReal), dimension(:), allocatable :: origin
    real(defReal)                            :: radius
    real(defReal), dimension(6)              :: boundingBox
    character(100), parameter                :: Here = 'init (sphere_class.f90)'

    ! Load id.
    call dict % get(id, 'id')
    call self % setId(id)
    
    ! Load origin.
    call dict % get(origin, 'origin')
    call self % setOrigin(origin, 3)

    ! Load radius and set surface tolerance.
    call dict % get(radius, 'radius')
    if (radius <= ZERO) call fatalError(Here, 'Radius of the sphere must be +ve. Is: '//numToChar(radius)//'.')
    self % radius = radius
    self % radiusSquared = radius * radius
    call self % setSurfTol(TWO * radius * SURF_TOL)

    ! Set bounding box.
    boundingBox(1:3) = origin - radius
    boundingBox(4:6) = origin + radius
    call self % setBoundingBox(boundingBox)
    call self % setType('sphere')

  end subroutine init

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(sphere), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: offsetCoords
    real(defReal)                           :: c

    ! Offset coordinates with respect to sphere's origin and compute c.
    offsetCoords = r - self % getOrigin()
    c = dot_product(offsetCoords, offsetCoords) - self % radiusSquared

  end function evaluate

  !! Function 'distance'
  !!
  !! Basic description:
  !!   Returns the distance to the surface of the sphere.
  !!
  !! Detailed description:
  !!   Consider a sphere of radius R. A point lying on the surface of this sphere satisfies:
  !!
  !!   (x - x0)² + (y - y0)² + (z - z0)² - R² = 0
  !!
  !!   Now consider a point p(d) = [p_x(d), p_y(d), p_z(d)] along a ray of direction u = [u_x, u_y, u_z].
  !!   Its parametric equation can be written as:
  !!
  !!   p(d) = r + d * u
  !!
  !!   Substituting this into the equation for the sphere, we get:
  !!
  !!   (r_x + d * u_x - x0)² + (r_y + d * u_y - y0)² + (r_z + d * u_z - z0)² - R² = 0
  !!=> (d * u_x + (r_x - x0))² + (d * u_y + (r_y - y0))² + (d * u_z + (r_z - z0))² - R² = 0
  !!=> (u_x² + u_y² + u_z²) * d² + 2 * ((r_x - x0) * u_x + (r_y - y0) * u_y + (r_z - z0) * u_z) * d 
  !!   + (r_x - x0)² + (r_y - y0)² + (r_z - z0)² - R² = 0
  !!
  !!   Now, define:
  !!
  !!   a = u_x² + u_y² + u_z²
  !!     = 1 (since norm(u) = 1)
  !!   k = (r_x - x0) * u_x + (r_y - y0) * u_y + (r_z - z0) * u_z
  !!   c = (r_x - x0)² + (r_y - y0)² + (r_z - z0)² - R² (note: this is the value returned by the 
  !!   'evaluate' function)
  !!
  !!   Then the equation becomes:
  !!
  !!   d² + 2kd + c = 0
  !!
  !!   which is quadratic in d. The determinant (delta) is given by:
  !!
  !!   delta = 4k² - 4c
  !!         = 4(k² - c)
  !!
  !!   And the two solutions for d are therefore:
  !!
  !!   d_1 = (-2k + 2sqrt(k² - c)) / 2
  !!       = -k + sqrt(k² - c)
  !!   d_2 = (-2k - 2sqrt(k² - c)) / 2
  !!       = -(k + sqrt(k² - c))
  !!
  !! Arguments:
  !!   r [in] -> Position of the particle.
  !!   u [in] -> Direction of the particle.
  !!
  pure function distance(self, r, u) result(d)
    class(sphere), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal)                           :: d, c, k, delta

    ! Initialise d = INF and calculate k and c.
    d = INF
    k = dot_product(r - self % getOrigin(), u)
    c = self % evaluate(r)
    
    ! Compute delta (technically, delta / 4).
    delta = k * k - c

    ! If delta < ZERO, the solutions are complex. There is no intersection and we can return early.
    if (delta < ZERO) return

    ! Check if particle is within surface tolerance of the sphere.
    if (abs(c) < self % getSurfTol()) then
      ! Update d only if k < ZERO (k >= ZERO corresponds to the particle moving away from the sphere). 
      ! Choose maximum distance and return.
      if (k < ZERO) d = -k + sqrt(delta)
      return

    end if

    ! If reached here, update d depending on the sign of c set d = INF if distance is negative.
    d = -(k + sign(sqrt(delta), c))
    if (d <= ZERO) d = INF

  end function distance

  !!
  !! Returns TRUE if particle is going into +ve halfspace
  !!
  !! See surface_inter for details
  !!
  pure function entersPositiveHalfspace(self, r, u) result(halfspace)
    class(sphere), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r, u
    logical(defBool)                        :: halfspace

    halfspace = dot_product(r - self % getOrigin(), u) >= ZERO

  end function entersPositiveHalfspace

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(sphere), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % radius = ZERO
    self % radiusSquared = ZERO

  end subroutine kill

  !!
  !! Returns radius of sphere
  !!
  elemental function getRadius(self) result(radius)
    class(sphere), intent(in) :: self
    real(defReal)             :: radius

    radius = self % radius

  end function getRadius

  !!
  !! Returns square of the radius of the sphere.
  !!
  elemental function getRadiusSquared(self) result(radiusSquared)
    class(sphere), intent(in) :: self
    real(defReal)             :: radiusSquared

    radiusSquared = self % radiusSquared

  end function getRadiusSquared

end module sphere_class