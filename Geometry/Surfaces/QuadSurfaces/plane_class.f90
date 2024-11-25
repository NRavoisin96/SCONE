module plane_class

  use numPrecision
  use universalVariables, only : X_AXIS, Y_AXIS, Z_AXIS, INF
  use genericProcedures,  only : fatalError, numToChar, isEqual
  use dictionary_class,   only : dictionary
  use quadSurface_inter,  only : quadSurface
  use surface_inter,      only : kill_super => kill
  
  implicit none
  private

  !!
  !! Generic plane surface
  !!
  !! F(r) = c1 * x + c2 * y + c3 * z - c4
  !!
  !! Surface tolerance: SURF_TOL
  !!
  !! Sample dictionary input:
  !!  pl { type plane; id 16; coeffs (1.0 2.0 1.0 0.0);}
  !!
  !! Private members:
  !!   norm   -> Normalised normal vector [c1, c2, c3]
  !!   offset -> Normalised offset [c4]
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(quadSurface) :: plane
    private
    real(defReal), dimension(3) :: norm = ZERO
    real(defReal)               :: offset = ZERO
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: evaluate
    procedure :: distance
    procedure :: entersPositiveHalfspace
    procedure :: kill
  end type plane

contains

  !!
  !! Initialise plane from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!   -fatalError if coefficients do not have 4 entries.
  !!   -fatalError if coefficients 1-3 are all ZERO (invalid normal).
  !!
  subroutine init(self, dict)
    class(plane), intent(inout)              :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id, N
    real(defReal), dimension(:), allocatable :: coeffs
    real(defReal)                            :: offset
    real(defReal), dimension(6)              :: boundingBox
    character(100), parameter                :: Here = 'init (plane_class.f90)'

    ! Load id.
    call dict % get(id, 'id')
    call self % setId(id)
    
    ! Load coefficients.
    call dict % get(coeffs, 'coeffs')
    N = size(coeffs)
    if (N /= 4) call fatalError(Here, 'Coefficients must have size 4. Has: '//numToChar(N)//'.')
    if (isEqual(coeffs(1:3), ZERO)) call fatalError(Here, 'Invalid plane normal. Coefficients 1-3 are all 0.')

    ! Set normal vector and offset.
    coeffs = coeffs / norm2(coeffs(1:3))
    offset = coeffs(4)
    self % norm = coeffs(1:3)
    self % offset = offset

    ! Initialise bounding box and check if plane is aligned with one axis. Check for
    ! perfect equality here.
    boundingBox(1:3) = -INF
    boundingBox(4:6) = INF
    if (coeffs(2) == ZERO .and. coeffs(3) == ZERO) boundingBox([1, 4]) = offset
    if (coeffs(1) == ZERO .and. coeffs(3) == ZERO) boundingBox([2, 5]) = offset
    if (coeffs(1) == ZERO .and. coeffs(2) == ZERO) boundingBox([3, 6]) = offset
    call self % setBoundingBox(boundingBox)
    call self % setType('plane')

  end subroutine init

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(plane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c

    c = dot_product(r, self % norm) - self % offset

  end function evaluate

  !! Function 'distance'
  !!
  !! Basic description:
  !!   Returns the distance to the plane.
  !!
  !! Detailed description:
  !!   Consider a generic plane. A point lying on this plane satisfies:
  !!
  !!   c1 * x + c2 * y + c3 * z - c4 = 0
  !!
  !!   Now consider a point p(d) = [p_x(d), p_y(d), p_z(d)] along a ray of direction u = [u_x, u_y, u_z].
  !!   Its parametric equation can be written as:
  !!
  !!   p(d) = r + d * u
  !!
  !!   Substituting this into the equation for the plane, we get:
  !!
  !!   c1 * (r_x + d * u_x) + c2 * (r_y + d * u_y) + c3 * (r_z + d * u_z) - c4 = 0
  !!=> (c1 * u_x + c2 * u_y + c3 * u_z) * d + c1 * r_x + c2 * u_y + c3 * u_z - c4 = 0
  !!
  !!   Now, define:
  !!
  !!   k = c1 * u_x + c2 * u_y + c3 * u_z
  !!   c = c1 * r_x + c2 * u_y + c3 * u_z - c4 (note: this is the value returned by the 'evaluate' function)
  !!
  !!   Then the equation becomes:
  !!
  !!   k * d + c = 0
  !!
  !!   which is linear in d. The solution is simply given by:
  !!
  !!   d = - c / k
  !!
  !! Arguments:
  !!   r [in] -> Position of the particle.
  !!   u [in] -> Direction of the particle.
  !!
  pure function distance(self, r, u) result(d)
    class(plane), intent(in)                :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal)                           :: d, k, c

    ! Initialise d = INF then compute k and c.
    d = INF
    k = dot_product(u, self % norm)
    c = self % evaluate(r)

    ! If the particle's direction is parallel to the plane (k = ZERO) or if the particle
    ! is within surface tolerance of the plane, there is no intersection and we can return
    ! early.
    if (isEqual(k, ZERO) .or. abs(c) < self % getSurfTol()) return
    
    ! If reached here, update d and cap resulting distance at infinity.
    d = -c / k
    if (d <= ZERO .or. d > INF) d = INF

  end function distance

  !!
  !! Returns TRUE if particle is going into +ve halfspace
  !!
  !! See surface_inter for details
  !!
  !! Note:
  !!   For parallel direction halfspace is asigned by the sign of `evaluate` result.
  !!
  pure function entersPositiveHalfspace(self, r, u) result(halfspace)
    class(plane), intent(in)                :: self
    real(defReal), dimension(3), intent(in) :: r, u
    logical(defBool)                        :: halfspace
    real(defReal)                           :: proj

    proj = dot_product(u, self % norm)

    ! Special case of parallel direction. Particle stays in its current halfspace.
    if (isEqual(proj, ZERO)) then
      halfspace = self % evaluate(r) >= ZERO
      return
      
    end if

    halfspace = proj > ZERO

  end function entersPositiveHalfspace

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(plane), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % norm = ZERO
    self % offset = ZERO

  end subroutine kill

end module plane_class