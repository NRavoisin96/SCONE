module aPlane_class

  use numPrecision
  use universalVariables, only : X_AXIS, Y_AXIS, Z_AXIS, INF
  use genericProcedures,  only : fatalError, isEqual
  use dictionary_class,   only : dictionary
  use quadSurface_inter,  only : quadSurface
  use surface_inter,      only : kill_super => kill
  
  implicit none
  private

  !!
  !! (Axis) Plane
  !!
  !! Planes with a normal which is one of the principal axis (x,y or z)
  !!
  !! F(r) = (r1 -a0)
  !!
  !! Surface tolerance: SURF_TOL
  !!
  !! Sample dictionary input:
  !!  x { type xPlane; id 16; x0 3.0;}
  !!  y { type yPlane; id 19, y0 -1.2;}
  !!
  !! Private members:
  !!   axis -> Axis indentifier
  !!   a0   -> Offset. Position of the plane on the axis
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(quadSurface) :: aPlane
    private
    integer(shortInt) :: axis = 0
    real(defReal)     :: a0 = ZERO
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: evaluate
    procedure :: distance
    procedure :: entersPositiveHalfspace
    procedure :: kill
  end type aPlane

contains

  !!
  !! Initialise aPlane from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!   -fatalError if id < 1.
  !!   -fatalError if plane type is not recognised.
  !!
  subroutine init(self, dict)
    class(aPlane), intent(inout)  :: self
    class(dictionary), intent(in) :: dict
    integer(shortInt)             :: id, axis
    character(nameLen)            :: type
    real(defReal), dimension(6)   :: boundingBox
    character(100), parameter     :: Here = 'init (aPlane_class.f90)'

    ! Load id.
    call dict % get(id, 'id')
    call self % setId(id)

    ! Load type.
    call dict % get(type, 'type')
    select case(type)
      case('xPlane')
        self % axis = X_AXIS
        call dict % get(self % a0, 'x0')

      case('yPlane')
        self % axis = Y_AXIS
        call dict % get(self % a0, 'y0')

      case('zPlane')
        self % axis = Z_AXIS
        call dict % get(self % a0, 'z0')

      case default
        call fatalError(Here, 'Unknown type of axis plane: '//type)
    end select

    ! Set bounding box.
    axis = self % axis
    boundingBox(1:3) = -INF
    boundingBox(4:6) = INF
    boundingBox([axis, axis + 3]) = self % a0
    call self % setBoundingBox(boundingBox)
    call self % setType(type)

  end subroutine init

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(aPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c

    c = r(self % axis) - self % a0

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! See surface_inter for details
  !!
  pure function distance(self, r, u) result(d)
    class(aPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r, u
    integer(shortInt)                       :: axis
    real(defReal)                           :: d, offsetCoord, uAxis

    ! Initialise d = INF and retrieve offset coordinates and direction components 
    ! along the plane axis.
    d = INF
    axis = self % axis
    offsetCoord = self % a0 - r(axis)
    uAxis = u(axis)

    ! If particle is within surface tolerance of the plane or parallel to the plane
    ! axis, there is no intersection and we can return early.
    if (abs(offsetCoord) < self % getSurfTol() .or. isEqual(uAxis, ZERO)) return
    
    ! Update d and cap distance at infinity.
    d = offsetCoord / uAxis
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
  pure function entersPositiveHalfspace(self, r, u) result(isHalfspacePositive)
    class(aPlane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r, u
    logical(defBool)                        :: isHalfspacePositive
    integer(shortInt)                       :: axis
    real(defReal)                           :: uAxis

    ! Retrieve plane axis and direction component along the axis.
    axis = self % axis
    uAxis = u(axis)

    ! If particle direction is parallel to the plane, halfspace is determined using
    ! position along this plane's axis.
    if (isEqual(uAxis, ZERO)) then
      isHalfspacePositive = r(axis) - self % a0 >= ZERO
      return

    end if

    isHalfspacePositive = uAxis > ZERO

  end function entersPositiveHalfspace

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(aPlane), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % axis = 0
    self % a0 = ZERO

  end subroutine kill

end module aPlane_class