module squareCylinder_class

  use numPrecision
  use universalVariables
  use genericProcedures,     only : fatalError, numToChar, areEqual
  use dictionary_class,      only : dictionary
  use surface_inter,         only : kill_super => kill
  use compoundSurface_inter, only : compoundSurface

  implicit none
  private

  !!
  !! Square cylinder aligned with x,y or z axis
  !!
  !! F(r) = maxval(abs(r - o) - a)
  !!
  !! Where: a -> 2D halfwidth vector, o-> 2D origin position
  !!        maxval -> maximum element (L_inf norm)
  !!
  !! Three different types are avaliable
  !!   xSquareCylinder -> aligned with X-axis
  !!   ySquareCylinder -> aligned with Y-axis
  !!   zSquareCylinder -> aligned with Z-axis
  !!
  !! Surface Tolerance: SURF_TOL
  !!
  !! Sample Dictionary Input:
  !!   x { type xSquareCylinder; id 92; origin (0.0 0.0 9.0); halfwidth(0.0 2.0 0.3);}
  !!   y { type ySquareCylinder; id 92; origin (0.0 0.0 9.0); halfwidth(2.0 0.0 0.3);}
  !!
  !!   Halfwidth and origin entry in axis along the cylinder direction is ignored.
  !!
  !! Boundary Conditions:
  !!   BCs order: x_min, x_max, y_min, y_max, z_min, z_max
  !!
  !!   Each face can have different BCs. Any combination is supported with co-ordinate transform.
  !!   BCs on all faces (even the infinate ones) must be provided. Of course BCs for planes normal
  !!   to the cylinder axis does not matter.
  !!
  !! Private Members:
  !!   origin    -> position of the middle of the squareCylinder
  !!   halfwidth -> Halfwidths (half-length) of the squareCylinder in each direction (must be > 0.0)
  !!   BCs       -> Boundary conditions flags for each face (x_min, x_max, y_min, y_max, z_min, z_max)
  !!   planes    -> Indices of the planes of the cylinder (e.g. for xSquareCylinder, y-z)
  !!   axis      -> Index of the axis of the cylinder.
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(compoundSurface) :: squareCylinder
    private
    integer(shortInt), dimension(2) :: planes = 0
    integer(shortInt)               :: axis = 0, nBCs = 4
  contains
    ! Superclass procedures
    procedure                       :: init
    procedure                       :: evaluate
    procedure                       :: distance
    procedure                       :: entersPositiveHalfspace
    procedure                       :: kill
    procedure                       :: setBCs
    procedure                       :: explicitBC
    procedure                       :: transformBC
  end type squareCylinder

contains
  !!
  !! Initialise squareCylinder from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!  - fatalError if id < 1.
  !!  - fatalError if cylinder type is unrecognised.
  !!  - fatalError if origin does not contain three entries.
  !!  - fatalError if halfwidth does not contain three entries.
  !!
  subroutine init(self, dict)
    class(squareCylinder), intent(inout)     :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt), dimension(2)          :: planes
    integer(shortInt)                        :: id, axis
    character(nameLen)                       :: type
    real(defReal), dimension(:), allocatable :: temp
    real(defReal), dimension(2)              :: origin, halfwidths
    real(defReal), dimension(6)              :: boundingBox
    character(100), parameter                :: Here = 'init (squareCylinder_class.f90)'

    ! Load id.
    call dict % get(id, 'id')
    call self % setId(id)

    ! Load type. Set axis and planes accordingly.
    call dict % get(type, 'type')
    select case(type)
      case('xSquareCylinder')
        axis = X_AXIS
        planes = [Y_AXIS, Z_AXIS]

      case('ySquareCylinder')
        axis = Y_AXIS
        planes = [X_AXIS, Z_AXIS]

      case('zSquareCylinder')
        axis = Z_AXIS
        planes = [X_AXIS, Y_AXIS]

      case default
        call fatalError(Here, 'Unknown type of square cylinder: '//type//'.')

    end select
    self % axis = axis
    self % planes = planes

    ! Load origin.
    call dict % get(temp,'origin')
    origin = temp(planes)
    call self % setOrigin(origin, 2)

    ! Load halfwidths.
    call dict % get(temp,'halfwidth')
    halfwidths = temp(planes)
    call self % setHalfwidths(halfwidths, origin, 2)

    ! Set bounding box.
    boundingBox(planes) = origin - halfwidths
    boundingBox(planes + 3) = origin + halfwidths
    boundingBox(axis) = -INF
    boundingBox(axis + 3) = INF
    call self % setBoundingBox(boundingBox)
    call self % setType(type)

    ! Initialise BCs.
    call self % setCompoundBCs(self % nBCs)

  end subroutine init

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(squareCylinder), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c

    ! Compute c from origin-centred coordinates and cylinder halfwidths along the
    ! bounding planes.
    c = self % evaluateCompound(r(self % planes))

  end function evaluate

  !! Function 'distance'
  !!
  !! Basic description:
  !!   Computes the distance to the nearest intersection between a ray of origin r and direction
  !!   u and the bounding planes of the square cylinder. Returns INF if the ray does not intersect 
  !!   the bounding planes.
  !!
  !! Detailed description:
  !!   Intersection test is based on the slab method. See box_class.f90 for details.
  !!
  !! Arguments:
  !!   r [in] -> Coordinates of the ray's origin.
  !!   u [in] -> Direction of the ray.
  !!
  !! Notes:
  !!   For a square cylinder, we only need to check for intersection in two dimensions instead of 
  !!   three for a box.
  !!
  pure function distance(self, r, u) result(d)
    class(squareCylinder), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal)                           :: d
    integer(shortInt), dimension(2)         :: planes

    ! Call compoundSurface procedure.
    planes = self % planes
    d = self % distancesCompound(r(planes), u(planes), -INF, INF, abs(self % evaluate(r)) < self % getSurfTol())

  end function distance

  !! Function 'entersPositiveHalfspace'
  !!
  !! Basic description:
  !!   Returns .true. if the particle is going into positive halfspace.
  !!
  !! Detailed description:
  !!   See box_class.f90.
  !!
  !! Arguments:
  !!   r [in] -> Coordinates of the ray's origin.
  !!   u [in] -> Direction of the ray.
  !!
  !! Notes:
  !!   For a square cylinder, we only need to check for intersection in two dimensions instead of 
  !!   three for a box.
  !!
  pure function entersPositiveHalfspace(self, r, u) result(isHalfspacePositive)
    class(squareCylinder), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r, u
    logical(defBool)                        :: isHalfspacePositive
    integer(shortInt), dimension(2)         :: planes

    planes = self % planes
    isHalfspacePositive = self % isHalfspacePositive(r(planes), u(planes))

  end function entersPositiveHalfspace

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(squareCylinder), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Compound.
    call self % killCompound()

    ! Local.
    self % planes = 0
    self % axis = 0

  end subroutine kill

  !!
  !! Set boundary conditions
  !!
  !! See surface_inter for details
  !!
  subroutine setBCs(self, BCs)
    class(squareCylinder), intent(inout)        :: self
    integer(shortInt), dimension(:), intent(in) :: BCs
    integer(shortInt), dimension(2)             :: planes
    integer(shortInt), dimension(self % nBCs)   :: idxArray

    ! Create BCs array.
    planes = self % planes
    idxArray([1, 3]) = 2 * planes - 1
    idxArray([2, 4]) = 2 * planes
    call self % setCompoundBCs(self % nBCs, BCs(idxArray))

  end subroutine setBCs

  !!
  !! Apply explicit BCs
  !!
  !! See surface_inter for details
  !!
  !! Note:
  !!   - Go through all directions in order to account for corners
  !!
  pure subroutine explicitBC(self, r, u)
    class(squareCylinder), intent(in)          :: self
    real(defReal), dimension(3), intent(inout) :: r, u

    call self % explicitCompoundBCs(self % planes, r, u)

  end subroutine explicitBC

  !!
  !! Apply co-ordinate transform BCs
  !!
  !! See surface_inter for details
  !!
  !! Note:
  !!   - Order of transformations does not matter
  !!   - Calculate distance (in # of transformations) for each direction and apply them
  !!
  pure subroutine transformBC(self, r, u)
    class(squareCylinder), intent(in)          :: self
    real(defReal), dimension(3), intent(inout) :: r, u

    call self % transformCompoundBCs(self % planes, r, u)

  end subroutine transformBC

end module squareCylinder_class