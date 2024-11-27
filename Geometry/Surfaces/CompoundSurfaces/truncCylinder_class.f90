module truncCylinder_class

  use numPrecision
  use universalVariables
  use genericProcedures,     only : fatalError, numToChar, swap, isEqual, computeQuadraticSolutions
  use dictionary_class,      only : dictionary
  use surface_inter,         only : kill_super => kill
  use compoundSurface_inter, only : compoundSurface

  implicit none
  private

  !!
  !! Finite length cylinder aligned with one of the co-ord axis (x,y or z)
  !!
  !! F(r) = max[ {(r1 - o1)^2 + (r2 - o2)^2 - R^2}/(2R); abs(r3 - o3) - a]
  !!
  !! r -> position; o -> origin; R -> radius; a -> halfwidth
  !! 1,2 -> planar axis, 3-> cylinder axis
  !! Denominator 2R is required for constant surface tolerance thickness
  !!
  !! Three different types are available
  !!   xTruncCylinder -> aligned with X-axis
  !!   yTruncCylinder -> aligned with Y-axis
  !!   zTruncCylinder -> aligned with Z-axis
  !!
  !! Surface Tolerance: SURF_TOL
  !!
  !! Sample Dictionary Input:
  !!  x {type xTruncCylinder; id 2; origin (1.0 -2.0 0.0); halfwidth 1.2; radius 0.5;}
  !!  y {type yTruncCylinder; id 2; origin (1.0 2.0 7.0); halfwidth 1.3; radius 1.5;}
  !!
  !! Boundary Conditions:
  !!   BC order: a_min, a_max
  !!   Where a is cylinder axis [x,y,z]
  !!   BC in radial direction of the cylinder is always VACUUM
  !!
  !! Private Members:
  !!   origin    -> Position of the centre of the cylinder
  !!   halfwidth -> Axial halfwidth (>0.0)
  !!   radius    -> Radius (>0.0)
  !!   axis      -> Cylinder axis specifier
  !!   plane     -> Planar axis specifiers
  !!   BCs       -> Boundary conditions flags [a_min, a_max]
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(compoundSurface) :: truncCylinder
    private
    real(defReal)                   :: radius = ZERO
    integer(shortInt)               :: axis = 0, nBCs = 2
    integer(shortInt), dimension(2) :: planes = 0

  contains
    ! Superclass procedures
    procedure :: init
    procedure :: evaluate
    procedure :: distance
    procedure :: entersPositiveHalfspace
    procedure :: kill
    procedure :: setBCs
    procedure :: explicitBC
    procedure :: transformBC

  end type truncCylinder

contains

  !!
  !! Initialise truncCylinder from a dictionary
  !!
  !! See surface_inter for more details
  !!
  subroutine init(self, dict)
    class(truncCylinder), intent(inout)      :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id, axis
    integer(shortInt), dimension(2)          :: planes
    real(defReal), dimension(:), allocatable :: origin
    real(defReal)                            :: radius, halfwidth
    character(nameLen)                       :: type
    real(defReal), dimension(6)              :: boundingBox
    character(100), parameter                :: Here = 'init (truncCylinder_class.f90)'

    ! Load id.
    call dict % get(id, 'id')
    call self % setId(id)

    ! Load origin.
    call dict % get(origin, 'origin')
    call self % setOrigin(origin, 3)

    ! Load radius.
    call dict % get(radius, 'radius')
    if (radius <= ZERO) call fatalError(Here, 'Radius must be +ve. Is: '//numToChar(radius)//'.')
    self % radius = radius

    ! Load halfwidth.
    call dict % get(halfwidth, 'halfwidth')
    call self % setHalfwidths([halfwidth], 1)

    ! Load type. Set axis and planes accordingly.
    call dict % get(type, 'type')
    select case (type)
      case ('xTruncCylinder')
        axis = X_AXIS
        planes = [Y_AXIS, Z_AXIS]

      case ('yTruncCylinder')
        axis = Y_AXIS
        planes = [X_AXIS, Z_AXIS]

      case ('zTruncCylinder')
        axis = Z_AXIS
        planes = [X_AXIS, Y_AXIS]

      case default
        call fatalError(Here, 'Unknown type of truncCylinder: '//type//'.')

    end select
    self % axis = axis
    self % planes = planes

    ! Set bounding box.
    boundingBox(axis) = origin(axis) - halfwidth
    boundingBox(axis + 3) = origin(axis) + halfwidth
    boundingBox(planes) = origin(planes) - radius
    boundingBox(planes + 3) = origin(planes) + radius
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
    class(truncCylinder), intent(in)        :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: offsetCoords
    real(defReal)                           :: offsetCoordAxis, radius, c
    integer(shortInt), dimension(2)         :: planes
    integer(shortInt)                       :: axis
    real(defReal), dimension(2)             :: offsetCoordsPlanes

    ! Offset coordinates with respect to cylinder's origin, retrieve indices of the
    ! planes and axis of the cylinder then retrieve offset coordinates components
    ! along the planes and axis of the cylinder and cylinder radius.
    offsetCoords = r - self % getOrigin()
    planes = self % planes
    axis = self % axis
    offsetCoordsPlanes = offsetCoords(planes)
    offsetCoordAxis = offsetCoords(axis)
    radius = self % radius

    ! Evaluate surface expressions and return overall maximum.
    c = HALF * (dot_product(offsetCoordsPlanes, offsetCoordsPlanes) / radius - radius)
    c = max(c, self % evaluateCompound([offsetCoordAxis]))

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! See surface_inter for details
  !!
  !! Uses intersection method analogous to box_class
  !! For min and max distance calculation for cylinder see cylinder_class
  !!
  pure function distance(self, r, u) result(d)
    class(truncCylinder), intent(in)         :: self
    real(defReal), dimension(3), intent(in)  :: r, u
    real(defReal)                            :: d, offsetCoordAxis, uAxis, radius, a, k, c, delta, &
                                                dMin, dMax, dBase, surfTol
    real(defReal), dimension(3)              :: offsetCoords
    integer(shortInt), dimension(2)          :: planes
    integer(shortInt)                        :: axis, nSolutions
    logical(defBool)                         :: cannotIntersect
    real(defReal), dimension(2)              :: offsetCoordsPlanes, uPlanes
    real(defReal), dimension(:), allocatable :: solutions

    ! Initialise d = INF.
    d = INF
    
    offsetCoords = r - self % getOrigin()
    axis = self % axis
    planes = self % planes
    offsetCoordAxis = offsetCoords(axis)
    offsetCoordsPlanes = offsetCoords(planes)
    uAxis = u(axis)
    uPlanes = u(planes)
    radius = self % radius
    a = ONE - uAxis * uAxis
    k = dot_product(offsetCoordsPlanes , uPlanes)
    c = dot_product(offsetCoordsPlanes, offsetCoordsPlanes) - radius * radius
    
    ! Compute delta (technically, delta / 4).
    delta = k * k - a * c
    if (delta < ZERO) return

    ! Compute distances to intersection with the cylinder surface. If there are no solutions (ray is
    ! parallel to the cylindrical surface) and the ray is outside the cylinder return early.
    solutions = computeQuadraticSolutions(a, k, c, delta)
    nSolutions = size(solutions)
    if (nSolutions == 0 .and. c > ZERO) return

    ! Initialise dMin = ZERO and dMax = INF. If nSolutions > 0 update dMin and dMax.
    dMin = ZERO
    dMax = INF
    if (nSolutions > 0) then
      dMin = max(minval(solutions), dMin)
      dMax = min(maxval(solutions), dMax)

    end if

    call self % distancesCompound([offsetCoordAxis], [uAxis], dMin, dMax, cannotIntersect)
    c = max(HALF * c / radius, self % evaluateCompound([offsetCoordAxis]))
    call self % chooseDistance(dMin, dMax, abs(c) < self % getSurfTol(), cannotIntersect, d)

  end function distance

  !!
  !! Returns TRUE if particle is going into +ve halfspace
  !!
  !! See surface_inter for details
  !!
  !! Selects closest face (by normal distance to the face -> value of c in surf-expresion).
  !! Then it uses the projection of direction on the normal of that face to select next
  !! halfspace.
  !!
  pure function entersPositiveHalfspace(self, r, u) result(isHalfspacePositive)
    class(truncCylinder), intent(in)        :: self
    real(defReal), dimension(3), intent(in) :: r, u
    logical(defBool)                        :: isHalfspacePositive
    real(defReal), dimension(3)             :: offsetCoords
    integer(shortInt), dimension(2)         :: planes
    integer(shortInt)                       :: axis
    real(defReal), dimension(2)             :: offsetCoordsPlanes, uPlanes
    real(defReal)                           :: offsetCoordAxis, uAxis, radius, surfTol, &
                                               cCylinder, proj
    
    offsetCoords = r - self % getOrigin()
    axis = self % axis
    offsetCoordAxis = offsetCoords(axis)
    uAxis = u(axis)

    isHalfspacePositive = self % isHalfspacePositive([offsetCoordAxis], [uAxis])
    if (isHalfspacePositive) return

    planes = self % planes
    offsetCoordsPlanes = offsetCoords(planes)
    uPlanes = u(planes)
    radius = self % radius

    cCylinder = HALF * (dot_product(offsetCoordsPlanes, offsetCoordsPlanes) / radius - radius)
    surfTol = self % getSurfTol()
    if (abs(cCylinder) < surfTol) then
      proj = dot_product(offsetCoordsPlanes, uPlanes)
      if (isEqual(proj, ZERO)) then
        isHalfspacePositive = cCylinder >= ZERO

      else
        isHalfspacePositive = proj > ZERO

      end if

    end if

  end function entersPositiveHalfspace

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(truncCylinder), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Compound surface.
    call self % killCompound()

    ! Local.
    self % radius = ZERO
    self % axis = 0
    self % planes = 0

  end subroutine kill

  !!
  !! Set boundary conditions
  !!
  !! See surface_inter for details
  !!
  subroutine setBCs(self, BCs)
    class(truncCylinder), intent(inout)         :: self
    integer(shortInt), dimension(:), intent(in) :: BCs

    call self % setCompoundBCs(self % nBCs, BCs)

  end subroutine setBCs

  !!
  !! Apply explicit BCs
  !!
  !! See surface_inter for details
  !!
  !! Approach analogous to box_class applied to axial direction only
  !!
  pure subroutine explicitBC(self, r, u)
    class(truncCylinder), intent(in)           :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    integer(shortInt)                          :: axis
    real(defReal), dimension(3)                :: origin

    ! Retrieve cylinder axis and origin then call parent procedure.
    axis = self % axis
    origin = self % getOrigin()
    call self % explicitCompoundBCs([axis], [origin(axis)], r, u)

  end subroutine explicitBC

  !!
  !! Apply co-ordinate transform BC
  !!
  !! See surface_inter for details
  !!
  !! Approach analogous to box_class applied to axial direction only
  !!
  pure subroutine transformBC(self, r, u)
    class(truncCylinder), intent(in)           :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    integer(shortInt)                          :: axis
    real(defReal), dimension(3)                :: origin

    ! Retrieve cylinder axis and origin then call parent procedure.
    axis = self % axis
    origin = self % getOrigin()
    call self % transformCompoundBCs([axis], [origin(axis)], r, u)

  end subroutine transformBC

end module truncCylinder_class