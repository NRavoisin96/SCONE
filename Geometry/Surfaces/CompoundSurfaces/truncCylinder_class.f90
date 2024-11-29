module truncCylinder_class

  use numPrecision
  use universalVariables
  use genericProcedures,      only : fatalError, numToChar, areEqual, solveQuadratic
  use dictionary_class,       only : dictionary
  use surface_inter,          only : kill_super => kill
  use compoundSurface_inter,  only : compoundSurface

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
    real(defReal)                        :: radius = ZERO
    integer(shortInt)                    :: axis = 0
    integer(shortInt), dimension(2)      :: planes = 0

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

    ! Load halfwidth.
    call dict % get(halfwidth, 'halfwidth')
    call self % setHalfwidths([halfwidth], [origin(axis)], 1)

    ! Set bounding box.
    boundingBox(axis) = origin(axis) - halfwidth
    boundingBox(axis + 3) = origin(axis) + halfwidth
    boundingBox(planes) = origin(planes) - radius
    boundingBox(planes + 3) = origin(planes) + radius
    call self % setBoundingBox(boundingBox)
    call self % setType(type)

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
    real(defReal)                           :: radius, c
    real(defReal), dimension(2)             :: offsetCoordsPlanes

    ! Offset coordinates with respect to cylinder's origin, retrieve indices of the
    ! planes and axis of the cylinder then retrieve offset coordinates components
    ! along the planes and axis of the cylinder and cylinder radius.
    offsetCoords = r - self % getOrigin()
    offsetCoordsPlanes = offsetCoords(self % planes)
    radius = self % radius

    ! Evaluate surface expressions and return overall maximum.
    c = HALF * (dot_product(offsetCoordsPlanes, offsetCoordsPlanes) / radius - radius)
    c = max(c, self % evaluateCompound(r(self % axis)))

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
    real(defReal)                            :: d, uAxis, radius, a, k, c, delta, &
                                                dMin, dMax, surfTol
    real(defReal), dimension(3)              :: offsetCoords
    integer(shortInt), dimension(2)          :: planes
    integer(shortInt)                        :: axis, nSolutions
    real(defReal), dimension(2)              :: offsetCoordsPlanes, uPlanes
    real(defReal), dimension(:), allocatable :: solutions
    logical(defBool)                         :: surfTolCondition

    ! Initialise d = INF.
    d = INF
    
    offsetCoords = r - self % getOrigin()
    axis = self % axis
    planes = self % planes
    offsetCoordsPlanes = offsetCoords(planes)
    uAxis = u(axis)
    uPlanes = u(planes)
    radius = self % radius
    a = ONE - uAxis * uAxis
    k = dot_product(offsetCoordsPlanes , uPlanes)
    c = dot_product(offsetCoordsPlanes, offsetCoordsPlanes) - radius * radius

    ! Compute delta (technically, delta / 4). Return early if delta < ZERO.
    delta = k * k - a * c
    if (delta < ZERO) return

    ! Compute distances to intersection with the cone surface and initialise dMin = ZERO and
    ! dMax = INF.
    solutions = solveQuadratic(a, k, c, delta)
    nSolutions = size(solutions)
    c = HALF * c / radius
    surfTol = self % getSurfTol()
    if (nSolutions == 0 .and. c > -surfTol) return
    dMin = -INF
    dMax = INF

    ! If there are no solutions then the ray is parallel to the cone opening and fully contained in its 
    ! surface. If there is only one solution then the ray is parallel to the cone opening. If this solution 
    ! is negative and the particle is inside the cone (c < ZERO) then the ray can only intersect with one of 
    ! the cone bases so compute this intersection and return.
    if (nSolutions > 0) then
      ! Pre-compute bound and loop over all solutions retrieved.
      dMin = minval(solutions)
      dMax = maxval(solutions)

    end if

    ! Update dMin and dMax, compute distances to intersection with the cone's halfwidth and choose correct
    ! distance.
    surfTolCondition = abs(max(c, self % evaluateCompound([r(axis)]))) < surfTol
    d = self % distancesCompound([r(axis)], [uAxis], dMin, dMax, surfTolCondition)

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
    real(defReal), dimension(2)             :: offsetCoordsPlanes
    real(defReal)                           :: radius, c, proj
    
    axis = self % axis
    isHalfspacePositive = self % isHalfspacePositive([r(axis)], [u(axis)])
    if (isHalfspacePositive) return

    offsetCoords = r - self % getOrigin()
    planes = self % planes
    offsetCoordsPlanes = offsetCoords(planes)
    radius = self % radius

    c = HALF * (dot_product(offsetCoordsPlanes, offsetCoordsPlanes) / radius - radius)
    if (abs(c) < self % getSurfTol()) then
      proj = dot_product(offsetCoordsPlanes, u(planes))
      if (areEqual(proj, ZERO)) then
        isHalfspacePositive = c >= ZERO

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

    ! Halfwidth surface.
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

    call self % setCompoundBCs(2, BCs)

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

    ! Call compoundSurface procedure.
    call self % explicitCompoundBCs([self % axis], r, u)

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

    ! Call compoundSurface procedure.
    call self % transformCompoundBCs([self % axis], r, u)

  end subroutine transformBC

end module truncCylinder_class