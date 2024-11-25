module squareCylinder_class

  use numPrecision
  use universalVariables
  use genericProcedures,     only : fatalError, numToChar, isEqual
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
    integer(shortInt), dimension(6) :: BCs = VACUUM_BC
    integer(shortInt), dimension(2) :: planes = 0
    integer(shortInt)               :: axis  = 0

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
    call self % setHalfwidths(halfwidths, 2)

    ! Set bounding box.
    boundingBox(planes) = origin - halfwidths
    boundingBox(planes + 3) = origin + halfwidths
    boundingBox(axis) = -INF
    boundingBox(axis + 3) = INF
    call self % setBoundingBox(boundingBox)
    call self % setType(type)

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
    c = self % evaluateHalfwidths(r(self % planes) - self % getOrigin())

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
    real(defReal)                           :: d, dMin, dMax
    integer(shortInt), dimension(2)         :: planes
    logical(defBool)                        :: cannotIntersect

    ! Initialise d = INF, dMin = ZERO and dMax = INF.
    d = INF
    dMin = ZERO
    dMax = INF
    planes = self % planes
    call self % closestDist(r(planes) - self % getOrigin(), u(planes), dMin, dMax, cannotIntersect)
    call self % chooseDistance(dMin, dMax, self % evaluate(r) < self % getSurfTol(), cannotIntersect, d)

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
    isHalfspacePositive = self % isHalfspacePositive(r(planes) - self % getOrigin(), u(planes))

  end function entersPositiveHalfspace

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(squareCylinder), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Compound.
    call self % killHalfwidths()

    ! Local.
    self % BCs = VACUUM_BC
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
    integer(shortInt)                           :: i
    character(100), parameter                   :: Here = 'setBCs (squareCylinder_inter.f90)'

    if(size(BCs) < 6) call fatalError(Here, 'Boundary conditions must have size 6. Has: '//numToChar(size(BCs))//'.')

    ! Verify BC flags.
    do i = 1, 6
      select case(BCs(i))
        case (VACUUM_BC, REFLECTIVE_BC, PERIODIC_BC)
          ! Do nothing, pass
        case default
          call fatalError(Here,'Unrecognised BC: '//numToChar(BCs(i))//' in position: '//numToChar(i)//'.')

      end select
    end do

    ! Verify periodic BCs.
    if(.not. all((BCs([1,3,5]) == PERIODIC_BC) .eqv. (BCs([2,4,6]) == PERIODIC_BC))) call fatalError(Here, &
                                                  'Periodic BCs need to be applied on opposite surfaces.')

    ! Load BC flags.
    self % BCs = BCs(1:6)

  end subroutine setBCs

  !!
  !! Apply explicit BCs
  !!
  !! See surface_inter for details
  !!
  !! Note:
  !!   - Go through all directions in order to account for corners
  !!
  subroutine explicitBC(self, r, u)
    class(squareCylinder), intent(in)          :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    integer(shortInt)                          :: i, plane, bc
    real(defReal)                              :: surfTol, offsetCoord, halfwidth
    real(defReal), dimension(2)                :: halfwidths, origin
    character(100), parameter                  :: Here = 'explicitBC (squareCylinder_class.f90)'

    ! Retrieve surface tolerance and loop over dimensions of the square cylinder's bounding planes.
    surfTol = self % getSurfTol()
    halfwidths = self % getHalfwidths()
    origin = self % getOrigin()
    do i = 1, 2
      ! Retrieve current plane and offset coordinate and halfwidth components along 
      ! the curent plane dimension. Note: because of the mix of 2D and 3D vectors to 
      ! get the correct components we use:
      ! i -> for 2D vectors (origin & halfwidths)
      ! plane -> for 3D vectors (r & u)
      plane = self % planes(i)
      offsetCoord = r(plane) - origin(i)
      halfwidth = halfwidths(i)

      ! Skip if particle is well inside the domain.
      if (abs(offsetCoord)  <= halfwidth - surfTol) cycle

      ! Retrieve appropriate BC based on whether offsetCoordinate < ZERO and apply it.
      bc = self % BCs(2 * plane + min(0, sign(1, floor(offsetCoord))))
      select case(bc)
        case (VACUUM_BC)
          ! Do nothing. Pass

        case (REFLECTIVE_BC)
          u(plane) = -u(plane)

        case (PERIODIC_BC)
          ! Calculate displacement and perform translation
          r(plane) = r(plane) - TWO * sign(halfwidth, offsetCoord)

        case default
          call fatalError(Here, 'Unrecognised BC: '//numToChar(bc)//'.')
  
      end select

    end do

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
  subroutine transformBC(self, r, u)
    class(squareCylinder), intent(in)          :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    real(defReal), dimension(2)                :: halfwidths, origin
    integer(shortInt), dimension(2)            :: planes, nTransforms
    integer(shortInt)                          :: i, j, plane, bc
    real(defReal)                              :: offsetCoord, originComponent, halfwidth
    character(100), parameter                  :: Here = 'transformBC (squareCylinder_class.f90)'

    ! Calculate number of transforms in each direction. Reduce halfwidths by surface
    ! tolerance to catch particles at the surface.
    planes = self % planes
    halfwidths = self % getHalfwidths()
    origin = self % getOrigin()
    nTransforms = ceiling(abs(r(planes) - origin) / (halfwidths - self % getSurfTol())) / 2

    ! Loop over dimensions of the square cylinder's bounding planes.
    do i = 1, 2
      ! Retrieve current plane and origin and halfwidth components along the current
      ! plane dimension. Note: because of the mix of 2D and 3D vectors to get the
      ! correct components we use:
      ! i -> for 2D vectors (origin & halfwidths)
      ! plane -> for 3D vectors (r & u)
      plane = self % planes(i)
      originComponent = origin(i)
      halfwidth = halfwidths(i)
      
      do j = 1, nTransforms(i)
        ! Offset coordinates with respect to box centre.
        offsetCoord = r(plane) - originComponent

        ! Retrieve appropriate BC based on whether offsetCoordinate < ZERO and apply it.
        bc = self % BCs(2 * plane + min(0, sign(1, floor(offsetCoord))))
        select case(bc)
          case (VACUUM_BC)
            ! Do nothing. Pass

          case (REFLECTIVE_BC)
            ! Perform reflection based on distance to the closest plane to offsetCoord.
            r(plane) = -r(plane) + TWO * (sign(halfwidth, offsetCoord) + originComponent)
            u(plane) = -u(plane)

          case (PERIODIC_BC)
            ! Perform translation.
            r(plane) = r(plane) - TWO * sign(halfwidth, offsetCoord)

          case default
            call fatalError(Here, 'Unrecognised BC: '// numToChar(bc)//'.')
        end select

      end do

    end do

  end subroutine transformBC

end module squareCylinder_class