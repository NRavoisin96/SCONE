module box_class

  use numPrecision
  use universalVariables
  use genericProcedures,     only : append, fatalError, numToChar, isEqual
  use dictionary_class,      only : dictionary
  use surface_inter,         only : kill_super => kill
  use compoundSurface_inter, only : compoundSurface

  implicit none
  private

  !!
  !! Axis Aligned Box
  !!
  !! F(r) = maxval(abs(r - o) - a)
  !!
  !! Where: a -> halfwidth vector, o-> origin position
  !!        maxval -> maximum element (L_inf norm)
  !!
  !! Surface Tolerance: SURF_TOL
  !!
  !! Sample Dictionary Input:
  !!   aab { type box; id 92; origin (0.0 0.0 9.0); halfwidth(1.0 2.0 0.3);}
  !!
  !! Boundary Conditions:
  !!   BC order: x_min, x_max, y_min, y_max, z_min, z_max
  !!   Each face can have diffrent BC. Any combination is supported with coordinates transform.
  !!
  !! Private Members:
  !!   origin     -> 3-D coordinates of the middle of the box
  !!   halfwidths -> Half-lengths of the box in each direction (must be > ZERO)
  !!   BCs        -> Boundary conditions flags for each face (x_min, x_max, y_min, y_max, z_min, z_max)
  !!
  !! Interface:
  !!   compoundSurface interface
  !!
  type, public, extends(compoundSurface) :: box
    private
    integer(shortInt), dimension(6) :: BCs = VACUUM_BC

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
  end type box

contains
  !!
  !! Initialise box from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!
  subroutine init(self, dict)
    class(box), intent(inout)                :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id
    real(defReal), dimension(:), allocatable :: origin, halfwidths
    real(defReal), dimension(6)              :: boundingBox

    ! Load id.
    call dict % get(id, 'id')
    call self % setId(id)

    ! Load origin.
    call dict % get(origin, 'origin')
    call self % setOrigin(origin, 3)

    ! Load halfwidths.
    call dict % get(halfwidths, 'halfwidth')
    call self % setHalfwidths(halfwidths, 3)

    ! Set bounding box.
    boundingBox(1:3) = origin - halfwidths
    boundingBox(4:6) = origin + halfwidths
    call self % setBoundingBox(boundingBox)
    call self % setType('box')

  end subroutine init
  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(box), intent(in)                  :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c

    ! Compute c from origin-centred coordinates and box halfwidths.
    c = self % evaluateHalfwidths(r - self % getOrigin())

  end function evaluate

  
  pure function distance(self, r, u) result(d)
    class(box), intent(in)                  :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal)                           :: d, dMin, dMax
    logical(defBool)                        :: cannotIntersect
    
    ! Initialise d = INF, dMin = ZERO and dMax = INF.
    d = INF
    dMin = ZERO
    dMax = INF
    call self % closestDist(r - self % getOrigin(), u, dMin, dMax, cannotIntersect)
    call self % chooseDistance(dMin, dMax, self % evaluate(r) < self % getSurfTol(), cannotIntersect, d)
  
  end function distance

  !! Function 'entersPositiveHalfspace'
  !!
  !! Basic description:
  !!   Returns .true. if the particle is going into positive halfspace.
  !!
  !! Detailed description:
  !!   Works similarly to the 'distance' function, except that here we are only interested
  !!   in the value of tMax (renamed t). The algorithm loops over all dimensions and checks
  !!   whether the particle is within surface tolerance of a halfwidth for the current 
  !!   dimension. If the direction component of a dimension for which the particle is within 
  !!   surface tolerance to a halfwidth is ZERO (ie, the ray is parallel to the halfwidth), 
  !!   then halfspace is determined based on which side of the halfwidth the particle lies. 
  !!   Else, t is updated for the current dimension. 
  !!
  !!   Then, the particle will cross into a positive halfspace if t < maxDist, where maxDist 
  !!   is the maximum possible distance that the particle can travel before either crossing 
  !!   one (or more) halfwidths or being outside surface tolerance of all three halfwidths. 
  !!   This distance is obtained by considering a particle lying within surface tolerance to 
  !!   all three halfwidths (ie, very close to a corner of the box) and having a direction 
  !!   vector u = [sqrt(3) / 3, sqrt(3) / 3, sqrt(3) / 3]; in this case applying Pythagoras 
  !!   theorem in 3-D gives maxDist = sqrt(3) * (surface tolerance). 
  !!
  !! Arguments:
  !!   r [in] -> Coordinates of the ray's origin.
  !!   u [in] -> Direction of the ray.
  !!
  pure function entersPositiveHalfspace(self, r, u) result(isHalfspacePositive)
    class(box), intent(in)                  :: self
    real(defReal), dimension(3), intent(in) :: r, u
    logical(defBool)                        :: isHalfspacePositive

    isHalfspacePositive = self % isHalfspacePositive(r - self % getOrigin(), u)

  end function entersPositiveHalfspace

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(box), intent(inout) :: self

    ! Superclass.
     call kill_super(self)

    ! Compound.
    call self % killHalfwidths()

    ! Local.
    self % BCs = VACUUM_BC

  end subroutine kill

  !!
  !! Set boundary conditions
  !!
  !! See surface_inter for details
  !!
  subroutine setBCs(self, BCs)
    class(box), intent(inout)                   :: self
    integer(shortInt), dimension(:), intent(in) :: BCs
    integer(shortInt)                           :: i
    character(100), parameter                   :: Here = 'setBCs (box_class.f90)'

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
  !! Apply explicit boundary conditions.
  !!
  !! See surface_inter for details
  !!
  !! Note:
  !!   - Go through all directions in order to account for corners
  !!
  subroutine explicitBC(self, r, u)
    class(box), intent(in)                     :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    integer(shortInt)                          :: i, bc
    real(defReal)                              :: surfTol, offsetCoord, halfwidth
    real(defReal), dimension(3)                :: halfwidths, origin
    character(100), parameter                  :: Here = 'explicitBC (box_class.f90)'

    ! Retrieve surface tolerance then loop over all dimensions.
    surfTol = self % getSurfTol()
    halfwidths = self % getHalfwidths()
    origin = self % getOrigin()
    do i = 1, 3
      ! Retrieve halfwidth and offset coordinate with respect to box centre
      ! for the current dimension.
      halfwidth = halfwidths(i)
      offsetCoord = r(i) - origin(i)

      ! Skip if particle is well inside the domain.
      if (abs(offsetCoord)  <= halfwidth - surfTol) cycle

      ! Retrieve appropriate BC based on whether offsetCoordinate < ZERO and apply it.
      bc = self % BCs(2 * i + min(0, sign(1, floor(offsetCoord))))
      select case(bc)
        case (VACUUM_BC)
          ! Do nothing. Pass

        case (REFLECTIVE_BC)
          u(i) = -u(i)

        case (PERIODIC_BC)
          ! Calculate displacement and perform translation
          r(i) = r(i) - TWO * sign(halfwidth, offsetCoord)

        case default
          call fatalError(Here, 'Unrecognised BC: '//numToChar(bc)//'.')
  
      end select

    end do

  end subroutine explicitBC

  !!
  !! Apply co-ordinate transform BC
  !!
  !! See surface_inter for details
  !!
  !! Note:
  !!   - Order of transformations does not matter
  !!   - Calculate distance (in # of transformations) for each direction and apply them
  !!
  subroutine transformBC(self, r, u)
    class(box), intent(in)                     :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    real(defReal), dimension(3)                :: halfwidths, origin
    integer(shortInt), dimension(3)            :: nTransforms
    integer(shortInt)                          :: i, j, bc
    real(defReal)                              :: offsetCoord, originComponent, halfwidth
    character(100), parameter                  :: Here = 'transformBC (box_class.f90)'

    ! Calculate number of transforms in each direction. Reduce halfwidths by surface
    ! tolerance to catch particles at the surface.
    halfwidths = self % getHalfwidths()
    origin = self % getOrigin()
    nTransforms = ceiling(abs(r - origin) / (halfwidths - self % getSurfTol())) / 2

    ! Loop over all dimensions.
    do i = 1, 3
      ! Retrieve origin and halfwidth and apply number of transforms for the current dimension.
      originComponent = origin(i)
      halfwidth = halfwidths(i)
      
      do j = 1, nTransforms(i)
        ! Offset coordinates with respect to box centre.
        offsetCoord = r(i) - originComponent

        ! Retrieve appropriate BC based on whether offsetCoordinate < ZERO and apply it.
        bc = self % BCs(2 * i + min(0, sign(1, floor(offsetCoord))))
        select case(bc)
          case (VACUUM_BC)
            ! Do nothing. Pass

          case (REFLECTIVE_BC)
            ! Perform reflection based on distance to the closest plane to offsetCoord.
            r(i) = -r(i) + TWO * (sign(halfwidth, offsetCoord) + originComponent)
            u(i) = -u(i)

          case (PERIODIC_BC)
            ! Perform translation.
            r(i) = r(i) - TWO * sign(halfwidth, offsetCoord)

          case default
            call fatalError(Here, 'Unrecognised BC: '// numToChar(bc)//'.')
        end select

      end do

    end do

  end subroutine transformBC

end module box_class