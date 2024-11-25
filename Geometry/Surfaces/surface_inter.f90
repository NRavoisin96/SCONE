module surface_inter

  use numPrecision
  use universalVariables, only : SURF_TOL, VACUUM_BC
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary

  implicit none
  private

  !!
  !! Extandable superclass procedures
  !!
  public :: kill

  !!
  !! Abstract interface for all surfaces
  !!
  !! All surfaces may be represented by some equation F(r) = 0
  !!
  !! By default all surfaces support only VACUUM boundary conditions given
  !! by first entry in BC string. OverrIde relevant procedures in subclasses to change
  !! this behaviour!
  !!
  !! Magnitude of surface tolerance is a property of the surface. By default it is
  !! equal to SURF_TOL parameter.
  !!
  !! Private Members:
  !!   surfId                  -> Id of the surface.
  !!   surfTol                 -> Surface tolerance of the surface.
  !!   boundingBox             -> Axis-aligned bounding box of the surface.
  !!   origin                  -> Origin or offset of the surface. It is an allocatable array as
  !!                              it can have a variable number of elements (3-D components) for
  !!                              different kinds of surfaces.
  !!   type                    -> Type of the surface.
  !!
  !! Interface:
  !!   getId                   -> Return surface Id
  !!   getOrigin               -> Returns the origin or offset of the surface.
  !!   getSurfTol              -> Get value of surface tolerance
  !!   getType                 -> Returns a string with surface type name
  !!   init                    -> Initialise surface from a dictionary
  !!   getBoundingBox          -> Returns definition of axis-aligned bounding box over the surface
  !!   kill                    -> Returns to unitinitialised state
  !!   halfspace               -> Return halfspace ocupied by a particle
  !!   evaluate                -> Return remainder of the surface equation c = F(r)
  !!   distance                -> Return distance to the surface
  !!   entersPositiveHalfspace -> Returns .true. if particle is crossing a boundary into positive
  !!                              halfspace.
  !!   explicitBC              -> Apply explicit BCs
  !!   transformBC             -> Apply transform BCs
  !!
  type, public, abstract :: surface
    private
    integer(shortInt)                            :: surfId = -1
    real(defReal)                                :: surfTol = SURF_TOL
    real(defReal), dimension(6)                  :: boundingBox = ZERO
    real(defReal), dimension(:), allocatable     :: origin
    character(:), allocatable                    :: type
  contains
    ! Initialisation procedures
    procedure                                    :: getBoundingBox
    procedure                                    :: getId
    procedure                                    :: getOrigin
    procedure, non_overridable                   :: getSurfTol
    procedure                                    :: getType
    procedure(init), deferred                    :: init
    procedure                                    :: kill
    procedure                                    :: setBCs
    procedure                                    :: setBoundingBox
    procedure                                    :: setId
    procedure                                    :: setOrigin
    procedure                                    :: setSurfTol
    procedure                                    :: setType

    ! Runtime procedures
    procedure                                    :: halfspace
    procedure(evaluate), deferred                :: evaluate
    procedure(distance), deferred                :: distance
    procedure(entersPositiveHalfspace), deferred :: entersPositiveHalfspace
    procedure                                    :: explicitBC
    procedure                                    :: isWithinSurfTol
    procedure                                    :: transformBC
  end type surface

  abstract interface
    !!
    !! Initialise surface from a dictionary
    !!
    !! Args:
    !!   dict [in] -> Dictionary with surface definition
    !!
    subroutine init(self, dict)
      import :: surface, dictionary
      class(surface), intent(inout) :: self
      class(dictionary), intent(in) :: dict
    end subroutine init

    !!
    !! Evaluate surface expression c = F(r)
    !!
    !! Args:
    !!  r [in] -> Position of the point
    !!
    !! Result:
    !!   Remainder of the surface expression c
    !!
    pure function evaluate(self, r) result(c)
      import :: surface, defReal
      class(surface), intent(in)              :: self
      real(defReal), dimension(3), intent(in) :: r
      real(defReal)                           :: c
    end function evaluate

    !!
    !! Return distance to the surface
    !!
    !! Respects surface transparency.
    !! If position r is such that |F(r)| = |c| < SURF_TOL,
    !! ignore closest crossing in terms of absolute distance |d|
    !! (so for d1 = -eps and d2 = 2*eps d1 would be ignored)
    !!
    !! Args:
    !!  r [in] -> Position of the particle
    !!  u [in] -> Direction of the particle. Assume norm2(u) = 1.0
    !!
    !! Result:
    !!   +ve distance to the next crossing. If there is no crossing
    !!   in +ve direction, returns INF
    !!
    pure function distance(self, r, u) result(d)
      import :: surface, defReal
      class(surface), intent(in)              :: self
      real(defReal), dimension(3), intent(in) :: r
      real(defReal), dimension(3), intent(in) :: u
      real(defReal)                           :: d
    end function distance

    !!
    !! Returns TRUE if particle is going into +ve halfspace
    !!
    !! Is used to determine to which halfspace particle is going when on a surface.
    !!
    !! Args:
    !!   r [in] -> Position of the particle. As the surface (F(r) ~= 0)
    !!   u [in] -> Direction of the partcle. Assume norm2(u) = 1.0
    !!
    !! Result:
    !!   If particle is moving into +ve halfspace (F(r+eps*u0 > 0)) return true
    !!   If particle is moving into -ve halfspace retunr false
    !!
    pure function entersPositiveHalfspace(self, r, u) result(hs)
      import :: surface, defReal, defBool
      class(surface), intent(in)              :: self
      real(defReal), dimension(3), intent(in) :: r, u
      logical(defBool)                        :: hs
    end function entersPositiveHalfspace

  end interface

contains
  !!
  !! Return axis-aligned bounding box
  !!
  !! Note:
  !!   If bounding box is infinite in any axis its smallest value is -INF and highest INF
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Smallest axis-aligned bounding box that contains the entire surface.
  !!   Real array of size 6
  !!   [x_min, y_min, z_min, x_max, y_max, z_max ]
  !!
  pure function getBoundingBox(self) result(boundingBox)
  class(surface), intent(in)  :: self
  real(defReal), dimension(6) :: boundingBox

    boundingBox = self % boundingBox

  end function getBoundingBox

  !!
  !! Return surface Id
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Id -> Id of the surface.
  !!
  elemental function getId(self) result(Id)
    class(surface), intent(in) :: self
    integer(shortInt)          :: Id

    Id = self % surfId

  end function getId

  pure function getOrigin(self) result(origin)
    class(surface), intent(in)                    :: self
    real(defReal), dimension(size(self % origin)) :: origin

    origin = self % origin

  end function getOrigin

  !!
  !! Return value of the surface tolerance for the surface
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Value of the surface tolerance for the surface
  !!
  elemental function getSurfTol(self) result(surfTol)
    class(surface), intent(in) :: self
    real(defReal)              :: surfTol

    surfTol = self % surfTol

  end function getSurfTol

  !!
  !! Return surface type name
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Allocatable string with surface type name without leading or trailing blanks
  !!
  pure function getType(self) result(type)
    class(surface), intent(in) :: self
    character(:), allocatable  :: type

    type = self % type
  end function getType

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(surface), intent(inout) :: self

    self % surfId = -1
    self % surfTol = SURF_TOL
    self % boundingBox = ZERO
    if (allocated(self % origin)) deallocate(self % origin)
    if (allocated(self % type)) deallocate(self % type)

  end subroutine kill

  !!
  !!
  !!
  pure function isWithinSurfTol(self, r) result(isIt)
    class(surface), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    logical(defBool)                        :: isIt

    isIt = abs(self % evaluate(r)) < self % surfTol

  end function

  !!
  !! Return true if particle is in +ve halfspace
  !!
  !! For c = F(r) halfspace is:
  !!   +ve  if c > 0
  !!   -ve  if c <= 0
  !!
  !! Use surface tolerance if |c| = |F(r)| < SURF_TOL then halfspace is
  !! determined by direction of the particle (+ve if is is moving into +ve halfspace;
  !! -ve otherwise; determined with 'going' procedure)
  !!
  !! Args:
  !!   r [in] -> Particle location
  !!   u [in] -> Particle direction (assume norm2(u) = 1.0)
  !!
  !! Result:
  !!   True if position is in +ve halfspace. False if it is in -ve
  !!
  pure function halfspace(self, r, u) result(hs)
    class(surface), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r, u
    logical(defBool)                        :: hs

    ! If |c| < surfTol, determine halfspace based on surface being crossed and particle
    ! direction then return.
    if (self % isWithinSurfTol(r)) then
      hs = self % entersPositiveHalfspace(r, u)
      return
      
    end if

    ! Else halfspace is given based on the sign of c.
    hs = self % evaluate(r) > ZERO

  end function halfspace

      !!
    !! Apply explicit BCs
    !!
    !! Vacuum by default. Override in a subclass to change it!
    !!
    !! Args:
    !!  r [inout] -> Position pre and post BC. Assume that (F(r) ~= 0)
    !!  u [inout] -> Direction pre and post BC. Assume that norm2(u) = 1.0
    !!
  subroutine explicitBC(self, r, u)
    class(surface), intent(in)                 :: self
    real(defReal), dimension(3), intent(inout) :: r, u

    ! Do nothing.
  end subroutine explicitBC

  !!
  !! Apply coordinates transform using BCs
  !!
  !! Vacuum by default. Override in a subclass to change it!
  !!
  !! Args:
  !!  r [inout] -> Position pre- and post-BC.
  !!  u [inout] -> Direction pre- and post-BC. Assume norm2(u) = 1.0
  !!
  subroutine transformBC(self, r, u)
    class(surface), intent(in)                 :: self
    real(defReal), dimension(3), intent(inout) :: r, u

    ! Do nothing.
  end subroutine transformBC

  !!
  !! Sets bounding box.
  !!
  !! Arguments:
  !!   boundingBox [in] -> Smallest axis-aligned bounding box containing the entire surface.
  !!                       Real array of size 6 [x_min, y_min, z_min, x_max, y_max, z_max ].
  !!
  pure subroutine setBoundingBox(self, boundingBox)
    class(surface), intent(inout)           :: self
    real(defReal), dimension(6), intent(in) :: boundingBox

    self % boundingBox = boundingBox

  end subroutine setBoundingBox

  !! Sets boundary conditions
  !!
  !! By default all surfaces support a single vacuum BC.
  !!
  !! Allowable BCs:
  !!   BC(1)   -> entire surface (only vacuum)
  !!
  !! Args:
  !!   BC [in] -> Integer array of arbitrary length with BC flags for different faces
  !!              of the surface. Order of BCs is determined by a surface subclass.
  !!
  !! Errors:
  !!   fatalError BC is not supported.
  !!
  subroutine setBCs(self, BCs)
    class(surface), intent(inout)               :: self
    integer(shortInt), dimension(:), intent(in) :: BCs
    character(100), parameter                   :: Here = 'setBCs (surface_inter.f90)'

    ! Return errors if necessary. Nothing else to be done.
    if (size(BCs) == 0) call fatalError(Here, 'BCs must contain at least one entry.')
    if (BCs(1) /= VACUUM_BC) call fatalError(Here, self % getType()//' only supports VACUUM BCs. Was given: '//&
                                             numToChar(BCs(1)))

  end subroutine setBCs

  !! Sets surface Id.
  !!
  !! Args:
  !!   id [in] -> Id of the surface.
  !!
  !! Errors:
  !!   fatalError if id is < 1.
  !!
  subroutine setId(self, id)
    class(surface), intent(inout) :: self
    integer(shortInt), intent(in) :: id
    character(100), parameter     :: here = 'setId (surface_inter.f90)'

    if (id < 1) call fatalError(here, 'Surface Id must be +ve. Is: '//numToChar(id)//'.')
    self % surfId = id
  end subroutine setId

  !! Sets origin or offset of the surface.
  !!
  !! Args:
  !!   origin [in]    -> Array containing the components of the surface's origin.
  !!   nRequired [in] -> Required number of origin components.
  !!
  !! Errors:
  !!   fatalError if number of origin components does not equal the number of required components.
  subroutine setOrigin(self, origin, nRequired)
    class(surface), intent(inout)           :: self
    real(defReal), dimension(:), intent(in) :: origin
    integer(shortInt), intent(in)           :: nRequired
    integer(shortInt)                       :: nOrigins
    character(*), parameter                 :: here = 'setOrigin (surface_inter.f90)'

    nOrigins = size(origin)
    if (nOrigins /= nRequired) call fatalError(here, &
    'Origin must have size '//numToChar(nRequired)//'. Has: '//numToChar(nOrigins)//'.')
    self % origin = origin

  end subroutine setOrigin

  !! Sets surface tolerance.
  !!
  !! Surface tolerance is a property of the surface. By default it is equal to the SURF_TOL
  !! parameter, which approximately represents a width of ambiguous range.
  !!
  !! Args:
  !!   surfTol [in] -> Tolerance of the surface.
  !!
  !! Errors:
  !!   fatalError is tolerance < ZERO.
  !!
  subroutine setSurfTol(self, surfTol)
    class(surface), intent(inout) :: self
    real(defReal), intent(in)     :: surfTol
    character(100), parameter     :: here = 'setSurfTol (surface_inter.f90)'

    if (surfTol <= ZERO) call fatalError(here, 'Surface tolerance must be +ve. Is: '//numToChar(surfTol)//'.')
    self % surfTol = surfTol
  end subroutine setSurfTol

  subroutine setType(self, type)
    class(surface), intent(inout) :: self
    character(*), intent(in)      :: type
    character(100), parameter     :: here = 'setType (surface_inter.f90)'

    if (len(type) == 0) call fatalError(here, 'Surface type must contain at least one character.')
    self % type = type
  end subroutine setType

end module surface_inter