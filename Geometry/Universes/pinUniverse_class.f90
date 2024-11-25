module pinUniverse_class

  use numPrecision
  use universalVariables, only : INF, targetNotFound
  use genericProcedures,  only : fatalError, numToChar, swap
  use dictionary_class,   only : dictionary
  use coord_class,        only : coord
  use charMap_class,      only : charMap
  use surfaceShelf_class, only : surfaceShelf
  use cylinder_class,     only : cylinder
  use cell_inter,         only : cell
  use cellShelf_class,    only : cellShelf
  use meshShelf_class,    only : meshShelf
  use universe_inter,     only : universe, kill_super => kill, charToFill
  implicit none
  private

  ! Parameters
  ! Are public for use in unit tests
  integer(shortInt), parameter, public :: MOVING_IN = -1, MOVING_OUT = -2

  !!
  !! Universe that represents a single pin.
  !!
  !! Is composed of concentric cylinders. Central cell has localId 1 and the id increases with
  !! subsequent rings.
  !!
  !! Sample Dictionary Input:
  !!   pinUni {
  !!     id 7;
  !!     type pinUniverse;
  !!     #origin (1.0 0.0 0.1);#
  !!     #rotation (30.0 0.0 0.0);#
  !!     radii (3.0 4.5 0.0 1.0 );
  !!     fills (u<3> void clad u<4>);
  !!   }
  !!
  !!  Radii must contain a 0.0 entry, which indicates the outermost annulus (infinite radius).
  !!  `fills` and `radii` are given as pairs by position in the input arrays. Thus, fills
  !!  are sorted together with the radii. As a result, in the example, local cell 1 is filled
  !!  with u<4>, cell 2 with u<3> etc.
  !!
  !! Public Members:
  !!  radiiSquared -> Array of radius squared for each annulus.
  !!  annuli       -> Array of cylinder surfaces corresponding to different annuli.
  !!
  !! Interface:
  !!   universe interface
  !!
  type, public, extends(universe) :: pinUniverse
    private
    real(defReal), dimension(:), allocatable  :: radiiSquared
    type(cylinder), dimension(:), allocatable :: annuli
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: findCell
    procedure :: distance
    procedure :: cross
    procedure :: cellOffset
  end type pinUniverse

contains

  !!
  !! Initialise Universe
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError for invalid input
  !!
  subroutine init(self, dict, mats, fills, cells, surfs, meshes)
    class(pinUniverse), intent(inout)                         :: self
    class(dictionary), intent(in)                             :: dict
    type(charMap), intent(in)                                 :: mats
    integer(shortInt), dimension(:), allocatable, intent(out) :: fills
    type(cellShelf), intent(inout)                            :: cells
    type(surfaceShelf), intent(inout)                         :: surfs
    type(meshShelf), intent(inout)                            :: meshes
    integer(shortInt)                                         :: idx, N, i
    real(defReal), dimension(:), allocatable                  :: radii
    character(nameLen), dimension(:), allocatable             :: fillNames
    character(100), parameter                                 :: Here = 'init (pinUniverse_class.f90)'

    ! Setup the base class
    ! With: id, origin rotations...
    call self % setupBase(dict)

    ! Load radii and fill data.
    call dict % get(radii, 'radii')
    call dict % get(fillNames, 'fills')

    ! Check values.
    if (size(radii) /= size(fillNames)) call fatalError(Here, 'Sizes of radii and fills do not match.')
    if (any(radii < ZERO)) call fatalError(Here, 'Found radius with -ve value.')

    ! Find radius with value 0.0 (representing outermost element) and swap it and corresponding fill to
    ! the last elements of each array. Set radius with 0.0 value to INF.
    N = size(radii)
    idx = minloc(radii, 1)
    if (radii(idx) /= ZERO) call fatalError(Here, 'Could not find outermost element with radius 0.0.')
    call swap(radii(idx), radii(N))
    call swap(fillNames(idx), fillNames(N))
    radii(N) = INF * 1.1_defReal

    ! Sort radii and fills into ascending order.
    do i = N - 1, 1, -1
      idx = maxloc(radii(1:i), 1)
      call swap(radii(idx), radii(i))
      call swap(fillNames(idx), fillNames(i))

    end do

    ! Check for duplicate values of radii.
    do i = 1, N - 1
      if (radii(i) == radii(i + 1)) call fatalError(Here, 'Duplicate value of radius: '//numToChar(radii(i))//'.')

    end do

    ! Load radii squared, build cylinders and create fill array.
    self % radiiSquared = radii * radii
    allocate(self % annuli(N))
    allocate(fills(N))
    do i = 1, N
      call self % annuli(i) % build(id = 1, origin = [ZERO, ZERO, ZERO], type = 'zCylinder', radius = radii(i))
      fills(i) = charToFill(fillNames(i), mats, Here)

    end do

  end subroutine init

  !!
  !! Find local cell ID given a point
  !!
  !! See universe_inter for details.
  !!
  subroutine findCell(self, r, u, localId, cellIdx, tetrahedronIdx)
    class(pinUniverse), intent(inout)       :: self
    real(defReal), dimension(3), intent(in) :: r, u
    integer(shortInt), intent(out)          :: localId, cellIdx, tetrahedronIdx
    real(defReal), dimension(2)             :: rPlanes, uPlanes
    real(defReal)                           :: rPlanesSquared, mul, surfTol

    ! Set cellIdx = 0, retrieve the particle's position and direction components in the
    ! cylinders' planes (all cylinders are zCylinders so these components are 1 and 2) and
    ! compute rPlanesSquared.
    cellIdx = 0
    tetrahedronIdx = 0
    rPlanes = r(1:2)
    uPlanes = u(1:2)
    rPlanesSquared = dot_product(rPlanes, rPlanes)

    ! Find local cell. Pre-compute multiplier based on particle direction.
    mul = sign(ONE, -dot_product(rPlanes, uPlanes))
    do localId = 1, size(self % radiiSquared)
      ! Retrieve surface tolerance of current cylinder and check if particle is inside it.
      surfTol = mul * self % annuli(localId) % getSurfTol()
      if(rPlanesSquared < self % radiiSquared(localId) + surfTol) return

    end do
    ! If reached here localID = size(self % r_sq) + 1

  end subroutine findCell

  !!
  !! Return distance to the next boundary between local cells in the universe
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError is localID is invalid
  !!
  subroutine distance(self, coords, d, surfIdx)
    class(pinUniverse), intent(inout)  :: self
    type(coord), intent(inout)         :: coords
    real(defReal), intent(out)         :: d
    integer(shortInt), intent(out)     :: surfIdx
    real(defReal), dimension(3)        :: r, u
    real(defReal)                      :: d_out, d_in
    integer(shortInt)                  :: localId, nAnnuli
    character(100), parameter          :: Here = 'distance (pinUniverse_class.f90)'

    ! Retrieve localId and number of annuli. Call fatalError if localId is out of bounds.
    localId = coords % localId
    nAnnuli = size(self % annuli)
    if (localId < 1 .or. localId > nAnnuli + 1) call fatalError(Here, 'Invalid local id: '//numToChar(localId)//'.')

    ! Retrieve particle's location and direction components and compute outer and inner distances.
    r = coords % r
    u = coords % dir

    ! Outer distance.
    if (localId > nAnnuli) then
      d_out = INF

    else
      d_out = self % annuli(localId) % distance(r, u)

    end if

    ! Inner distance.
    if (localId == 1) then
      d_in = INF

    else
      d_in = self % annuli(localId - 1) % distance(r, u)

    end if

    ! Select distance and surface.
    d = min(d_in, d_out)
    surfIdx = MOVING_IN
    if (d_out <= d_in) surfIdx = MOVING_OUT

  end subroutine distance

  !!
  !! Cross between local cells
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError if surface from distance is not MOVING_IN or MOVING_OUT
  !!
  subroutine cross(self, coords, surfIdx)
    class(pinUniverse), intent(inout)  :: self
    type(coord), intent(inout)         :: coords
    integer(shortInt), intent(in)      :: surfIdx
    character(100), parameter          :: Here = 'cross (pinUniverse_class.f90)'

    if (surfIdx == MOVING_IN) then
      coords % localID = coords % localID - 1

    else if (surfIdx == MOVING_OUT) then
      coords % localID = coords % localID + 1

    else
      call fatalError(Here, 'Unknown surface memento: '//numToChar(surfIdx)//'.')

    end if

  end subroutine cross

  !!
  !! Return offset for the current cell
  !!
  !! See universe_inter for details.
  !!
  function cellOffset(self, coords) result (offset)
    class(pinUniverse), intent(in) :: self
    type(coord), intent(in)         :: coords
    real(defReal), dimension(3)     :: offset

    ! There is no cell offset
    offset = ZERO

  end function cellOffset

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(pinUniverse), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill local
    if(allocated(self % radiiSquared)) deallocate(self % radiiSquared)
    if(allocated(self % annuli)) deallocate(self % annuli)

  end subroutine kill

end module pinUniverse_class