module rootUniverse_class

  use cell_inter,         only : cell
  use cellShelf_class,    only : cellShelf
  use charMap_class,      only : charMap
  use cylinder_class,     only : cylinder
  use dictionary_class,   only : dictionary
  use genericProcedures,  only : fatalError, numToChar
  use meshShelf_class,    only : meshShelf
  use numPrecision
  use publicObjects,      only : coordData
  use surface_inter,      only : surface
  use surfaceShelf_class, only : surfaceShelf
  use universalVariables, only : OUTSIDE_MAT
  use universe_inter,     only : universe, kill_super => kill, charToFill

  implicit none
  private

  ! Parameters
  integer(shortInt), parameter :: INSIDE_ID = 1, OUTSIDE_ID = 2

  !!
  !! A top level (root) universe of geometry
  !!
  !! Is composed of two regions. Inside and outside separated by a single surface.
  !! Inside is the -ve halfspace of the boundary surface
  !! +ve halfspace is OUTSIDE
  !! Filling can be universe given by ID (`u<6>` syntax) or a material given by name (e.g. 'fuel')
  !!
  !! Local ID 1 is inside. 2 is outside.
  !!
  !! Sample Input Dictionary:
  !!   root { type rootUniverse;
  !!          id 7;
  !!          border 78;   // Boundary surface
  !!          fill u<17>;  // Inside filling
  !!        }
  !!
  !!
  !! Public Members:
  !!   surf    -> Pointer to the boundary surface
  !!   surfIdx -> Index of the boundary surface
  !!
  !! Interface:
  !!   Universe interface
  !!   border -> Return surfIdx of the boundary surface
  !!
  type, public, extends(universe) :: rootUniverse
    class(surface), pointer       :: surf => null()
    integer(shortInt)             :: surfIdx = 0
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: findCell
    procedure :: distance
    procedure :: cross

    ! Subclass procedures
    procedure :: border

  end type rootUniverse


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
    class(rootUniverse), intent(inout)                        :: self
    class(dictionary), intent(in)                             :: dict
    type(charMap), intent(in)                                 :: mats
    integer(shortInt), dimension(:), allocatable, intent(out) :: fills
    type(cellShelf), intent(inout)                            :: cells
    type(surfaceShelf), intent(inout)                         :: surfs
    type(meshShelf), intent(inout)                            :: meshes
    integer(shortInt)                                         :: id
    character(nameLen)                                        :: name
    character(100), parameter                                 :: Here = 'init (rootUniverse_class.f90)'

    ! Setup the base class
    ! With: id, origin rotations...
    call self % setupBase(dict)

    ! Make sure root does not contain origin nor rotation.
    if (dict % isPresent('origin')) call fatalError(Here, 'Origin is not allowed. Centre of the root universe is &
                                                    &always (0.0 0.0 0.0).')
    if (dict % isPresent('rotation')) call fatalError(Here, 'Rotation is not allowed. Root universe cannot be rotated.')

    ! Get boundary surface.
    call dict % get(id, 'border')
    if (id < 1) call fatalError(Here, 'Border id must be +ve. Inside is always in &
                                &-ve halfspace. Was given: '//numToChar(id)//'.')

    self % surfIdx = surfs % getIdx(id)
    self % surf => surfs % getPtr(self % surfIdx)

    ! Create fill array.
    allocate(fills(2))
    fills(OUTSIDE_ID) = OUTSIDE_MAT
    call dict % get(name, 'fill')
    fills(INSIDE_ID) = charToFill(name, mats, Here)

  end subroutine init

  !!
  !! Find local cell ID given a point
  !!
  !! See universe_inter for details.
  !!
  pure subroutine findCell(self, data)
    class(rootUniverse), intent(inout) :: self
    type(coordData), intent(inout)     :: data

    ! Initialise localId = INSIDE_ID then check halfspace.
    data % localId = merge(OUTSIDE_ID, INSIDE_ID, self % surf % halfspace(data % r, data % u))

  end subroutine findCell

  !!
  !! Return distance to the next boundary between local cells in the universe
  !!
  !! See universe_inter for details.
  !!
  pure subroutine distance(self, data)
    class(rootUniverse), intent(inout) :: self
    type(coordData), intent(inout)     :: data

    data % surfaceIdx = self % surfIdx
    data % d = self % surf % distance(data % r, data % u)

  end subroutine distance

  !!
  !! Cross between local cells
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError if surface from distance is not MOVING_IN or MOVING_OUT
  !!
  subroutine cross(self, data)
    class(rootUniverse), intent(inout) :: self
    type(coordData), intent(inout)     :: data

    ! Cross by cell finding in case of significant undershoots
    call self % findCell(data)

  end subroutine cross

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(rootUniverse), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill local
    self % surf => null()
    self % surfIdx = 0

  end subroutine kill

  !!
  !! Return Index of the boundary surface
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Integer IDX of the boundary (border) surface
  !!
  pure function border(self) result(idx)
    class(rootUniverse), intent(in) :: self
    integer(shortInt)               :: idx

    idx = self % surfIdx

  end function border

end module rootUniverse_class