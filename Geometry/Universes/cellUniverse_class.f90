module cellUniverse_class

  use numPrecision
  use universalVariables, only : UNDEF_MAT, NUDGE
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use coord_class,        only : coord
  use charMap_class,      only : charMap
  use surfaceShelf_class, only : surfaceShelf
  use cell_inter,         only : cell
  use cellShelf_class,    only : cellShelf
  use meshShelf_class,    only : meshShelf
  use universe_inter,     only : universe, kill_super => kill

  implicit none
  private

  !!
  !! Local helper class to group cell data
  !!
  !! Public members:
  !!   idx           -> cellIdx of the cell in cellShelf.
  !!   ptr           -> Pointer to the cell.
  !!
  type, private :: localCell
    integer(shortInt)                            :: idx = 0
    class(cell), pointer                         :: ptr => null()
  end type localCell

  !!
  !! Representation of a universe via cells
  !!
  !! Each local cell in the universe corresponds to a cell given by an id.
  !! An extra local cell is always defined inside the cellUniverse with UNDEF_MAT
  !! (undefined material) filling. If position is not in any user-defined cell, it is in this
  !! extra cell. Extra cell exists to enable plotting of geometry without fatalErrors.
  !!
  !! Sample Input Dictionary:
  !!   uni { type cellUniverse;
  !!         id 7;
  !!         # origin (2.0 0.0 0.0);    #
  !!         # rotation (23.0 0.0 0.0); #
  !!         cells ( 1 3 4);         }
  !!
  !! Note:
  !!   - Local ids are assigned in order of definition. In the example above local ids would map
  !!     to the following cell ids  [localID: cellID] {1: 1, 2: 3, 3: 4, 4: UNDEF }
  !!   - Cell overlaps are forbidden, but there is no check to find overlaps.
  !! 
  !!   TODO: add check to find overlaps?
  !!
  !! Public Members:
  !!   cells -> Structure that stores cellIdx and pointers to the cells
  !!
  !! Interface:
  !!   universe interface
  !!
  type, public, extends(universe) :: cellUniverse
    type(localCell), dimension(:), allocatable :: cells
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: findCell
    procedure :: distance
    procedure :: cross
    procedure :: cellOffset
  end type cellUniverse

contains

  !!
  !! Initialise Universe
  !!
  !! See universe_inter for details.
  !!
  subroutine init(self, dict, mats, fills, cells, surfs, meshes)
    class(cellUniverse), intent(inout)                        :: self
    class(dictionary), intent(in)                             :: dict
    type(charMap), intent(in)                                 :: mats
    integer(shortInt), dimension(:), allocatable, intent(out) :: fills
    type(cellShelf), intent(inout)                            :: cells
    type(surfaceShelf), intent(inout)                         :: surfs
    type(meshShelf), intent(inout)                            :: meshes
    integer(shortInt), dimension(:), allocatable              :: cellTemp
    integer(shortInt)                                         :: N, i

    ! Setup the base class
    ! With: id, origin rotations...
    call self % setupBase(dict)

    ! Load Cells by ID
    call dict % get(cellTemp, 'cells')
    N = size(cellTemp)

    ! Allocate cells and fill arrays.
    allocate(self % cells(N))
    allocate(fills(N + 1))

    ! Convert cell ids to idxs, load pointers and create fill array.
    self % cells % idx = cellTemp
    do i = 1, N
      self % cells(i) % idx = cells % getIdx(self % cells(i) % idx)
      self % cells(i) % ptr => cells % getPtr(self % cells(i) % idx)
      fills(i) = cells % getFill(self % cells(i) % idx)
      
    end do
    fills(N + 1) = UNDEF_MAT

  end subroutine init

  !!
  !! Find local cell ID given a point
  !!
  !! See universe_inter for details.
  !!
  pure subroutine findCell(self, r, u, localId, cellIdx, elementIdx)
    class(cellUniverse), intent(inout)      :: self
    real(defReal), dimension(3), intent(in) :: r, u
    integer(shortInt), intent(out)          :: localId, cellIdx, elementIdx

    elementIdx = 0
    ! Search all cells
    do localId = 1, size(self % cells)
      if (self % cells(localId) % ptr % inside(r, u)) then
        cellIdx = self % cells(localId) % idx
        return

      end if

    end do

    ! If not found return undefined cell.
    ! Already set to localId (== size(self % cells) + 1) by the do loop
    cellIdx = 0

  end subroutine findCell

  !!
  !! Return distance to the next boundary between local cells in the universe
  !!
  !! See universe_inter for details.
  !!
  !! Errors:
  !!   fatalError if in UNDEFINED cell
  !!
  subroutine distance(self, coords, d, surfIdx)
    class(cellUniverse), intent(inout) :: self
    type(coord), intent(inout)         :: coords
    real(defReal), intent(out)         :: d
    integer(shortInt), intent(out)     :: surfIdx
    integer(shortInt)                  :: localId
    character(100), parameter          :: Here = 'distance (cellUniverse_class.f90)'

    localId = coords % localId

    if (localId > size(self % cells)) call fatalError(Here, &
    'Particle is in undefined cell with local id: '//numToChar(localId)//'.')

    ! Calculate distance
    call self % cells(localId) % ptr % distance(d, surfIdx, coords % r, coords % dir)

  end subroutine distance

  !!
  !! Cross between local cells
  !!
  !! See universe_inter for details.
  !!
  !! Note: Introduces extra movement to the particle to push it over boundary
  !!   for more efficient search. Distance is NUGDE.
  !!
  pure subroutine cross(self, coords, surfIdx)
    class(cellUniverse), intent(inout) :: self
    type(coord), intent(inout)         :: coords
    integer(shortInt), intent(in)      :: surfIdx

    ! NUDGE position slightly forward to escape surface tolerance
    ! and avoid calculating normal and extra dot-products
    coords % r = coords % r + coords % dir * NUDGE

    ! Find cell
    ! TODO: Some cell neighbour list
    call self % findCell(coords % r, coords % dir, coords % localId, coords % cellIdx, coords % elementIdx)

  end subroutine cross

  !!
  !! Return offset for the current cell
  !!
  !! See universe_inter for details.
  !!
  pure function cellOffset(self, coords) result (offset)
    class(cellUniverse), intent(in) :: self
    type(coord), intent(in)         :: coords
    real(defReal), dimension(3)     :: offset

    ! There is no cell offset
    offset = ZERO

  end function cellOffset

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(cellUniverse), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    if(allocated(self % cells)) deallocate(self % cells)

  end subroutine kill

end module cellUniverse_class