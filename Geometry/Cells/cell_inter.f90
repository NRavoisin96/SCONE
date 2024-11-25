module cell_inter

  use numPrecision
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use surfaceShelf_class, only : surfaceShelf

  implicit none
  private

  ! Extendable procedures
  public :: kill

  !!
  !! Abstract interface for all cells
  !!
  !! Cell is intented to represent a volume of space. It is not intended to be
  !! a independent entity, but rather to be used as a component of universes.
  !!
  !! Private Members:
  !!   cellId   -> Id for a particular cell instance
  !!
  !! Interface:
  !!   init     -> Initialise from dictionary and surface Shelf
  !!   inside   -> Return true is a given position is contained inside the cell
  !!   distance -> Assuming the point is inside the cell, calculate distance to the boundary
  !!               and give surfIdx for the surface that will be crossed
  !!  setId     -> Set ID of a cell
  !!  id        -> Return id of a cell
  !!  kill      -> Returns to uninitialised state
  !!
  type, public, abstract :: cell
    private
    integer(shortInt)             :: cellId = 0
  contains
    procedure(init), deferred     :: init
    procedure(inside), deferred   :: inside
    procedure(distance), deferred :: distance
    procedure, non_overridable    :: setId
    procedure, non_overridable    :: id
    procedure                     :: kill
  end type cell

  abstract interface

    !!
    !! Initialise cell
    !!
    !! Args:
    !!   dict [in] -> Dictionary with definition of the cell
    !!   surfs [inout] -> Surface shelf with all user-defined surfaces
    !!
    !! Note: surfs is intent in/out some subclasses of a cell may require adding new surfaces
    !!   to the surface shelf.
    !!
    subroutine init(self, dict, surfs)
      import :: cell, dictionary, surfaceShelf
      class(cell), intent(inout)        :: self
      class(dictionary), intent(in)     :: dict
      type(surfaceShelf), intent(inout) :: surfs
    end subroutine init

    !!
    !! Return TRUE is position is inside the cell
    !!
    !! Args:
    !!   r [in] -> Position
    !!   u [in] -> Normalised direction (norm2(u) = 1.0)
    !!
    !! Result:
    !!   True if position is inside the cell. False otherwise
    !!
    pure function inside(self, r, u, idx) result(isIt)
      import :: cell, defReal, shortInt, defBool
      class(cell), intent(in)                 :: self
      real(defReal), dimension(3), intent(in) :: r, u
      integer(shortInt), intent(in), optional :: idx
      logical(defBool)                        :: isIt
    end function inside

    !!
    !! Return distance to cell boundary
    !!
    !! Assumes that particle is INSIDE the cell.
    !!
    !! Args:
    !!   d [out]       -> Distance to the boundary
    !!   surfIdx [out] -> Index of a surface that will be crossed. If the surface is not defined on
    !!     the surface shelf its value should be -ve. If no surface is hit return 0.
    !!   r [in]        -> Position
    !!   u [in]        -> Normalised direction (norm2(u) = 1.0)
    !!
    !! Note: If cell has some internally defined surfaces it should return -ve values of the
    !!   surfIdx, with the meaning for different -ve numbers being dependant on particular
    !!   subclass of the cell
    !!
    pure subroutine distance(self, d, surfIdx, r, u)
      import :: cell, defReal, shortInt
      class(cell), intent(in)                 :: self
      real(defReal), intent(out)              :: d
      integer(shortInt), intent(out)          :: surfIdx
      real(defReal), dimension(3), intent(in) :: r, u
    end subroutine distance

  end interface

contains

  !!
  !! Sets Id of the cell.
  !!
  !! Args:
  !!   cellId [in] -> Id of the cell.
  !!
  !! Errors:
  !!   fatalError if cellId < 1.
  !!
  subroutine setId(self, cellId)
    class(cell), intent(inout)    :: self
    integer(shortInt), intent(in) :: cellId
    character(100), parameter     :: Here = 'setId (cell_inter.f90)'

    if (cellId < 1) call fatalError(Here, 'Cell ID must be +ve. Is: '//numToChar(cellId))
    self % cellId = cellId

  end subroutine setId

  !!
  !! Return the id of the cell
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Cell id
  !!
  pure function id(self)
    class(cell), intent(in) :: self
    integer(shortInt)       :: id

    id = self % cellId

  end function id

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(cell), intent(inout) :: self

    self % cellId = 0

  end subroutine kill

end module cell_inter