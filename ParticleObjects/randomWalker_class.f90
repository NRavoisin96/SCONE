module randomWalker_class

  use coordList_class, only : coordList
  use numPrecision
  use RNG_class,       only : RNG

  implicit none
  private

  ! Public procedures.
  public :: newRandomWalker

  !!
  !!
  !!
  type, public :: randomWalker
    private
    class(RNG), pointer :: RNGPtr => null()
    logical(defBool)    :: isDead = .false.
    integer(shortInt)   :: geometryIdx = 0
    real(defReal)       :: accumulatedValue = ZERO, time = ZERO, weight = ZERO
    type(coordList)     :: coords
  contains
    procedure :: accumulateValue
    procedure :: getAccumulatedValue
    procedure :: getCoordsPtr
    procedure :: setAccumulatedValue
  end type randomWalker

contains
  !!
  !!
  !!
  elemental subroutine accumulateValue(self, value)
    class(randomWalker), intent(inout) :: self
    real(defReal), intent(in)          :: value

    self % accumulatedValue = self % accumulatedValue + value

  end subroutine accumulateValue

  !!
  !!
  !!
  elemental function getAccumulatedValue(self) result(accumulatedValue)
    class(randomWalker), intent(in) :: self
    real(defReal)                   :: accumulatedValue

    accumulatedValue = self % accumulatedValue

  end function getAccumulatedValue

  !!
  !!
  !!
  function getCoordsPtr(self) result(coordsPtr)
    class(randomWalker), target, intent(in) :: self
    type(coordList), pointer                :: coordsPtr

    coordsPtr => self % coords

  end function getCoordsPtr

  !!
  !!
  !!
  function newRandomWalker() result(new)
    type(randomWalker) :: new

  end function newRandomWalker

  !!
  !!
  !!
  elemental subroutine setAccumulatedValue(self, value)
    class(randomWalker), intent(inout) :: self
    real(defReal), intent(in)          :: value

    self % accumulatedValue = value

  end subroutine setAccumulatedValue

end module randomWalker_class