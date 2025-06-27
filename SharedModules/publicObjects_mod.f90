module publicObjects

  use numPrecision
  use RNG_class,          only : RNG
  use universalVariables, only : ZERO

  implicit none
  public

  !!
  !!
  !!
  type :: particleData
    integer(shortInt)   :: matIdx = 0
    real(defReal)       :: E = ZERO
    class(RNG), pointer :: rand => null()
  end type particleData

end module publicObjects