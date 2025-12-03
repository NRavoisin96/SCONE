module massFraction_class

  use atomicDensitiesCalculator_inter, only : atomicDensitiesCalculator
  use nuclideInfo_class,               only : nuclideInfo
  use numPrecision
  use universalVariables,              only : CENTIMETRES_SQUARED_PER_BARN, CUBIC_METRES_PER_CUBIC_CENTIMETRE, &
                                              KILOGRAMS_PER_ATOMIC_MASS_UNIT

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(atomicDensitiesCalculator) :: massFraction
    private
  contains
    procedure :: computeAtomicDensities
  end type massFraction

contains
  !!
  !!
  !!
  pure subroutine computeAtomicDensities(self, density, array, nuclides)
    class(massFraction), intent(in)                          :: self
    real(defReal), intent(in)                                :: density
    real(defReal), dimension(:), intent(in)                  :: array
    type(nuclideInfo), dimension(size(array)), intent(inout) :: nuclides
    integer(shortInt)                                        :: i
    real(defReal), dimension(size(array))                    :: atomicDensities

    call nuclides % setDensity(array * density * CUBIC_METRES_PER_CUBIC_CENTIMETRE * CENTIMETRES_SQUARED_PER_BARN / &
                               (nuclides % getMassNumber() * KILOGRAMS_PER_ATOMIC_MASS_UNIT * sum(array)))

  end subroutine computeAtomicDensities

end module massFraction_class