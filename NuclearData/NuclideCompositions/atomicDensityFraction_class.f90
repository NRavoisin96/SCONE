module atomicDensityFraction_class

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
  type, public, extends(atomicDensitiesCalculator) :: atomicDensityFraction
    private
  contains
    procedure :: computeAtomicDensities
  end type atomicDensityFraction

contains
  !!
  !!
  !!
  pure subroutine computeAtomicDensities(self, density, array, nuclides)
    class(atomicDensityFraction), intent(in)                 :: self
    real(defReal), intent(in)                                :: density
    real(defReal), dimension(:), intent(in)                  :: array
    type(nuclideInfo), dimension(size(array)), intent(inout) :: nuclides

    call nuclides % setDensity(array * density * CUBIC_METRES_PER_CUBIC_CENTIMETRE * CENTIMETRES_SQUARED_PER_BARN / &
                               (dot_product(array, nuclides % getMassNumber()) * KILOGRAMS_PER_ATOMIC_MASS_UNIT))

  end subroutine computeAtomicDensities

end module atomicDensityFraction_class