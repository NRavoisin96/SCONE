module rawAtomicDensities_class

  use atomicDensitiesCalculator_inter, only : atomicDensitiesCalculator
  use nuclideInfo_class,               only : nuclideInfo
  use numPrecision

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(atomicDensitiesCalculator) :: rawAtomicDensities
    private
  contains
    procedure :: computeAtomicDensities
  end type rawAtomicDensities

contains
  !!
  !!
  !!
  pure subroutine computeAtomicDensities(self, density, array, nuclides)
    class(rawAtomicDensities), intent(in)                    :: self
    real(defReal), intent(in)                                :: density
    real(defReal), dimension(:), intent(in)                  :: array
    type(nuclideInfo), dimension(size(array)), intent(inout) :: nuclides

    call nuclides % setDensity(array)

  end subroutine computeAtomicDensities

end module rawAtomicDensities_class