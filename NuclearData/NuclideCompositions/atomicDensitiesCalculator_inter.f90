module atomicDensitiesCalculator_inter

  use dictionary_class,  only : dictionary
  use nuclideInfo_class, only : nuclideInfo
  use numPrecision

  implicit none
  private

  !!
  !!
  !!
  type, public, abstract :: atomicDensitiesCalculator
    private
  contains
    procedure(computeAtomicDensities), deferred :: computeAtomicDensities
  end type atomicDensitiesCalculator

  abstract interface
    !!
    !!
    !!
    pure subroutine computeAtomicDensities(self, density, array, nuclides)
      import                                                   :: atomicDensitiesCalculator, defReal, nuclideInfo
      class(atomicDensitiesCalculator), intent(in)             :: self
      real(defReal), intent(in)                                :: density
      real(defReal), dimension(:), intent(in)                  :: array
      type(nuclideInfo), dimension(size(array)), intent(inout) :: nuclides
    end subroutine computeAtomicDensities

  end interface

end module atomicDensitiesCalculator_inter