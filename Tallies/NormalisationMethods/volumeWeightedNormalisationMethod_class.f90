module volumeWeightedNormalisationMethod_class

  use normalisationMethod_inter, only : normalisationMethod
  use numPrecision
  use scoreMemory_class,         only : scoreMemory
  use tallyClerk_inter,          only : tallyClerk

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(normalisationMethod) :: volumeWeightedNormalisationMethod
    private
  contains
    procedure :: computeNormalisationFactor
  end type volumeWeightedNormalisationMethod

contains
  !!
  !!
  !!
  function computeNormalisationFactor(self, memory, clerk) result(normalisationFactor)
    class(volumeWeightedNormalisationMethod), intent(in) :: self
    type(scoreMemory), intent(in)                        :: memory
    class(tallyClerk), intent(in)                        :: clerk
    real(defReal)                                        :: normalisationFactor, volumeWeightedSum

    ! Request the clerk to compute the volume-weighted sum of its scores.
    volumeWeightedSum = clerk % computeVolumeWeightedSum(memory)
    if (volumeWeightedSum == ZERO) then
      normalisationFactor = ONE

    else
      normalisationFactor = self % getTargetValue() / volumeWeightedSum

    end if

  end function computeNormalisationFactor

end module volumeWeightedNormalisationMethod_class