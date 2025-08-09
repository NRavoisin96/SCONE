module constantNormalisationMethod_class

  use genericProcedures,         only : fatalError
  use normalisationMethod_inter, only : normalisationMethod
  use numPrecision
  use scoreMemory_class,         only : scoreMemory
  use tallyClerk_inter,          only : tallyClerk

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(normalisationMethod) :: constantNormalisationMethod
    private
  contains
    procedure :: computeNormalisationFactor
  end type constantNormalisationMethod

contains
  !!
  !!
  !!
  function computeNormalisationFactor(self, memory, clerk) result(normalisationFactor)
    class(constantNormalisationMethod), intent(in) :: self
    type(scoreMemory), intent(in)                  :: memory
    class(tallyClerk), intent(in)                  :: clerk
    real(defReal)                                  :: normalisationFactor, normalisationScore
    character(*), parameter                        :: here = 'computeNormalisationFactor (constantNormalisationMethod_class.f90)'

    normalisationScore = memory % getScore(clerk % getMemAddress())
    if (normalisationScore == ZERO) &
    call fatalError(Here, 'Normalisation score from clerk:' //self % getNormalisationClerkName()// 'is 0.')
    normalisationFactor = self % getTargetValue() / normalisationScore

  end function computeNormalisationFactor

end module constantNormalisationMethod_class