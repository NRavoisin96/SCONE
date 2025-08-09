module normalisationMethod_inter

  use dictionary_class,   only : dictionary
  use numPrecision
  use scoreMemory_class,  only : scoreMemory
  use tallyClerk_inter,   only : tallyClerk

  implicit none
  private

  ! Public procedures.
  public :: init, kill

  !!
  !!
  !!
  type, public, abstract :: normalisationMethod
    private
    character(nameLen) :: normalisationClerkName = ''
    real(defReal)      :: targetValue = ZERO
  contains
    procedure(computeNormalisationFactor), deferred :: computeNormalisationFactor
    procedure                                       :: getNormalisationClerkName
    procedure                                       :: getTargetValue
    Procedure                                       :: init
    procedure                                       :: kill
  end type normalisationMethod

  abstract interface
    !!
    !!
    !!
    function computeNormalisationFactor(self, memory, clerk) result(normalisationFactor)
      import                                 :: defReal, normalisationMethod, scoreMemory, tallyClerk
      class(normalisationMethod), intent(in) :: self
      type(scoreMemory), intent(in)          :: memory
      class(tallyClerk), intent(in)          :: clerk
      real(defReal)                          :: normalisationFactor
    end function computeNormalisationFactor

  end interface

contains
  !!
  !!
  !!
  elemental function getNormalisationClerkName(self) result(normalisationClerkName)
    class(normalisationMethod), intent(in) :: self
    character(nameLen)                     :: normalisationClerkName

    normalisationClerkName = self % normalisationClerkName

  end function getNormalisationClerkName

  !!
  !!
  !!
  elemental function getTargetValue(self) result(targetValue)
    class(normalisationMethod), intent(in) :: self
    real(defReal)                          :: targetValue

    targetValue = self % targetValue

  end function getTargetValue

  !!
  !!
  !!
  subroutine init(self, dict)
    class(normalisationMethod), intent(inout) :: self
    class(dictionary), intent(in)             :: dict

    ! Load normalisationClerkName and targetValue.
    call dict % get(self % normalisationClerkName, 'normalisationClerk')
    call dict % get(self % targetValue, 'targetValue')

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(normalisationMethod), intent(inout) :: self

    ! Local.
    self % normalisationClerkName = ''
    self % targetValue = ZERO

  end subroutine kill

end module normalisationMethod_inter