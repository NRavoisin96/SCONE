module CEPhoton_class

  use CEParticle_inter,           only : CEParticle
  use CEPhotonState_class,        only : CEPhotonState
  use numPrecision
  use transportObjectState_class, only : transportObjectState
  use universalVariables,         only : lightSpeed, P_PHOTON_CE

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(CEParticle) :: CEPhoton
    private
  contains
    procedure :: allocateState
    procedure :: getSpeed
    procedure :: getType
  end type CEPhoton

contains
  !!
  !!
  !!
  subroutine allocateState(self, state)
    class(CEPhoton), intent(inout)                        :: self
    class(transportObjectState), allocatable, intent(out) :: state

    allocate(CEPhotonState :: state)

  end subroutine allocateState

  !!
  !!
  !!
  function getSpeed(self) result(speed)
    class(CEPhoton), intent(in) :: self
    real(defReal)               :: speed

    speed = lightSpeed

  end function getSpeed

  !!
  !!
  !!
  elemental function getType(self) result(type)
    class(CEPhoton), intent(in) :: self
    integer(shortInt)           :: type

    type = P_PHOTON_CE

  end function getType

end module CEPhoton_class