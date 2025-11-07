module testPhysicalParticle_class

  use numPrecision
  use physicalParticle_inter,      only : physicalParticle
  use physicalParticleState_class, only : physicalParticleState
  use transportObjectState_class,  only : transportObjectState
  use universalVariables,          only : P_TEST_TRANSPORT_OBJECT

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(physicalParticle) :: testPhysicalParticle
    private
  contains
    procedure :: allocateState
    procedure :: getSpeed
    procedure :: getType
  end type testPhysicalParticle

contains
  !!
  !!
  !!
  subroutine allocateState(self, state)
    class(testPhysicalParticle), intent(inout)            :: self
    class(transportObjectState), allocatable, intent(out) :: state

    allocate(physicalParticleState :: state)

  end subroutine allocateState

  !!
  !!
  !!
  function getSpeed(self) result(speed)
    class(testPhysicalParticle), intent(in) :: self
    real(defReal)                           :: speed

    speed = ZERO

  end function getSpeed

  !!
  !!
  !!
  elemental function getType(self) result(type)
    class(testPhysicalParticle), intent(in) :: self
    integer(shortInt)                       :: type

    type = P_TEST_TRANSPORT_OBJECT

  end function getType

end module testPhysicalParticle_class