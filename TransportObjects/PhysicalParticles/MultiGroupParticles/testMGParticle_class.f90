module testMGParticle_class

  use MGParticle_inter,           only : MGParticle
  use MGParticleState_class,      only : MGParticleState
  use numPrecision
  use transportObjectState_class, only : transportObjectState
  use universalVariables,         only : P_TEST_TRANSPORT_OBJECT

  implicit none
  private
  
  !!
  !!
  !!
  type, public, extends(MGParticle) :: testMGParticle
    private
  contains
    procedure :: allocateState
    procedure :: getType
  end type testMGParticle

contains
  !!
  !!
  !!
  subroutine allocateState(self, state)
    class(testMGParticle), intent(inout)                  :: self
    class(transportObjectState), allocatable, intent(out) :: state

    allocate(MGParticleState :: state)

  end subroutine allocateState

  !!
  !!
  !!
  elemental function getType(self) result(type)
    class(testMGParticle), intent(in) :: self
    integer(shortInt)                 :: type

    type = P_TEST_TRANSPORT_OBJECT

  end function getType

end module testMGParticle_class