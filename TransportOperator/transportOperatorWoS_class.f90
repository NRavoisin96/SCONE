module transportOperatorWoS_class

  use particle_class,          only : particle
  use tallyAdmin_class,        only : tallyAdmin
  use transportOperator_inter, only : transportOperator

  !!
  !!
  !!
  type, public, extends(transportOperator) :: transportOperatorWoS
    private
  contains
    procedure :: transit
  end type transportOperatorWoS

contains
  !!
  !!
  !!
  subroutine transit(self, p, tally)
    class(transportOperatorWoS), intent(inout) :: self
    class(particle), intent(inout)             :: p
    type(tallyAdmin), intent(inout)            :: tally

  end subroutine transit

end module transportOperatorWoS_class