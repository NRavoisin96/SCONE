module transportOperatorWoS_class

  use tallyAdmin_class,        only : tallyAdmin
  use transportObject_inter,   only : transportObject
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
  subroutine transit(self, object, tally)
    class(transportOperatorWoS), intent(inout) :: self
    class(transportObject), intent(inout)      :: object
    type(tallyAdmin), intent(inout)            :: tally

  end subroutine transit

end module transportOperatorWoS_class