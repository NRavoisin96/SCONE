module constantPropertyLaw_class

  use dictionary_class,          only : dictionary
  use numPrecision
  use physicalPropertyLaw_inter, only : kill_super => kill, physicalPropertyLaw

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(physicalPropertyLaw) :: constantPropertyLaw
    private
    real(defReal) :: value = ZERO
  contains
    procedure :: computeProperty
    procedure :: init
  end type constantPropertyLaw

contains
  !!
  !!
  !!
  pure function computeProperty(self, stateVariable) result(property)
    class(constantPropertyLaw), intent(in) :: self
    real(defReal), intent(in)              :: stateVariable
    real(defReal)                          :: property

    property = self % value

  end function computeProperty

  !!
  !!
  !!
  subroutine init(self, dict)
    class(constantPropertyLaw), intent(inout) :: self
    type(dictionary), intent(in)              :: dict

    call dict % get(self % value, 'value')

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(constantPropertyLaw), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % value = ZERO

  end subroutine kill

end module constantPropertyLaw_class