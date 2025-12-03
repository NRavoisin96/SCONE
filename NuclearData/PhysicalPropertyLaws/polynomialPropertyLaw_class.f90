module polynomialPropertyLaw_class

  use dictionary_class,          only : dictionary
  use errors_mod,                only : fatalError
  use numPrecision
  use physicalPropertyLaw_inter, only : kill_super => kill, physicalPropertyLaw

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(physicalPropertyLaw) :: polynomialPropertyLaw
    private
    real(defReal), dimension(:), allocatable :: coefficients
  contains
    procedure :: computeProperty
    procedure :: init
    procedure :: kill
  end type polynomialPropertyLaw

contains
  !!
  !!
  !!
  function computeProperty(self, stateVariable) result(property)
    class(polynomialPropertyLaw), intent(in) :: self
    real(defReal), intent(in)                :: stateVariable
    integer(shortInt)                        :: i
    real(defReal)                            :: property
    character(*), parameter                  :: HERE = 'computeProperty (polynomialPropertyLaw_class.f90)'

    if (.not. allocated(self % coefficients)) call fatalError(HERE, 'Unallocated coefficients.')
    
    property = ZERO
    do i = 1, size(self % coefficients)
      property = property + self % coefficients(i) * stateVariable ** (i - 1)

    end do

  end function computeProperty

  !!
  !!
  !!
  subroutine init(self, dict)
    class(polynomialPropertyLaw), intent(inout) :: self
    type(dictionary), intent(in)                :: dict

    call dict % get(self % coefficients, 'coefficients')

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(polynomialPropertyLaw), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    if (allocated(self % coefficients)) deallocate(self % coefficients)

  end subroutine kill

end module polynomialPropertyLaw_class