module testFilter_class

  use dictionary_class,           only : dictionary
  use errors_mod,                 only : fatalError
  use genericProcedures,          only : numToChar
  use numPrecision
  use tallyFilter_inter,          only : tallyFilter
  use transportObjectState_class, only : transportObjectState

  implicit none
  private

  !!
  !! Very simple filter used for testing of other components only
  !!   Returns true if state % matIdx (mi) :  minIdx <= mi <= maxIdx
  !!
  type, public, extends(tallyFilter) :: testFilter
    private
    integer(shortInt) :: maxIdx = 0, minIdx = 0
  contains
    procedure :: init
    procedure :: isPass
  end type testFilter

contains

  !!
  !! Initialise testFilter from dictionary
  !!
  subroutine init(self, dict)
    class(testFilter), intent(inout) :: self
    class(dictionary), intent(in)    :: dict
    character(*), parameter          :: here = 'init (testFilter_class.f90)'

    call dict % get(self % minIdx, 'minIdx')
    call dict % get(self % maxIdx, 'maxIdx')

    ! Verify
    if (self % minIdx > self % maxIdx) call fatalError(Here, 'maxIdx < minIdx.')

  end subroutine init

  !!
  !! Returns true if energy value is between specified bounds
  !!
  function isPass(self, state) result(passed)
    class(testFilter), intent(in)           :: self
    class(transportObjectState), intent(in) :: state
    integer(shortInt)                       :: materialIdx
    logical(defBool)                        :: passed

    materialIdx = state % getMaterialIdx()
    passed = self % minIdx <= materialIdx .and. materialIdx <= self % maxIdx

  end function isPass

end module testFilter_class
