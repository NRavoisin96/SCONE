module testResponse_class

  use dictionary_class,      only : dictionary
  use nuclearDatabase_inter, only : nuclearDatabase
  use numPrecision
  use tallyResponse_inter,   only : tallyResponse
  use transportObject_inter, only : transportObject

  implicit none
  private

  !!
  !! tallyResponse to be used for testing
  !!
  !! Always returns a constant
  !!
  !! Private Members:
  !!   value -> Value that is always returned
  !!
  !! Interface:
  !!   tallyResponse Interface
  !!
  !!
  type, public, extends(tallyResponse) :: testResponse
    private
    real(defReal) :: value = ONE
  contains
    procedure :: init
    procedure :: get
    procedure :: kill
  end type testResponse

contains

  !!
  !! Initialise Response from dictionary
  !!
  !! See tallyResponse_inter for details
  !!
  subroutine init(self, dict)
    class(testResponse), intent(inout) :: self
    class(dictionary), intent(in)      :: dict

    call dict % get(self % value,'value')

  end subroutine init

  !!
  !! Get Response
  !!
  !! Returns "value" Member
  !!
  !! See tallyResponse_inter for details
  !!
  subroutine get(self, object, val, xsData)
    class(testResponse), intent(in)                 :: self
    class(transportObject), intent(in)              :: object
    real(defReal), intent(out)                      :: val
    class(nuclearDatabase), intent(inout), optional :: xsData

    val = self % value

  end subroutine get

  !!
  !! Return to uninitialised State
  !!
  elemental subroutine kill(self)
    class(testResponse), intent(inout) :: self

    self % value = ONE

  end subroutine kill

end module testResponse_class
