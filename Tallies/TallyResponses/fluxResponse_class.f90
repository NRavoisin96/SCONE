module fluxResponse_class

  use dictionary_class,      only : dictionary
  use nuclearDatabase_inter, only : nuclearDatabase
  use numPrecision
  use tallyResponse_inter,   only : tallyResponse
  use transportObject_inter, only : transportObject

  implicit none
  private

  !!
  !! tallyResponse to score flux contribution
  !!
  !! Always returns ONE
  !!
  !! Interface:
  !!   tallyResponse Interface
  !!
  type, public, extends(tallyResponse) :: fluxResponse
    private
  contains
    procedure :: init
    procedure :: get
    procedure :: kill
  end type fluxResponse

contains

  !!
  !! Initialise Response from dictionary
  !!
  !! See tallyResponse_inter for details
  !!
  subroutine init(self, dict)
    class(fluxResponse), intent(inout) :: self
    class(dictionary), intent(in)      :: dict

    ! Do nothing

  end subroutine init

  !!
  !! Get 1.0 (Response to score flux)
  !!
  !! See tallyResponse_inter for details
  !!
  subroutine get(self, object, val, xsData)
    class(fluxResponse), intent(in)                 :: self
    class(transportObject), intent(in)              :: object
    real(defReal), intent(out)                      :: val
    class(nuclearDatabase), intent(inout), optional :: xsData

    val = ONE

  end subroutine get

  !!
  !! Return to uninitialised State
  !!
  elemental subroutine kill(self)
    class(fluxResponse), intent(inout) :: self

    ! Do nothing for nothing can be done

  end subroutine kill

end module fluxResponse_class
