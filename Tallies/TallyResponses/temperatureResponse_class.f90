module temperatureResponse_class

  use dictionary_class,      only : dictionary
  use nuclearDatabase_inter, only : nuclearDatabase
  use numPrecision
  use randomWalker_class,    only : castRandomWalkerPtr, randomWalker
  use tallyResponse_inter,   only : tallyResponse
  use transportObject_inter, only : transportObject

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(tallyResponse) :: temperatureResponse
    private
  contains
    procedure :: get
    procedure :: init
    procedure :: kill
  end type temperatureResponse

contains
  !!
  !!
  !!
  subroutine get(self, object, val, xsData)
    class(temperatureResponse), intent(in)          :: self
    class(transportObject), intent(in)              :: object
    real(defReal), intent(out)                      :: val
    class(nuclearDatabase), intent(inout), optional :: xsData
    type(randomWalker), pointer                     :: randomWalkerPtr

    ! Downcast transportObject to randomWalker.
    randomWalkerPtr => castRandomWalkerPtr(object, .true.)
    val = randomWalkerPtr % getAccumulatedValue()

  end subroutine get

  !!
  !!
  !!
  subroutine init(self, dict)
    class(temperatureResponse), intent(inout) :: self
    class(dictionary), intent(in)             :: dict

    ! Do nothing.

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(temperatureResponse), intent(inout) :: self

    ! Do nothing.

  end subroutine kill

end module temperatureResponse_class