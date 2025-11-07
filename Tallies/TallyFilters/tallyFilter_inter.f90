module tallyFilter_inter

  use dictionary_class,           only : dictionary
  use numPrecision
  use transportObjectState_class, only : transportObjectState

  implicit none
  private

  !!
  !! Abstract interface of tallyFilters
  !!
  !! All tallyFilters given a particle state return .true. or .false.
  !! This determines whether an event should be tallied or not
  !!
  type, public, abstract :: tallyFilter
    private
  contains
    procedure(init), deferred   :: init
    procedure(isPass), deferred :: isPass
    procedure                   :: isFail
  end type tallyFilter

  abstract interface
    !!
    !! Initialise filter from dictionary
    !!
    subroutine init(self, dict)
      import                            :: dictionary, tallyFilter
      class(tallyFilter), intent(inout) :: self
      class(dictionary), intent(in)     :: dict
    end subroutine init

    !!
    !! Return .true. if state passes filter test
    !! Return .false. otherwise or if test is undefined
    !!
    function isPass(self, state) result(passed)
      import                                  :: defBool, tallyFilter, transportObjectState
      class(tallyFilter), intent(in)          :: self
      class(transportObjectState), intent(in) :: state
      logical(defBool)                        :: passed
    end function isPass

  end interface

contains

  !!
  !! Shorthand for [.not.isPass ] for semantic clarity
  !!
  function isFail(self, state) result(failed)
    class(tallyFilter), intent(in)          :: self
    class(transportObjectState), intent(in) :: state
    logical(defBool)                        :: failed

    failed = .not. self % isPass(state)

  end function isFail

end module tallyFilter_inter