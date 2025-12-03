module physicalPropertyLaw_inter

  use dictionary_class, only : dictionary
  use numPrecision

  implicit none
  private

  ! Public procedures.
  public :: kill

  !!
  !!
  !!
  type, public, abstract :: physicalPropertyLaw
    private
  contains
    procedure(computeProperty), deferred :: computeProperty
    procedure(init), deferred            :: init
    procedure                            :: kill
  end type physicalPropertyLaw

  abstract interface
    !!
    !!
    !!
    function computeProperty(self, stateVariable) result(property)
      import                                 :: defReal, physicalPropertyLaw
      class(physicalPropertyLaw), intent(in) :: self
      real(defReal), intent(in)              :: stateVariable
      real(defReal)                          :: property
    end function computeProperty

    !!
    !!
    !!
    subroutine init(self, dict)
      import                                    :: dictionary, physicalPropertyLaw
      class(physicalPropertyLaw), intent(inout) :: self
      type(dictionary), intent(in)              :: dict
    end subroutine init

  end interface

contains
  !!
  !!
  !!
  elemental subroutine kill(self)
    class(physicalPropertyLaw), intent(inout) :: self

    ! Do nothing by default.

  end subroutine kill

end module physicalPropertyLaw_inter