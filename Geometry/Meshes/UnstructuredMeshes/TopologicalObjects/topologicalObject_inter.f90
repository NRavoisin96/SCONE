module topologicalObject_inter

  use numPrecision

  implicit none
  private

  ! Extendable procedures.
  public :: kill

  !!
  !! Small, local container to store polymorphic topological objects in an array.
  !!
  !! Public members:
  !!   ptr  -> Pointer to the topological object.
  !!
  type, public :: topologicalObjectBox
    class(topologicalObject), pointer :: ptr => null()
  end type topologicalObjectBox

  type, public, abstract :: topologicalObject
    private
    integer(shortInt)    :: idx = 0
  contains
    procedure, non_overridable :: getIdx
    procedure                  :: kill
    procedure, non_overridable :: setIdx
  end type topologicalObject

contains
  !!
  !!
  !!
  elemental function getIdx(self) result(idx)
    class(topologicalObject), intent(in) :: self
    integer(shortInt)                    :: idx

    idx = self % idx

  end function getIdx

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(topologicalObject), intent(inout) :: self

    ! Local.
    self % idx = 0

  end subroutine kill

  !!
  !!
  !!
  elemental subroutine setIdx(self, idx)
    class(topologicalObject), intent(inout) :: self
    integer(shortInt), intent(in)           :: idx

    self % idx = idx

  end subroutine setIdx

end module topologicalObject_inter