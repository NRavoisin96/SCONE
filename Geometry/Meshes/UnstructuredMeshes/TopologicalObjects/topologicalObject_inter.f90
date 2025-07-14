module topologicalObject_inter

  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
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

  !!
  !!
  !!
  type, public, abstract :: topologicalObject
    private
    integer(shortInt)    :: idx = 0
  contains
    ! Build procedures.
    procedure, non_overridable                  :: setIdx
    ! Runtime procedures.
    procedure(distanceSquared), deferred        :: distanceSquared
    procedure(getBoundingBoxBounds), deferred   :: getBoundingBoxBounds
    procedure(getCentroid), deferred            :: getCentroid
    procedure(getElements), deferred            :: getElements
    procedure, non_overridable                  :: getIdx
    generic                                     :: intersects => intersects_BoundingBox
    procedure(intersects_BoundingBox), deferred :: intersects_BoundingBox
    procedure                                   :: kill
  end type topologicalObject

  abstract interface
    !!
    !!
    !!
    function distanceSquared(self, r) result(dSquared)
      import                                  :: defReal, topologicalObject
      class(topologicalObject), intent(in)    :: self
      real(defReal), dimension(3), intent(in) :: r
      real(defReal)                           :: dSquared
    end function distanceSquared

    !!
    !!
    !!
    pure function getBoundingBoxBounds(self) result(bounds)
      import                               :: defReal, topologicalObject
      class(topologicalObject), intent(in) :: self
      real(defReal), dimension(3, 2)       :: bounds
    end function getBoundingBoxBounds

    !!
    !!
    !!
    pure function getCentroid(self) result(centroid)
      import                               :: defReal, topologicalObject
      class(topologicalObject), intent(in) :: self
      real(defReal), dimension(3)          :: centroid
    end function getCentroid

    !!
    !!
    !!
    function getElements(self) result(elements)
      import                                                :: topologicalObject, topologicalObjectBox
      class(topologicalObject), target, intent(in)          :: self
      type(topologicalObjectBox), dimension(:), allocatable :: elements
    end function getElements

    !!
    !!
    !!
    elemental function intersects_BoundingBox(self, boundingBox) result(doesIt)
      import                                   :: axisAlignedBoundingBox, defBool, topologicalObject
      class(topologicalObject), intent(in)     :: self
      type(axisAlignedBoundingBox), intent(in) :: boundingBox
      logical(defBool)                         :: doesIt
    end function intersects_BoundingBox

  end interface

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