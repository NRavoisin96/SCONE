module topologicalObject_inter

  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use numPrecision
  use publicObjects,                only : intersectionTestPayload, intersectionTestResult

  implicit none
  private

  ! Public procedures.
  public :: getUniqueSharingElements, init, kill

  !!
  !!
  !!
  type, public :: buildTopologicalObjectPayload
    integer(shortInt)         :: idx = 0
  end type buildTopologicalObjectPayload

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
    logical(defBool)     :: isActive = .true.
  contains
    ! Build procedures.
    procedure(build), deferred                  :: build
    procedure                                   :: init
    procedure                                   :: setIdx
    ! Runtime procedures.
    procedure, non_overridable                  :: deactivate
    procedure(distanceSquared), deferred        :: distanceSquared
    procedure(getBoundingBoxBounds), deferred   :: getBoundingBoxBounds
    procedure(getCentroid), deferred            :: getCentroid
    procedure, non_overridable                  :: getIdx
    procedure, non_overridable                  :: getIsActive
    procedure(getSharingElements), deferred     :: getSharingElements
    procedure(intersects_BoundingBox), deferred :: intersects_BoundingBox
    procedure                                   :: kill
  end type topologicalObject

  abstract interface
    !!
    !!
    !!
    subroutine build(self, payload)
      import                                              :: buildTopologicalObjectPayload, topologicalObject
      class(topologicalObject), intent(inout)             :: self
      class(buildTopologicalObjectPayload), intent(inout) :: payload
    end subroutine build

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
    function getSharingElements(self) result(sharingElements)
      import                                                :: topologicalObject, topologicalObjectBox
      class(topologicalObject), target, intent(in)          :: self
      type(topologicalObjectBox), dimension(:), allocatable :: sharingElements
    end function getSharingElements

    !!
    !!
    !!
    subroutine intersects_BoundingBox(self, boundingBox, doesIt)
      import                                   :: axisAlignedBoundingBox, defBool, topologicalObject
      class(topologicalObject), intent(in)     :: self
      type(axisAlignedBoundingBox), intent(in) :: boundingBox
      logical(defBool), intent(out)            :: doesIt
    end subroutine intersects_BoundingBox

  end interface

contains
  !!
  !!
  !!
  elemental subroutine deactivate(self)
    class(topologicalObject), intent(inout) :: self

    self % isActive = .false.

  end subroutine deactivate

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
  elemental function getIsActive(self) result(isActive)
    class(topologicalObject), intent(in) :: self
    logical(defBool)                     :: isActive

    isActive = self % isActive

  end function getIsActive

  !!
  !!
  !!
  function getUniqueSharingElements(objects) result(uniqueSharingElements)
    type(topologicalObjectBox), dimension(:), intent(in)  :: objects
    type(topologicalObjectBox), dimension(:), allocatable :: sharingElements, tempUniqueSharingElements,&
                                                             uniqueSharingElements
    integer(shortInt)                                     :: i, j, k, nUniqueSharingElements
    logical(defBool)                                      :: alreadyFound

    do i = 1, size(objects)
      sharingElements = objects(i) % ptr % getSharingElements()
      do j = 1, size(sharingElements)
        alreadyFound = .false.
        if (allocated(uniqueSharingElements)) then
          do k = 1, size(uniqueSharingElements)
            if (associated(sharingElements(j) % ptr, uniqueSharingElements(k) % ptr)) then
              alreadyFound = .true.
              exit

            end if

          end do

        end if

        if (.not. alreadyFound) then
          if (allocated(uniqueSharingElements)) then
            nUniqueSharingElements = size(uniqueSharingElements)
            allocate(tempUniqueSharingElements(nUniqueSharingElements + 1))
            tempUniqueSharingElements(1:nUniqueSharingElements) = uniqueSharingElements
            tempUniqueSharingElements(nUniqueSharingElements + 1) = sharingElements(j)
            call move_alloc(tempUniqueSharingElements, uniqueSharingElements)

          else
            allocate(uniqueSharingElements(1))
            uniqueSharingElements(1) = sharingElements(j)

          end if

        end if

      end do

    end do

  end function getUniqueSharingElements

  !!
  !!
  !!
  subroutine init(self, payload)
    class(topologicalObject), intent(inout)             :: self
    class(buildTopologicalObjectPayload), intent(inout) :: payload

    ! Set index from payload then build.
    self % idx = payload % idx
    call self % build(payload)

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(topologicalObject), intent(inout) :: self

    ! Local.
    self % idx = 0
    self % isActive = .true.

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