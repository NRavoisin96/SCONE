module extentTopologicalObject_inter

  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use genericProcedures,            only : fatalError
  use numPrecision
  use publicObjects,                only : intersectionTestPayload, intersectionTestResult
  use ratint
  use topologicalObject_inter,      only : buildTopologicalObjectPayload, init_super => init, kill_super => kill, &
                                           topologicalObject
  use universalVariables,           only : INF
  use vertex_class,                 only : vertexBox

  implicit none
  private

  ! Extendable procedures.
  public :: kill

  !!
  !!
  !!
  type, public, extends(buildTopologicalObjectPayload) :: buildExtentTopologicalObjectPayload
    real(defReal), dimension(3)                 :: centroid
    real(defReal), dimension(:, :), allocatable :: allCoords
    type(ratint_t), dimension(3)                :: rationalCentroid
    type(vertexBox), dimension(:), allocatable  :: vertices
  end type buildExtentTopologicalObjectPayload

  !!
  !!
  !!
  type, public, abstract, extends(topologicalObject) :: extentTopologicalObject
    private
    real(defReal), dimension(3)  :: centroid
    type(axisAlignedBoundingBox) :: boundingBox
    type(ratint_t), dimension(3) :: rationalCentroid
  contains
    ! Build procedures.
    procedure(connectComponents), deferred :: connectComponents
    procedure                              :: init
    procedure                              :: kill
    ! Runtime procedures.
    procedure                              :: getBoundingBoxBounds
    procedure                              :: getBoundingBoxPtr
    procedure                              :: getCentroid
    procedure                              :: getRationalCentroid
    generic                                :: intersectsBoundingBox => intersectsBoundingBox_BoundingBox
    procedure, private                     :: intersectsBoundingBox_BoundingBox
  end type extentTopologicalObject

  abstract interface
    !!
    !!
    !!
    subroutine connectComponents(self)
      import                                                :: extentTopologicalObject
      class(extentTopologicalObject), target, intent(inout) :: self
    end subroutine connectComponents

  end interface

contains
  !!
  !!
  !!
  pure function getBoundingBoxBounds(self) result(bounds)
    class(extentTopologicalObject), intent(in) :: self
    real(defReal), dimension(3, 2)             :: bounds

    bounds = self % boundingBox % getBounds()

  end function getBoundingBoxBounds

  !!
  !!
  !!
  function getBoundingBoxPtr(self) result(ptr)
    class(extentTopologicalObject), target, intent(in) :: self
    type(axisAlignedBoundingBox), pointer              :: ptr

    ptr => self % boundingBox

  end function getBoundingBoxPtr

  !!
  !!
  !!
  pure function getCentroid(self) result(centroid)
    class(extentTopologicalObject), intent(in) :: self
    real(defReal), dimension(3)                :: centroid

    centroid = self % centroid

  end function getCentroid

  !!
  !!
  !!
  pure function getRationalCentroid(self) result(rationalCentroid)
    class(extentTopologicalObject), intent(in) :: self
    type(ratint_t), dimension(3)               :: rationalCentroid

    rationalCentroid = self % rationalCentroid

  end function getRationalCentroid

  !!
  !!
  !!
  subroutine init(self, payload)
    class(extentTopologicalObject), intent(inout)       :: self
    class(buildTopologicalObjectPayload), intent(inout) :: payload
    class(buildExtentTopologicalObjectPayload), pointer :: payloadPtr
    character(*), parameter                             :: here = 'init (extentTopologicalObject_inter.f90)'

    ! Downcast payload to correct type.
    select type(ptr => payload)
      class is(buildExtentTopologicalObjectPayload)
        payloadPtr => ptr

        ! Initialise superclass and update connectivity.
        if (allocated(payloadPtr % allCoords)) deallocate(payloadPtr % allCoords)
        call init_super(self, payloadPtr)
        call self % connectComponents()

        ! Set centroid and compute bounding box.
        if (.not. allocated(payloadPtr % allCoords)) &
        call fatalError(here, 'Unallocated coordinates array for bounding box computation.')
        
        self % centroid = payloadPtr % centroid
        self % rationalCentroid = payloadPtr % rationalCentroid
        call self % boundingBox % computeBounds(payloadPtr % allCoords)

      class default
        call fatalError(here, 'Invalid payload type.')

    end select

  end subroutine

  !!
  !!
  !!
  elemental function intersectsBoundingBox_BoundingBox(self, boundingBox) result(doesIt)
    class(extentTopologicalObject), intent(in) :: self
    type(axisAlignedBoundingBox), intent(in)   :: boundingBox
    logical(defBool)                           :: doesIt

    call self % boundingBox % intersects(boundingBox, doesIt)

  end function intersectsBoundingBox_BoundingBox

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(extentTopologicalObject), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % centroid = ZERO
    call self % boundingBox % kill()

  end subroutine kill

end module extentTopologicalObject_inter