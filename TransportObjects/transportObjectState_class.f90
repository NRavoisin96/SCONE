module transportObjectState_class

  use errors_mod,    only : fatalError
  use numPrecision
  use publicObjects, only : transportObjectStateCoordUpdateData

  implicit none
  private

  ! Public procedures.
  public :: display, init, kill, preparePayload

  !!
  !!
  !!
  type, public :: buildTransportObjectStatePayload
    integer(shortInt)           :: geometryIdx = 0, lowestCellIdx = 0, lowestElementIdx = 0, materialIdx = 0, uniqueId = 0
    real(defReal)               :: time = ZERO, weight = ZERO
    real(defReal), dimension(3) :: rGlobal = ZERO, uGlobal = ZERO
  end type buildTransportObjectStatePayload

  !!
  !!
  !!
  type, public :: transportObjectState
    private
    integer(shortInt)           :: geometryIdx = 0, lowestCellIdx = 0, lowestElementIdx = 0, materialIdx = 0, uniqueId = 0
    real(defReal)               :: time = ZERO, weight = ONE
    real(defReal), dimension(3) :: rGlobal = ZERO, uGlobal = ZERO
  contains
    procedure          :: display
    procedure          :: getCellIdx
    procedure          :: getCoordsUpdateData
    procedure          :: getElementIdx
    procedure          :: getGeometryIdx
    generic            :: getGlobalDirection => getGlobalDirection_defReal, getGlobalDirection_defRealArray
    procedure, private :: getGlobalDirection_defReal
    procedure, private :: getGlobalDirection_defRealArray
    generic            :: getGlobalPosition => getGlobalPosition_defReal, getGlobalPosition_defRealArray
    procedure, private :: getGlobalPosition_defReal
    procedure, private :: getGlobalPosition_defRealArray
    procedure          :: getMaterialIdx
    procedure          :: getUniqueId
    procedure          :: getWeight
    procedure          :: init
    procedure          :: kill
    procedure          :: preparePayload
    procedure          :: setGeometryIdx
    procedure          :: setGlobalDirection
    generic            :: setGlobalPosition => setGlobalPosition_defReal, setGlobalPosition_defRealArray
    procedure, private :: setGlobalPosition_defReal
    procedure, private :: setGlobalPosition_defRealArray
    procedure          :: setLowestCellIdx
    procedure          :: setLowestElementIdx
    procedure          :: setMaterialIdx
    procedure          :: setTime
    procedure          :: setUniqueId
    procedure          :: setWeight
    procedure          :: updateFromCoords
  end type transportObjectState

  !!
  !!
  !!
  type, public :: transportObjectStateBox
    class(transportObjectState), pointer :: ptr => null()
  end type transportObjectStateBox

contains
  !!
  !!
  !!
  subroutine display(self)
    class(transportObjectState), intent(in) :: self

    print *, 'Position: ', self % rGlobal
    print *, 'Direction: ', self % uGlobal
    print *, 'Weight: ', self % weight
    print *, 'Time: ', self % time
    print *, 'Geometry index: ', self % geometryIdx
    print *, 'Cell index: ', self % lowestCellIdx
    print *, 'Element index: ', self % lowestElementIdx
    print *, 'Unique id: ', self % uniqueId
    print *, 'Material index: ', self % materialIdx

  end subroutine display

  !!
  !!
  !!
  elemental function getCellIdx(self) result(cellIdx)
    class(transportObjectState), intent(in) :: self
    integer(shortInt)                       :: cellIdx

    cellIdx = self % lowestCellIdx

  end function getCellIdx

  !!
  !!
  !!
  elemental function getCoordsUpdateData(self) result(data)
    class(transportObjectState), intent(in)   :: self
    type(transportObjectStateCoordUpdateData) :: data

    data % geometryIdx = self % geometryIdx
    data % lowestCellIdx = self % lowestCellIdx
    data % lowestElementIdx = self % lowestElementIdx
    data % materialIdx = self % materialIdx
    data % uniqueId = self % uniqueId
    data % rGlobal = self % rGlobal
    data % uGlobal = self % uGlobal

  end function getCoordsUpdateData

  !!
  !!
  !!
  elemental function getElementIdx(self) result(elementIdx)
    class(transportObjectState), intent(in) :: self
    integer(shortInt)                       :: elementIdx

    elementIdx = self % lowestElementIdx

  end function getElementIdx

  !!
  !!
  !!
  elemental function getGeometryIdx(self) result(geometryIdx)
    class(transportObjectState), intent(in) :: self
    integer(shortInt)                       :: geometryIdx

    geometryIdx = self % geometryIdx

  end function getGeometryIdx

  !!
  !!
  !!
  elemental function getGlobalDirection_defReal(self, dimension) result(uGlobal)
    class(transportObjectState), intent(in) :: self
    integer(shortInt), intent(in)           :: dimension
    real(defReal)                           :: uGlobal

    uGlobal = self % uGlobal(dimension)

  end function getGlobalDirection_defReal

  !!
  !!
  !!
  pure function getGlobalDirection_defRealArray(self) result(uGlobal)
    class(transportObjectState), intent(in) :: self
    real(defReal), dimension(3)             :: uGlobal

    uGlobal = self % uGlobal

  end function getGlobalDirection_defRealArray

  !!
  !!
  !!
  elemental function getGlobalPosition_defReal(self, dimension) result(rGlobal)
    class(transportObjectState), intent(in) :: self
    integer(shortInt), intent(in)           :: dimension
    real(defReal)                           :: rGlobal

    rGlobal = self % rGlobal(dimension)

  end function getGlobalPosition_defReal

  !!
  !!
  !!
  pure function getGlobalPosition_defRealArray(self) result(rGlobal)
    class(transportObjectState),  intent(in) :: self
    real(defReal), dimension(3)              :: rGlobal

    rGlobal = self % rGlobal

  end function getGlobalPosition_defRealArray

  !!
  !!
  !!
  elemental function getMaterialIdx(self) result(materialIdx)
    class(transportObjectState), intent(in) :: self
    integer(shortInt)                       :: materialIdx

    materialIdx = self % materialIdx

  end function getMaterialIdx

  !!
  !!
  !!
  elemental function getUniqueId(self) result(uniqueId)
    class(transportObjectState), intent(in) :: self
    integer(shortInt)                       :: uniqueId

    uniqueId = self % uniqueId

  end function getUniqueId

  !!
  !!
  !!
  elemental function getWeight(self) result(weight)
    class(transportObjectState), intent(in) :: self
    real(defReal)                           :: weight

    weight = self % weight

  end function getWeight

  !!
  !!
  !!
  subroutine init(self, payload)
    class(transportObjectState), intent(inout)          :: self
    class(buildTransportObjectStatePayload), intent(in) :: payload

    self % geometryIdx = payload % geometryIdx
    self % lowestCellIdx = payload % lowestCellIdx
    self % lowestElementIdx = payload % lowestElementIdx
    self % materialIdx = payload % materialIdx
    self % uniqueId = payload % uniqueId
    self % time = payload % time
    self % weight = payload % weight
    self % rGlobal = payload % rGlobal
    self % uGlobal = payload % uGlobal

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(transportObjectState), intent(inout) :: self

    ! Local.
    self % geometryIdx = 0
    self % lowestCellIdx = 0
    self % lowestElementIdx = 0
    self % materialIdx = 0
    self % uniqueId = 0
    self % time = ZERO
    self % weight = ZERO
    self % rGlobal = ZERO
    self % uGlobal = ZERO

  end subroutine kill

  !!
  !!
  !!
  subroutine preparePayload(self, payload)
    class(transportObjectState), intent(in)                :: self
    class(buildTransportObjectStatePayload), intent(inout) :: payload

    ! Copy everything into payload.
    payload % geometryIdx = self % geometryIdx
    payload % lowestCellIdx = self % lowestCellIdx
    payload % lowestElementIdx = self % lowestElementIdx
    payload % materialIdx = self % materialIdx
    payload % uniqueId = self % uniqueId
    payload % time = self % time
    payload % weight = self % weight
    payload % rGlobal = self % rGlobal
    payload % uGlobal = self % uGlobal

  end subroutine preparePayload

  !!
  !!
  !!
  elemental subroutine setGeometryIdx(self, geometryIdx)
    class(transportObjectState), intent(inout) :: self
    integer(shortInt), intent(in)              :: geometryIdx

    self % geometryIdx = geometryIdx

  end subroutine setGeometryIdx

  !!
  !!
  !!
  pure subroutine setGlobalDirection(self, uGlobal)
    class(transportObjectState), intent(inout)      :: self
    real(defReal), dimension(3), target, intent(in) :: uGlobal

    self % uGlobal = uGlobal

  end subroutine setGlobalDirection

  !!
  !!
  !!
  elemental subroutine setGlobalPosition_defReal(self, dimension, rGlobal)
    class(transportObjectState), intent(inout) :: self
    integer(shortInt), intent(in)              :: dimension
    real(defReal), intent(in)                  :: rGlobal

    self % rGlobal(dimension) = rGlobal

  end subroutine setGlobalPosition_defReal

  !!
  !!
  !!
  pure subroutine setGlobalPosition_defRealArray(self, rGlobal)
    class(transportObjectState), intent(inout)      :: self
    real(defReal), dimension(3), target, intent(in) :: rGlobal

    self % rGlobal = rGlobal

  end subroutine setGlobalPosition_defRealArray

  !!
  !!
  !!
  elemental subroutine setLowestCellIdx(self, lowestCellIdx)
    class(transportObjectState), intent(inout) :: self
    integer(shortInt), intent(in)              :: lowestCellIdx

    self % lowestCellIdx = lowestCellIdx

  end subroutine setLowestCellIdx

  !!
  !!
  !!
  elemental subroutine setLowestElementIdx(self, lowestElementIdx)
    class(transportObjectState), intent(inout) :: self
    integer(shortInt), intent(in)              :: lowestElementIdx

    self % lowestElementIdx = lowestElementIdx

  end subroutine setLowestElementIdx

  !!
  !!
  !!
  elemental subroutine setMaterialIdx(self, materialIdx)
    class(transportObjectState), intent(inout) :: self
    integer(shortInt), intent(in)              :: materialIdx

    self % materialIdx = materialIdx

  end subroutine setMaterialIdx

  !!
  !!
  !!
  elemental subroutine setTime(self, time)
    class(transportObjectState), intent(inout) :: self
    real(defReal), intent(in)                  :: time

    self % time = time

  end subroutine setTime

  !!
  !!
  !!
  elemental subroutine setUniqueId(self, uniqueId)
    class(transportObjectState), intent(inout) :: self
    integer(shortInt), intent(in)              :: uniqueId

    self % uniqueId = uniqueId

  end subroutine setUniqueId

  !!
  !!
  !!
  elemental subroutine setWeight(self, weight)
    class(transportObjectState), intent(inout) :: self
    real(defReal), intent(in)                  :: weight

    self % weight = weight

  end subroutine setWeight

  !!
  !!
  !!
  elemental subroutine updateFromCoords(self, data)
    class(transportObjectState), intent(inout)            :: self
    type(transportObjectStateCoordUpdateData), intent(in) :: data

    ! Copy everything from data.
    self % geometryIdx = data % geometryIdx
    self % lowestCellIdx = data % lowestCellIdx
    self % lowestElementIdx = data % lowestElementIdx
    self % materialIdx = data % materialIdx
    self % uniqueId = data % uniqueId
    self % rGlobal = data % rGlobal
    self % uGlobal = data % uGlobal

  end subroutine updateFromCoords

end module transportObjectState_class