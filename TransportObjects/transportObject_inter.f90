module transportObject_inter

  use coordList_class,            only : coordList
  use errors_mod,                 only : fatalError
  use numPrecision
  use RNG_class,                  only : RNG
  use transportObjectState_class, only : buildTransportObjectStatePayload, transportObjectState

  implicit none
  private

  ! Public procedures.
  public :: init_base, kill

  !!
  !!
  !!
  type, public, abstract :: transportObject
    private
    class(RNG), pointer                      :: RNGPtr => null()
    class(transportObjectState), allocatable :: currentState, preTransitionState
    integer(shortInt)                        :: fate = 0
    logical(defBool)                         :: isDead = .false.
    real(defReal)                            :: initialWgt = ONE
    type(coordList)                          :: coords
  contains
    procedure(allocateState), deferred :: allocateState
    procedure                          :: copyCurrentState
    procedure                          :: displayCurrentState
    procedure                          :: generateDistance
    procedure                          :: getCellIdx
    procedure                          :: getCoordsLevel
    procedure                          :: getCoordsPtr
    procedure                          :: getCurrentStatePtr
    procedure                          :: getFate
    procedure                          :: getGeometryIdx
    procedure                          :: getGlobalDirection
    procedure                          :: getGlobalPosition
    procedure                          :: getIsDead
    procedure                          :: getLocalDirection
    procedure                          :: getLocalPosition
    procedure                          :: getMaterialIdx
    procedure                          :: getMeshIdx
    procedure                          :: getPreTransitionStatePtr
    procedure                          :: getRNGPtr
    procedure(getType), deferred       :: getType
    procedure                          :: getUniverseIdx
    procedure                          :: getWeight
    generic                            :: init => init_base, init_fromPayload, init_fromState
    procedure                          :: init_base
    procedure, private                 :: init_fromPayload
    procedure, private                 :: init_fromState
    procedure                          :: kill
    procedure                          :: moveGlobal
    procedure                          :: moveLocal
    procedure                          :: point
    procedure                          :: preparePayload
    procedure                          :: rotate
    procedure, non_overridable         :: savePreTransitionState
    procedure                          :: setFate
    procedure                          :: setGeometryIdx
    procedure                          :: setGlobalDirection
    procedure                          :: setGlobalPosition
    procedure                          :: setIsDead
    procedure                          :: setMaterialIdx
    procedure                          :: setRNGPtr
    procedure                          :: setWeight
    procedure                          :: strideRNG
    procedure, private                 :: synchroniseCurrentStateWithCoords
    procedure                          :: takeAboveGeometry
    procedure                          :: teleport
    procedure                          :: updateAndGetCurrentStatePtr
  end type transportObject

  abstract interface
    !!
    !!
    !!
    subroutine allocateState(self, state)
      import                                                :: transportObject, transportObjectState
      class(transportObject), intent(inout)                 :: self
      class(transportObjectState), allocatable, intent(out) :: state
    end subroutine allocateState

    !!
    !!
    !!
    elemental function getType(self) result(type)
      import                             :: shortInt, transportObject
      class(transportObject), intent(in) :: self
      integer(shortInt)                  :: type
    end function getType

  end interface

contains
  !!
  !!
  !!
  function copyCurrentState(self) result(currentStateCopy)
    class(transportObject), intent(in)       :: self
    class(transportObjectState), allocatable :: currentStateCopy

    currentStateCopy = self % currentState
    call currentStateCopy % updateFromCoords(self % coords % getTransportObjectStateUpdateData())

  end function copyCurrentState

  !!
  !!
  !!
  subroutine generateDistance(self, mult, distance)
    class(transportObject), intent(inout) :: self
    real(defReal), intent(in)             :: mult
    real(defReal), intent(out)            :: distance

    call self % RNGPtr % generateDistance(mult, distance)

  end subroutine generateDistance

  !!
  !!
  !!
  subroutine displayCurrentState(self)
    class(transportObject), intent(in) :: self

    call self % currentState % display()

  end subroutine displayCurrentState

  !!
  !!
  !!
  elemental function getCellIdx(self, l) result(cellIdx)
    class(transportObject), intent(in)      :: self
    integer(shortInt), intent(in), optional :: l
    integer(shortInt)                       :: cellIdx

    cellIdx = self % coords % getCellIdx(self % getCoordsLevel(l))

  end function getCellIdx

  !!
  !!
  !!
  elemental function getCoordsLevel(self, l) result(lvl)
    class(transportObject), intent(in)      :: self
    integer(shortInt), intent(in), optional :: l
    integer(shortInt)                       :: lvl

    if (present(l)) then
      lvl = l

    else
      lvl = self % coords % getNesting()

    end if

  end function getCoordsLevel

  !!
  !!
  !!
  function getCoordsPtr(self) result(coordsPtr)
    class(transportObject), target, intent(in) :: self
    type(coordList), pointer                   :: coordsPtr

    coordsPtr => self % coords

  end function getCoordsPtr

  !!
  !!
  !!
  function getCurrentStatePtr(self) result(currentStatePtr)
    class(transportObject), target, intent(in) :: self
    class(transportObjectState), pointer       :: currentStatePtr

    currentStatePtr => self % currentState

  end function getCurrentStatePtr

  !!
  !!
  !!
  elemental function getFate(self) result(fate)
    class(transportObject), intent(in) :: self
    integer(shortInt)                  :: fate

    fate = self % fate

  end function getFate

  !!
  !!
  !!
  elemental function getGeometryIdx(self) result(geometryIdx)
    class(transportObject), intent(in) :: self
    integer(shortInt)                  :: geometryIdx

    geometryIdx = self % coords % getGeometryIdx()

  end function getGeometryIdx

  !!
  !!
  !!
  pure function getGlobalDirection(self) result(uGlobal)
    class(transportObject), intent(in) :: self
    real(defReal), dimension(3)        :: uGlobal

    uGlobal = self % coords % getDirection(1)

  end function getGlobalDirection

  !!
  !!
  !!
  pure function getGlobalPosition(self) result(rGlobal)
    class(transportObject), intent(in) :: self
    real(defReal), dimension(3)        :: rGlobal

    rGlobal = self % coords % getPosition(1)

  end function getGlobalPosition

  !!
  !!
  !!
  elemental function getIsDead(self) result(isDead)
    class(transportObject), intent(in) :: self
    logical(defBool)                   :: isDead

    isDead = self % isDead

  end function getIsDead

  !!
  !!
  !!
  pure function getLocalDirection(self, l) result(uLocal)
    class(transportObject), intent(in)      :: self
    integer(shortInt), intent(in), optional :: l
    real(defReal), dimension(3)             :: uLocal

    uLocal = self % coords % getDirection(self % getCoordsLevel(l))

  end function getLocalDirection

  !!
  !!
  !!
  pure function getLocalPosition(self, l) result(rLocal)
    class(transportObject), intent(in)      :: self
    integer(shortInt), intent(in), optional :: l
    real(defReal), dimension(3)             :: rLocal

    rLocal = self % coords % getPosition(self % getCoordsLevel(l))

  end function getLocalPosition

  !!
  !!
  !!
  elemental function getMaterialIdx(self) result(materialIdx)
    class(transportObject), intent(in) :: self
    integer(shortInt)                  :: materialIdx

    materialIdx = self % coords % getMaterialIdx()

  end function getMaterialIdx

  !!
  !!
  !!
  elemental function getMeshIdx(self, l) result(meshIdx)
    class(transportObject), intent(in)      :: self
    integer(shortInt), intent(in), optional :: l
    integer(shortInt)                       :: meshIdx

    meshIdx = self % coords % getMeshIdx(self % getCoordsLevel(l))

  end function getMeshIdx

  !!
  !!
  !!
  function getPreTransitionStatePtr(self) result(preTransitionStatePtr)
    class(transportObject), target, intent(in) :: self
    class(transportObjectState), pointer       :: preTransitionStatePtr

    preTransitionStatePtr => self % preTransitionState

  end function getPreTransitionStatePtr

  !!
  !!
  !!
  function getRNGPtr(self) result(RNGPtr)
    class(transportObject), target, intent(in) :: self
    class(RNG), pointer                        :: RNGPtr

    RNGPtr => self % RNGPtr

  end function getRNGPtr

  !!
  !!
  !!
  elemental function getUniverseIdx(self, l) result(universeIdx)
    class(transportObject), intent(in)      :: self
    integer(shortInt), intent(in), optional :: l
    integer(shortInt)                       :: universeIdx

    universeIdx = self % coords % getUniverseIdx(self % getCoordsLevel(l))

  end function getUniverseIdx

  !!
  !!
  !!
  elemental function getWeight(self) result(weight)
    class(transportObject), intent(in) :: self
    real(defReal)                      :: weight

    weight = self % currentState % getWeight()

  end function getWeight

  !!
  !!
  !!
  subroutine init_base(self)
    class(transportObject), intent(inout) :: self

    ! Allocate memory for states.
    call self % allocateState(self % currentState)
    call self % allocateState(self % preTransitionState)

  end subroutine init_base

  !!
  !!
  !!
  subroutine init_fromPayload(self, payload)
    class(transportObject), intent(inout)               :: self
    class(buildTransportObjectStatePayload), intent(in) :: payload

    ! Initialise base attributes.
    call self % init_base()

    ! Build from payload.
    call self % coords % init(payload % rGlobal, payload % uGlobal)
    self % initialWgt = payload % weight
    call self % currentState % init(payload)

  end subroutine init_fromPayload

  !!
  !!
  !!
  subroutine init_fromState(self, state)
    class(transportObject), intent(inout)   :: self
    class(transportObjectState), intent(in) :: state

    ! Initialise base attributes.
    call self % init_base()

    ! Set current state.
    self % currentState = state
    self % initialWgt = state % getWeight()

    ! Set informations in coordList.
    call self % coords % init(state % getCoordsUpdateData())

  end subroutine init_fromState

  !!
  !!
  !!
  subroutine kill(self)
    class(transportObject), intent(inout) :: self

    ! Local.
    self % RNGPtr => null()
    if (allocated(self % currentState)) then
      call self % currentState % kill()
      deallocate(self % currentState)

    end if
    if (allocated(self % preTransitionState)) then
      call self % preTransitionState % kill()
      deallocate(self % preTransitionState)

    end if
    self % fate = 0
    self % isDead = .false.
    self % initialWgt = ZERO
    call self % coords % kill()

  end subroutine kill

  !!
  !!
  !!
  elemental subroutine moveGlobal(self, d)
    class(transportObject), intent(inout) :: self
    real(defReal), intent(in)             :: d

    call self % coords % moveGlobal(d)

  end subroutine moveGlobal

  !!
  !!
  !!
  subroutine moveLocal(self, lvl, d)
    class(transportObject), intent(inout) :: self
    integer(shortInt), intent(in)         :: lvl
    real(defReal), intent(in)             :: d

    call self % coords % moveLocal(d, lvl)

  end subroutine moveLocal

  !!
  !!
  !!
  pure subroutine point(self, u)
    class(transportObject), intent(inout)   :: self
    real(defReal), dimension(3), intent(in) :: u

    call self % coords % assignDirection(u)

  end subroutine point

  !!
  !!
  !!
  subroutine preparePayload(self, payload)
    class(transportObject), intent(inout)                  :: self
    class(buildTransportObjectStatePayload), intent(inout) :: payload

    ! Update current state then copy everything into payload.
    call self % currentState % updateFromCoords(self % coords % getTransportObjectStateUpdateData())
    call self % currentState % preparePayload(payload)

  end subroutine preparePayload

  !!
  !!
  !!
  elemental subroutine rotate(self, mu, phi)
    class(transportObject), intent(inout) :: self
    real(defReal), intent(in)             :: mu, phi

    call self % coords % rotate(mu, phi)

  end subroutine rotate

  !!
  !!
  !!
  subroutine savePreTransitionState(self)
    class(transportObject), intent(inout) :: self

    self % preTransitionState = self % copyCurrentState()

  end subroutine savePreTransitionState

  !!
  !!
  !!
  elemental subroutine setFate(self, fate)
    class(transportObject), intent(inout) :: self
    integer(shortInt), intent(in)         :: fate

    self % fate = fate

  end subroutine setFate

  !!
  !!
  !!
  elemental subroutine setGeometryIdx(self, geometryIdx)
    class(transportObject), intent(inout) :: self
    integer(shortInt), intent(in)         :: geometryIdx

    call self % coords % setGeometryIdx(geometryIdx)

  end subroutine setGeometryIdx

  !!
  !!
  !!
  pure subroutine setGlobalDirection(self, uGlobal)
    class(transportObject), intent(inout)   :: self
    real(defReal), dimension(3), intent(in) :: uGlobal

    ! Set in coordinates and current state.
    call self % coords % setDirection(uGlobal, 1)

  end subroutine setGlobalDirection

  !!
  !!
  !!
  pure subroutine setGlobalPosition(self, rGlobal)
    class(transportObject), intent(inout)   :: self
    real(defReal), dimension(3), intent(in) :: rGlobal

    ! Set in coordinates and current state.
    call self % coords % setPosition(rGlobal, 1)

  end subroutine setGlobalPosition

  !!
  !!
  !!
  elemental subroutine setIsDead(self, isDead)
    class(transportObject), intent(inout) :: self
    logical(defBool), intent(in)          :: isDead

    self % isDead = isDead

  end subroutine setIsDead

  !!
  !!
  !!
  elemental subroutine setMaterialIdx(self, materialIdx)
    class(transportObject), intent(inout) :: self
    integer(shortInt), intent(in)         :: materialIdx

    call self % coords % setMaterialIdx(materialIdx)

  end subroutine setMaterialIdx

  !!
  !!
  !!
  subroutine setRNGPtr(self, rand)
    class(transportObject), intent(inout) :: self
    type(RNG), target, intent(in)         :: rand

    self % RNGPtr => rand

  end subroutine setRNGPtr

  !!
  !!
  !!
  elemental subroutine setWeight(self, weight)
    class(transportObject), intent(inout) :: self
    real(defReal), intent(in)             :: weight

    call self % currentState % setWeight(weight)

  end subroutine setWeight

  !!
  !!
  !!
  subroutine strideRNG(self, n)
    class(transportObject), intent(inout) :: self
    integer(shortInt), intent(in)         :: n

    call self % RNGPtr % stride(n)

  end subroutine strideRNG

  !!
  !!
  !!
  elemental subroutine synchroniseCurrentStateWithCoords(self)
    class(transportObject), intent(inout) :: self

    call self % currentState % updateFromCoords(self % coords % getTransportObjectStateUpdateData())

  end subroutine synchroniseCurrentStateWithCoords

  !!
  !!
  !!
  elemental subroutine takeAboveGeometry(self)
    class(transportObject), intent(inout) :: self

    call self % coords % takeAboveGeom()

  end subroutine takeAboveGeometry

  !!
  !!
  !!
  pure subroutine teleport(self, r)
    class(transportObject), intent(inout)   :: self
    real(defReal), dimension(3), intent(in) :: r

    call self % coords % assignPosition(r)

  end subroutine teleport

  !!
  !!
  !!
  function updateAndGetCurrentStatePtr(self) result(currentStatePtr)
    class(transportObject), target, intent(in) :: self
    class(transportObjectState), pointer       :: currentStatePtr

    ! Update then get pointer.
    currentStatePtr => self % currentState
    call currentStatePtr % updateFromCoords(self % coords % getTransportObjectStateUpdateData())

  end function updateAndGetCurrentStatePtr

end module transportObject_inter