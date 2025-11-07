module physicalParticle_inter

  use errors_mod,                  only : fatalError
  use numPrecision
  use physicalParticleState_class, only : castPhysicalParticleStatePtr, physicalParticleState
  use transportObject_inter,       only : init_base_super => init_base, kill_super => kill, transportObject
  use transportObjectState_class,  only : transportObjectState

  implicit none
  private

  ! Public procedures.
  public :: castPhysicalParticlePtr, init_base, kill

  !!
  !!
  !!
  type, public, abstract, extends(transportObject) :: physicalParticle
    private
    class(transportObjectState), allocatable :: preCollisionState, preHistoryState, prePathState
    integer(shortInt)                        :: nCollisions = 0, nSplits = 0
    real(defReal)                            :: k_eff = ONE, timeMax = ZERO
  contains
    procedure                     :: getBroodId
    procedure                     :: getKEff
    procedure                     :: getPreCollisionStatePtr
    procedure                     :: getPreHistoryStatePtr
    procedure                     :: getPrePathStatePtr
    procedure(getSpeed), deferred :: getSpeed
    procedure                     :: getSplitsNumber
    procedure                     :: incrementCollisionsNumber
    procedure                     :: init_base
    procedure                     :: kill
    procedure, non_overridable    :: savePreCollisionState
    procedure, non_overridable    :: savePreHistoryState
    procedure, non_overridable    :: savePrePathState
    procedure                     :: setBroodId
    procedure                     :: setKEff
    procedure                     :: setMass
    procedure                     :: setSplitsNumber
  end type physicalParticle

  abstract interface
    !!
    !!
    !!
    function getSpeed(self) result(speed)
      import                              :: defReal, physicalParticle
      class(physicalParticle), intent(in) :: self
      real(defReal)                       :: speed
    end function getSpeed

  end interface

contains
  !!
  !!
  !!
  function castPhysicalParticlePtr(source, fatal) result(ptr)
    class(transportObject), intent(in)     :: source
    logical(defBool), intent(in), optional :: fatal
    class(physicalParticle), pointer       :: ptr
    logical(defBool)                       :: throwError
    character(*), parameter                :: here = 'castPhysicalParticlePtr (physicalParticle_inter.f90)'

    ! Downcast source to correct class.
    select type(temp => source)
      class is(physicalParticle)
        ptr => temp

      class default
        ptr => null()

    end select

    ! Throw error if requested.
    throwError = .false.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) call fatalError(here, "Unable to downcast source to class 'physicalParticle'.")

  end function castPhysicalParticlePtr

  !!
  !!
  !!
  function getBroodId(self) result(broodId)
    class(physicalParticle), intent(in)   :: self
    class(physicalParticleState), pointer :: physicalParticleStatePtr
    integer(shortInt)                     :: broodId

    physicalParticleStatePtr => castPhysicalParticleStatePtr(self % getCurrentStatePtr(), .true.)
    broodId = physicalParticleStatePtr % getBroodId()

  end function getBroodId

  !!
  !!
  !!
  elemental function getKEff(self) result(k_eff)
    class(physicalParticle), intent(in) :: self
    real(defReal)                       :: k_eff

    k_eff = self % k_eff

  end function getKEff

  !!
  !!
  !!
  function getPreCollisionStatePtr(self) result(preCollisionStatePtr)
    class(physicalParticle), target, intent(in)  :: self
    class(transportObjectState), pointer         :: preCollisionStatePtr

    preCollisionStatePtr => self % preCollisionState

  end function getPreCollisionStatePtr

  !!
  !!
  !!
  function getPreHistoryStatePtr(self) result(preHistoryStatePtr)
    class(physicalParticle), target, intent(in) :: self
    class(transportObjectState), pointer        :: preHistoryStatePtr

    preHistoryStatePtr => self % preHistoryState

  end function getPreHistoryStatePtr

  !!
  !!
  !!
  function getPrePathStatePtr(self) result(prePathStatePtr)
    class(physicalParticle), target, intent(in) :: self
    class(transportObjectState), pointer        :: prePathStatePtr

    prePathStatePtr => self % prePathState

  end function getPrePathStatePtr

  !!
  !!
  !!
  elemental function getSplitsNumber(self) result(nSplits)
    class(physicalParticle), intent(in) :: self
    integer(shortInt)                   :: nSplits

    nSplits = self % nSplits

  end function getSplitsNumber

  !!
  !!
  !!
  elemental subroutine incrementCollisionsNumber(self, nCollisions)
    class(physicalParticle), intent(inout) :: self
    integer(shortInt), intent(in)          :: nCollisions

    self % nCollisions = self % nCollisions + nCollisions

  end subroutine incrementCollisionsNumber

  !!
  !!
  !!
  subroutine init_base(self)
    class(physicalParticle), intent(inout) :: self

    ! Superclass.
    call init_base_super(self)

    ! Allocate physical states.
    call self % allocateState(self % preCollisionState)
    call self % allocateState(self % preHistoryState)
    call self % allocateState(self % prePathState)

  end subroutine init_base

  !!
  !!
  !!
  subroutine kill(self)
    class(physicalParticle), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    if (allocated(self % preCollisionState)) then
      call self % preCollisionState % kill()
      deallocate(self % preCollisionState)

    end if
    
    if (allocated(self % preHistoryState)) then
      call self % preHistoryState % kill()
      deallocate(self % preHistoryState)

    end if

    if (allocated(self % prePathState)) then
      call self % prePathState % kill()
      deallocate(self % prePathState)

    end if

    self % nCollisions = 0
    self % nSplits = 0
    self % k_eff = ONE
    self % timeMax = ZERO

  end subroutine kill

  !!
  !!
  !!
  subroutine savePreCollisionState(self)
    class(physicalParticle), intent(inout) :: self

    self % preCollisionState = self % copyCurrentState()

  end subroutine savePreCollisionState

  !!
  !!
  !!
  subroutine savePreHistoryState(self)
    class(physicalParticle), intent(inout) :: self

    self % preHistoryState = self % copyCurrentState()

  end subroutine savePreHistoryState

  !!
  !!
  !!
  subroutine savePrePathState(self)
    class(physicalParticle), intent(inout) :: self

    self % prePathState = self % copyCurrentState()

  end subroutine savePrePathState

  !!
  !!
  !!
  subroutine setBroodId(self, broodId)
    class(physicalParticle), intent(inout) :: self
    integer(shortInt), intent(in)          :: broodId
    class(physicalParticleState), pointer  :: physicalParticleStatePtr

    physicalParticleStatePtr => castPhysicalParticleStatePtr(self % getCurrentStatePtr(), .true.)
    call physicalParticleStatePtr % setBroodId(broodId)

  end subroutine setBroodId

  !!
  !!
  !!
  elemental subroutine setKEff(self, k_eff)
    class(physicalParticle), intent(inout) :: self
    real(defReal), intent(in)              :: k_eff

    self % k_eff = k_eff

  end subroutine setKEff

  !!
  !!
  !!
  subroutine setMass(self, mass)
    class(physicalParticle), intent(inout) :: self
    real(defReal), intent(in)              :: mass
    class(physicalParticleState), pointer  :: currentStatePtr

    currentStatePtr => castPhysicalParticleStatePtr(self % getCurrentStatePtr(), .true.)
    call currentStatePtr % setMass(mass)

  end subroutine setMass

  !!
  !!
  !!
  elemental subroutine setSplitsNumber(self, nSplits)
    class(physicalParticle), intent(inout) :: self
    integer(shortInt), intent(in)          :: nSplits

    self % nSplits = nSplits

  end subroutine setSplitsNumber

end module physicalParticle_inter