module collisionProcessor_inter

  use collisionData_class,         only : collisionData
  use dictionary_class,            only : dictionary
  use endfConstants
  use errors_mod,                  only : fatalError
  use genericProcedures,           only : numToChar
  use materialHandle_inter,        only : materialHandle
  use neutronMaterial_inter,       only : neutronMaterial
  use neutronXsPackages_class,     only : neutronMacroXSs
  use nuclearDatabase_inter,       only : nuclearDatabase
  use nuclideHandle_inter,         only : nuclideHandle
  use numPrecision
  use particleDungeon_class,       only : particleDungeon
  use physicalParticle_inter,      only : castPhysicalParticlePtr, physicalParticle
  use physicalParticleState_class, only : physicalParticleState
  use reactionHandle_inter,        only : reactionHandle
  use RNG_class,                   only : RNG
  use tallyCodes
  use tallyAdmin_class,            only : tallyAdmin
  use transportObject_inter,       only : transportObject
  use transportObjectState_class,  only : buildTransportObjectStatePayload

  implicit none
  private

  ! Public procedures.
  public :: fission, implicit, init, kill

  !!
  !! This is an abstract interface for all types of collision processing
  !!  -> This interface only deals with SCALAR processing of collisions
  !!  -> Note that it is NOT just a collection of deferred function. There is a master
  !!     non_overridable function that controls the flow of the calculation and calls a number
  !!     of customisable deferred function.
  !!  -> Interesed user can refer to http://www.gotw.ca/publications/mill18.htm for justification
  !!     for this approach.
  !!
  !! Public interface:
  !!   collide(p, tally, thisCycle, nextCycle) -> Given particle, tallyAdmin, and particleDungeons
  !!     for particles in this cycle, or next cycle performs particle collision. Sends pre and post
  !!     collision events report to the tally admin. Reports end of history if particle was absorbed.
  !!
  !!   init(dict) -> initialises collisionProcessor from a dictionary
  !!
  !!
  !! Customisable procedures or collision actions (implemented in subclasses):
  !!   sampleCollision -> should determine collision target and type. Sets approperiate data in
  !!                      collisionData package.
  !!   implicit        -> preforms any implicit treatment to be done before reaction processing,
  !!                      e.g. implicit fission site production
  !!   scatter         -> defines behaviour for a scattering event (macroscopic scatter or any
  !!                      nuclide scattering reaction (elastic, inealastic, campton, NXN etc. )
  !!   capture         -> defines behaviour for any capture reaction (e.g. gamma capture, photoelectric)
  !!   fission         -> defines bahaviour for any fission reaction
  !!   cutoffs         -> Any post collision implicit treatments i.e. energy cutoffs
  !!
  type, public, abstract :: collisionProcessor
    private
    ! Nuclear Data block pointers.
    class(neutronMaterial), pointer :: currentMaterialPtr => null()
    class(nuclearDatabase), pointer :: xsData => null()
  contains
    procedure(allocateCollisionData), deferred    :: allocateCollisionData
    procedure                                     :: capture
    procedure, non_overridable                    :: collide
    procedure                                     :: computeNumberOfSecondaryParticles
    procedure                                     :: cutoffs
    procedure                                     :: elastic
    procedure                                     :: fission
    procedure                                     :: getCurrentMaterialIsFissile
    procedure                                     :: getCurrentMaterialPtr
    procedure(getFissionHandlePtr), deferred      :: getFissionHandlePtr
    procedure                                     :: getMaterialPtr
    procedure                                     :: getNuclearDatabasePtr
    procedure                                     :: getNuclidePtr
    procedure                                     :: getReaction
    procedure                                     :: implicit
    procedure(inelastic), deferred                :: inelastic
    procedure                                     :: init
    procedure(getImplicitCondition), deferred     :: getImplicitCondition
    procedure                                     :: kill
    procedure(newPhysicalParticleState), deferred :: newPhysicalParticleState
    procedure(sampleCollision), deferred          :: sampleCollision
    procedure                                     :: setCurrentMaterialPtr
    procedure                                     :: setNuclearDatabasePtr
  end type collisionProcessor

  abstract interface
    !!
    !!
    !!
    subroutine allocateCollisionData(self, collDat)
      import                                         :: collisionData, collisionProcessor
      class(collisionProcessor), intent(in)          :: self
      class(collisionData), allocatable, intent(out) :: collDat
    end subroutine allocateCollisionData

    !!
    !!
    !!
    function getFissionHandlePtr(self, collDat) result(fissionHandlePtr)
      import                                :: collisionData, collisionProcessor, reactionHandle
      class(collisionProcessor), intent(in) :: self
      class(collisionData), intent(in)      :: collDat
      class(reactionHandle), pointer        :: fissionHandlePtr
    end function getFissionHandlePtr

    !!
    !!
    !!
    subroutine inelastic(self, collDat, p)
      import                                 :: collisionData, collisionProcessor, physicalParticle
      class(collisionProcessor), intent(in)  :: self
      class(collisionData), intent(inout)    :: collDat
      class(physicalParticle), intent(inout) :: p
    end subroutine inelastic

    !!
    !!
    !!
    elemental function getImplicitCondition(self) result(isIt)
      import                                :: collisionProcessor, defBool
      class(collisionProcessor), intent(in) :: self
      logical(defBool)                      :: isIt
    end function getImplicitCondition

    !!
    !!
    !!
    subroutine newPhysicalParticleState(self, collDat, fissionHandlePtr, payload, newState)
      import :: buildTransportObjectStatePayload, collisionData, collisionProcessor, physicalParticleState, reactionHandle
      class(collisionProcessor), intent(in)                  :: self
      class(collisionData), intent(in)                       :: collDat
      class(reactionHandle), pointer, intent(in)             :: fissionHandlePtr
      class(buildTransportObjectStatePayload), intent(inout) :: payload
      class(physicalParticleState), allocatable, intent(out) :: newState
    end subroutine newPhysicalParticleState

    !!
    !!
    !!
    subroutine sampleCollision(self, p, collDat)
      import                                   :: collisionProcessor, collisionData, physicalParticle
      class(collisionProcessor), intent(inout) :: self
      class(physicalParticle), intent(in)      :: p
      class(collisionData), intent(inout)      :: collDat
    end subroutine sampleCollision
    
  end interface

contains
  !!
  !!
  !!
  subroutine capture(self, p)
    class(collisionProcessor), intent(in)  :: self
    class(physicalParticle), intent(inout) :: p

    ! Kill particle by default.
    call p % setIsDead(.true.)

  end subroutine capture

  !!
  !! Generic flow of collision processing
  !!
  subroutine collide(self, object, tally, thisCycle, nextCycle)
    class(collisionProcessor), intent(inout) :: self
    class(transportObject), intent(inout)    :: object
    type(tallyAdmin), intent(inout)          :: tally
    class(particleDungeon), intent(inout)    :: thisCycle, nextCycle
    class(collisionData), allocatable        :: collDat
    class(physicalParticle), pointer         :: p
    logical(defBool)                         :: virtual
    character(*), parameter                  :: here = 'collide (collisionProcessor.f90)'

    ! Downcast transportObject to physicalParticle.
    p => castPhysicalParticlePtr(object, .true.)

    ! Allocate collisionData object to correct type then copy particle data into it.
    call self % allocateCollisionData(collDat)
    call p % prepareCollisionData(collDat)

    ! Choose collision nuclide and general type (Scatter, Capture or Fission)
    call self % sampleCollision(p, collDat)

    ! In case of a TMS rejection, set collision as virtual
    virtual = collDat % MT == noInteraction

    ! Report in-collision & save pre-collison state
    ! Note: the ordering must not be changed between feeding the particle to the tally
    ! and updating the particle's preCollision state, otherwise this may cause certain
    ! tallies (e.g., collisionProbability) to return dubious results
    call tally % reportInColl(p, virtual)
    call p % savePreCollisionState()

    ! Perform implicit treatment
    if (.not. virtual) then
      call self % implicit(self % getImplicitCondition(), collDat % initialWeight * collDat % sigma_tot, &
                           collDat, p, nextCycle, tally)
    
      ! Select physics to be processed based on MT number
      select case(collDat % MT)
        case(N_N_elastic, macroAllScatter)
          call self % elastic(collDat, p)

        case(N_N_inelastic, macroIEScatter)
          call self % inelastic(collDat, p)

        case(N_DISAP, macroCapture)
          call self % capture(p)

        case(N_FISSION, macroFission)
          call self % fission(collDat, p, nextCycle, tally)

        case default
          call fatalError(here, 'Unsupported MT number: '//numToChar(collDat % MT)//'.')

      end select

    end if

    ! Apply post collision implicit treatments if particle survived collision.
    if (.not. p % getIsDead()) call self % cutoffs(p, thisCycle, tally)

    ! Update particle collision counter
    call p % incrementCollisionsNumber(merge(0, 1, virtual))

    ! Report out-of-collision
    call tally % reportOutColl(p, collDat % MT, collDat % muL)

    ! Report end-of-history if particle was killed
    if (p % getIsDead()) then
      call p % setFate(ABS_FATE)
      call tally % reportHist(p)

    end if

  end subroutine collide

  !!
  !!
  !!
  subroutine computeNumberOfSecondaryParticles(self, factor, collDat, p)
    class(collisionProcessor), intent(in)  :: self
    real(defReal), intent(in)              :: factor
    class(collisionData), intent(inout)    :: collDat
    class(physicalParticle), intent(inout) :: p
    real(defReal)                          :: randomNumber

    ! Sample number of fission sites generated. Supports negative weights.
    call collDat % RNGPtr % generate(randomNumber)
    collDat % n = int(abs((collDat % weight * collDat % sigma_nuFiss) / (factor * collDat % k_eff)) + randomNumber, shortInt)
    collDat % implicitWeight = sign(collDat % initialWeight, collDat % weight)

  end subroutine computeNumberOfSecondaryParticles

  !!
  !!
  !!
  subroutine cutoffs(self, p, dungeon, tally)
    class(collisionProcessor), intent(in)  :: self
    class(physicalParticle), intent(inout) :: p
    type(particleDungeon), intent(inout)   :: dungeon
    type(tallyAdmin), intent(inout)        :: tally

    ! Do nothing by default.

  end subroutine cutoffs

  !!
  !!
  !!
  subroutine elastic(self, collDat, p)
    class(collisionProcessor), intent(in)  :: self
    class(collisionData), intent(inout)    :: collDat
    class(physicalParticle), intent(inout) :: p

    ! Do nothing by default.

  end subroutine elastic

  !!
  !!
  !!
  subroutine fission(self, collDat, p, dungeon, tally)
    class(collisionProcessor), intent(in)  :: self
    class(collisionData), intent(inout)    :: collDat
    class(physicalParticle), intent(inout) :: p
    type(particleDungeon), intent(inout)   :: dungeon
    type(tallyAdmin), intent(inout)        :: tally

    ! Kill particle by default.
    call p % setIsDead(.true.)

  end subroutine fission

  !!
  !!
  !!
  elemental function getCurrentMaterialIsFissile(self) result(isFissile)
    class(collisionProcessor), intent(in) :: self
    logical(defBool)                      :: isFissile

    isFissile = self % currentMaterialPtr % isFissile()

  end function getCurrentMaterialIsFissile

  !!
  !!
  !!
  function getCurrentMaterialPtr(self) result(currentMaterialPtr)
    class(collisionProcessor), intent(in) :: self
    class(neutronMaterial), pointer       :: currentMaterialPtr

    currentMaterialPtr => self % currentMaterialPtr

  end function getCurrentMaterialPtr

  !!
  !!
  !!
  function getMaterialPtr(self, materialIdx) result(materialPtr)
    class(collisionProcessor), intent(in) :: self
    integer(shortInt), intent(in)         :: materialIdx
    class(materialHandle), pointer        :: materialPtr

    materialPtr => self % xsData % getMaterial(materialIdx)

  end function getMaterialPtr

  !!
  !!
  !!
  function getNuclearDatabasePtr(self) result(xsDataPtr)
    class(collisionProcessor), intent(in) :: self
    class(nuclearDatabase), pointer       :: xsDataPtr

    xsDataPtr => self % xsData

  end function getNuclearDatabasePtr

  !!
  !!
  !!
  function getNuclidePtr(self, nuclideIdx) result(nuclidePtr)
    class(collisionProcessor), intent(in) :: self
    integer(shortInt), intent(in)         :: nuclideIdx
    class(nuclideHandle), pointer         :: nuclidePtr

    nuclidePtr => self % xsData % getNuclide(nuclideIdx)

  end function getNuclidePtr

  !!
  !!
  !!
  function getReaction(self, reactionChannel, idx) result(reaction)
    class(collisionProcessor), intent(in) :: self
    integer(shortInt), intent(in)         :: reactionChannel, idx
    class(reactionHandle), pointer        :: reaction

    reaction => self % xsData % getReaction(reactionChannel, idx)

  end function getReaction

  !!
  !!
  !!
  subroutine implicit(self, implicitCondition, factor, collDat, p, dungeon, tally)
    class(collisionProcessor), intent(in)                :: self
    logical(defBool), intent(in)                         :: implicitCondition
    real(defReal), intent(in)                            :: factor
    class(collisionData), intent(inout)                  :: collDat
    class(physicalParticle), intent(inout)               :: p
    type(particleDungeon), intent(inout)                 :: dungeon
    type(tallyAdmin), intent(inout)                      :: tally
    class(buildTransportObjectStatePayload), allocatable :: payload
    class(physicalParticleState), allocatable            :: newState
    class(reactionHandle), pointer                       :: fissionHandlePtr
    integer(shortInt)                                    :: i
    character(*), parameter                              :: HERE = 'implicit (collisionProcessor_inter.f90)'

    ! By default perform implicit fission.
    if (.not. implicitCondition) return
    
    ! Compute number of secondary particles generated.
    call self % computeNumberOfSecondaryParticles(factor, collDat, p)

    ! Add secondary particles to dungeon.
    if (collDat % n < 1) return
    
    ! Get pointer to fission reaction.
    fissionHandlePtr => self % getFissionHandlePtr(collDat)
    if (.not. associated(fissionHandlePtr)) call fatalError(HERE, 'Failed to retrieve fission reaction handle.')

    ! Store new sites in the next cycle dungeon.
    call p % preparePayload(payload)
    payload % weight = collDat % implicitWeight
    do i = 1, collDat % n
      call self % newPhysicalParticleState(collDat, fissionHandlePtr, payload, newState)
      call dungeon % detain(newState)

      ! Report birth of new particle
      call tally % reportSpawn(N_FISSION, p, newState)

    end do

  end subroutine implicit

  !!
  !! Extendable initialisation procedure
  !!
  subroutine init(self, dict)
    class(collisionProcessor), intent(inout) :: self
    class(dictionary), intent(in)            :: dict
    character(*), parameter                  :: HERE = 'init (collisionProcessor_inter.f90)'

    ! For now does nothing.

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(collisionProcessor), intent(inout) :: self

    self % currentMaterialPtr => null()
    self % xsData => null()

  end subroutine kill

  !!
  !!
  !!
  subroutine setCurrentMaterialPtr(self, materialPtr)
    class(collisionProcessor), intent(inout)    :: self
    class(neutronMaterial), pointer, intent(in) :: materialPtr
    character(*), parameter                     :: HERE = 'setMaterialPtr (collisionProcessor_inter.f90)'

    self % currentMaterialPtr => materialPtr
    if (.not. associated(self % currentMaterialPtr)) &
    call fatalError(HERE, 'Unable to retrieve material pointer.')

  end subroutine setCurrentMaterialPtr

  !!
  !!
  !!
  subroutine setNuclearDatabasePtr(self, xsDataPtr)
    class(collisionProcessor), intent(inout)    :: self
    class(nuclearDatabase), pointer, intent(in) :: xsDataPtr
    character(*), parameter                     :: HERE = 'setNuclearDatabasePtr'

    self % xsData => xsDataPtr
    if (.not. associated(self % xsData)) &
    call fatalError(HERE, 'Unable to retrieve active nuclear database pointer.')

  end subroutine setNuclearDatabasePtr

end module collisionProcessor_inter
