module collisionProcessor_inter

  use dictionary_class,       only : dictionary
  use endfConstants
  use errors_mod,             only : fatalError
  use genericProcedures,      only : numToChar
  use numPrecision
  use particleDungeon_class,  only : particleDungeon
  use physicalParticle_inter, only : castPhysicalParticlePtr, physicalParticle
  use RNG_class,              only : RNG
  use tallyCodes
  use tallyAdmin_class,       only : tallyAdmin
  use transportObject_inter,  only : transportObject

  implicit none
  private

  !!
  !! Data package with all relevant data about the collision to move beetween customisable
  !! procedures
  !!
  type, public :: collisionData
    integer(shortInt) :: matIdx  = -1   !! Material Index at collision
    integer(shortInt) :: nucIdx  = -1   !! Nuclide Index of target
    integer(shortInt) :: MT      = 0    !! MT Number of Realction
    real(defReal)     :: muL     = ONE  !! Cosine of deflection angle in LAB-frame
    real(defReal)     :: A       = ZERO !! Target Mass [Neutron Mass]
    real(defReal)     :: kT      = ZERO !! Target temperature [MeV]
    real(defReal)     :: E       = ZERO !! Collision energy (could be relative to target) [MeV]
  end type


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
  contains
    ! Master non-overridable procedures
    procedure, non_overridable :: collide

    ! Extendable procedures
    procedure :: init

    ! Customisable deffered procedures
    procedure(sampleCollision), deferred :: sampleCollision
    procedure(collisionAction), deferred :: implicit
    procedure(collisionAction), deferred :: elastic
    procedure(collisionAction), deferred :: inelastic
    procedure(collisionAction), deferred :: capture
    procedure(collisionAction), deferred :: fission
    procedure(collisionAction), deferred :: cutoffs

  end type collisionProcessor

  !! Extandable procedures
  public :: init


  abstract interface
    !!
    !! Procedure interface for all customisable actions associated with
    !! processing of sollision event (scatter, fission etc.)
    !!
    subroutine collisionAction(self, p, tally, collDat, thisCycle, nextCycle)
      import :: collisionProcessor, collisionData, tallyAdmin, particleDungeon, physicalParticle
      class(collisionProcessor), intent(inout) :: self
      class(physicalParticle), intent(inout)   :: p
      type(tallyAdmin), intent(inout)          :: tally
      type(collisionData), intent(inout)       :: collDat
      class(particleDungeon), intent(inout)    :: thisCycle, nextCycle
    end subroutine collisionAction

    !!
    !!
    !!
    subroutine sampleCollision(self, p, collDat)
      import :: collisionProcessor, collisionData, defReal, physicalParticle
      class(collisionProcessor), intent(inout) :: self
      class(physicalParticle), intent(inout)   :: p
      type(collisionData), intent(inout)       :: collDat
    end subroutine sampleCollision
    
  end interface

contains

  !!
  !! Generic flow of collision processing
  !!
  subroutine collide(self, object, tally, thisCycle, nextCycle)
    class(collisionProcessor), intent(inout) :: self
    class(transportObject), intent(inout)    :: object
    type(tallyAdmin), intent(inout)          :: tally
    class(particleDungeon), intent(inout)    :: thisCycle, nextCycle
    class(physicalParticle), pointer         :: p
    logical(defBool)                         :: virtual
    type(collisionData)                      :: collDat
    character(*), parameter                  :: here = 'collide (collisionProcessor.f90)'

    ! Downcast transportObject to physicalParticle.
    p => castPhysicalParticlePtr(object, .true.)

    ! Load material index into data package
    collDat % matIdx = p % getMaterialIdx()

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
    if (collDat % MT /= noInteraction) call self % implicit(p, tally, collDat, thisCycle, nextCycle)

    ! Select physics to be processed based on MT number
    select case(collDat % MT)
      case(N_N_elastic, macroAllScatter)
        call self % elastic(p, tally, collDat, thisCycle, nextCycle)

      case(N_N_inelastic, macroIEScatter)
        call self % inelastic(p, tally, collDat, thisCycle, nextCycle)

      case(N_DISAP, macroCapture)
        call self % capture(p, tally, collDat, thisCycle, nextCycle)

      case(N_FISSION, macroFission)
        call self % fission(p, tally, collDat, thisCycle, nextCycle)

      case(noInteraction)
        ! Do nothing

      case default
        call fatalError(here, 'Unsupported MT number: '//numToChar(collDat % MT)//'.')

    end select

    ! Apply post collision implicit treatments
    call self % cutoffs(p, tally, collDat, thisCycle, nextCycle)

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
  !! Extendable initialisation procedure
  !!
  subroutine init(self, dict)
    class(collisionProcessor), intent(inout) :: self
    class(dictionary), intent(in)            :: dict

    ! For now does nothing

  end subroutine init

end module collisionProcessor_inter
