module neutronMGCollisionProcessor_inter

  use collisionData_class,         only : collisionData
  use collisionProcessor_inter,    only : collisionProcessor, init_super => init
  use dictionary_class,            only : dictionary
  use endfConstants,               only : macroFission, macroIEscatter
  use errors_mod,                  only : fatalError
  use fissionMG_class,             only : fissionMG, fissionMG_TptrCast
  use genericProcedures,           only : rotateVector
  use MGCollisionData_class,       only : castMGCollisionDataPtr, MGCollisionData
  use MGNeutron_class,             only : castMGNeutronPtr, MGNeutron
  use mgNeutronMaterial_inter,     only : mgNeutronMaterial, mgNeutronMaterial_CptrCast
  use MGNeutronState_class,        only : MGNeutronState
  use MGParticleState_class,       only : buildMGParticleStatePayload, castBuildMGParticleStatePayloadPtr
  use multiScatterMG_class,        only : multiScatterMG, multiScatterMG_CptrCast
  use neutronXsPackages_class,     only : neutronMacroXSs
  use nuclearDataReg_mod,          only : getNeutronMG
  use numPrecision
  use physicalParticle_inter,      only : physicalParticle
  use physicalParticleState_class, only : physicalParticleState
  use reactionHandle_inter,        only : reactionHandle
  use RNG_class,                   only : RNG
  use transportObjectState_class,  only : buildTransportObjectStatePayload

  implicit none
  private

  ! Public procedures.
  public :: init

  !!
  !!
  !!
  type, public, abstract, extends(collisionProcessor) :: neutronMGCollisionProcessor
    private
  contains
    procedure :: allocateCollisionData
    procedure :: getCurrentMaterialMacroXSs
    procedure :: getFissionHandlePtr
    procedure :: inelastic
    procedure :: init
    procedure :: getImplicitCondition
    procedure :: newPhysicalParticleState
    procedure :: sampleCollision
  end type neutronMGCollisionProcessor

contains
  !!
  !!
  !!
  subroutine allocateCollisionData(self, collDat)
    class(neutronMGCollisionProcessor), intent(in) :: self
    class(collisionData), allocatable, intent(out) :: collDat

    allocate(MGCollisionData :: collDat)

  end subroutine allocateCollisionData

  !!
  !!
  !!
  subroutine getCurrentMaterialMacroXSs(self, energyGroup, rand, macroXSs)
    class(neutronMGCollisionProcessor), intent(in) :: self
    integer(shortInt), intent(in)                  :: energyGroup
    type(RNG), intent(inout)                       :: rand
    type(neutronMacroXSs), intent(out)             :: macroXSs
    class(mgNeutronMaterial), pointer              :: MGNeutronMaterialPtr
    character(*), parameter :: HERE = 'getCurrentMaterialMacroXSs (neutronMGCollisionProcessor_inter.f90)'

    ! Get pointer to MG neutron material.
    MGNeutronMaterialPtr => mgNeutronMaterial_CptrCast(self % getCurrentMaterialPtr())
    if (.not. associated(MGNeutronMaterialPtr)) &
    call fatalError(HERE, 'Failed to retrieve MG neutron material.')
    call MGNeutronMaterialPtr % getMacroXSs(energyGroup, macroXSs, rand)

  end subroutine getCurrentMaterialMacroXSs

  !!
  !!
  !!
  function getFissionHandlePtr(self, collDat) result(fissionHandlePtr)
    class(neutronMGCollisionProcessor), intent(in) :: self
    class(collisionData), intent(in)               :: collDat
    class(reactionHandle), pointer                 :: fissionHandlePtr

    fissionHandlePtr => fissionMG_TptrCast(self % getReaction(macroFission, collDat % matIdx))

  end function getFissionHandlePtr

  !!
  !!
  !!
  subroutine inelastic(self, collDat, p)
    class(neutronMGCollisionProcessor), intent(in) :: self
    class(collisionData), intent(inout)            :: collDat
    class(physicalParticle), intent(inout)         :: p
    class(multiScatterMG), pointer                 :: scatter
    real(defReal)                                  :: phi
    type(MGCollisionData), pointer                 :: MGCollisionDataPtr
    type(MGNeutron), pointer                       :: MGNeutronPtr
    character(*), parameter                        :: HERE = 'inelastic (neutronMGCollisionProcessor_inter.f90)'

    ! First downcast physical particle to MG neutron.
    MGCollisionDataPtr => castMGCollisionDataPtr(collDat, .true.)
    MGNeutronPtr => castMGNeutronPtr(p, .true.)

    ! Assign MT number and get scattering object.
    MGCollisionDataPtr % MT = macroIEscatter
    scatter => multiScatterMG_CptrCast(self % getReaction(MGCollisionDataPtr % MT, MGCollisionDataPtr % matIdx))
    if (.not. associated(scatter)) call fatalError(HERE, "Failed to get scattering reaction object for MG neutron")

    ! Sample Mu and G_out
    call scatter % sampleOut(collDat % muL, phi, MGCollisionDataPtr % G_out, MGCollisionDataPtr % G_in, &
                             MGCollisionDataPtr % RNGPtr)

    ! Update neutron state
    call MGNeutronPtr % setEnergyGroup(MGCollisionDataPtr % G_out)
    MGCollisionDataPtr % weight = MGCollisionDataPtr % weight * scatter % production(MGCollisionDataPtr % G_in, &
                                                                                     MGCollisionDataPtr % G_out)
    call MGNeutronPtr % setWeight(MGCollisionDataPtr % weight)
    call MGNeutronPtr % rotate(collDat % muL, phi)

  end subroutine inelastic

  !!
  !!
  !!
  subroutine init(self, dict)
    class(neutronMGCollisionProcessor), intent(inout) :: self
    class(dictionary), intent(in)                     :: dict

    ! Initialise superclass.
    call init_super(self, dict)

    ! Load and verify nuclear data pointer.
    call self % setNuclearDatabasePtr(getNeutronMG())

  end subroutine init

  !!
  !!
  !!
  elemental function getImplicitCondition(self) result(isIt)
    class(neutronMGCollisionProcessor), intent(in) :: self
    logical(defBool)                               :: isIt

    isIt = self % getCurrentMaterialIsFissile()

  end function getImplicitCondition

  !!
  !!
  !!
  subroutine newPhysicalParticleState(self, collDat, fissionHandlePtr, payload, newState)
    class(neutronMGCollisionProcessor), intent(in)         :: self
    class(collisionData), intent(in)                       :: collDat
    class(reactionHandle), pointer, intent(in)             :: fissionHandlePtr
    class(buildTransportObjectStatePayload), intent(inout) :: payload
    class(physicalParticleState), allocatable, intent(out) :: newState
    type(buildMGParticleStatePayload), pointer             :: payloadPtr
    type(MGCollisionData), pointer                         :: MGCollisionDataPtr
    type(fissionMG), pointer                               :: fissionMGPtr
    real(defReal)                                          :: mu, phi

    ! Downcast classes to correct types.
    payloadPtr => castBuildMGParticleStatePayloadPtr(payload)
    MGCollisionDataPtr => castMGCollisionDataPtr(collDat)
    fissionMGPtr => fissionMG_TptrCast(fissionHandlePtr)

    ! Allocate new state.
    allocate(MGNeutronState :: newState)

    ! Finalise payload then initialise state.
    call fissionMGPtr % sampleOut(mu, phi, payloadPtr % energyGroup, MGCollisionDataPtr % G_in, MGCollisionDataPtr % RNGPtr)
    payloadPtr % uGlobal = rotateVector(MGCollisionDataPtr % u, mu, phi)
    call newState % init(payload)

  end subroutine newPhysicalParticleState

  !!
  !!
  !!
  subroutine sampleCollision(self, p, collDat)
    class(neutronMGCollisionProcessor), intent(inout) :: self
    class(physicalParticle), intent(in)               :: p
    class(collisionData), intent(inout)               :: collDat
    class(mgNeutronMaterial), pointer                 :: MGNeutronMaterialPtr
    real(defReal)                                     :: randomNumber
    type(MGCollisionData), pointer                    :: MGCollisionDataPtr
    type(neutronMacroXSs)                             :: macroXSs

    ! Downcast physical particle to an MG neutron.
    MGCollisionDataPtr => castMGCollisionDataPtr(collDat)

    ! Get and verify material pointer.
    MGNeutronMaterialPtr => mgNeutronMaterial_CptrCast(self % getMaterialPtr(MGCollisionDataPtr % matIdx))
    call self % setCurrentMaterialPtr(MGNeutronMaterialPtr)

    ! Select Main reaction channel
    call MGNeutronMaterialPtr % getMacroXSs(MGCollisionDataPtr % G_in, macroXSs, MGCollisionDataPtr % RNGPtr)
    call MGCollisionDataPtr % RNGPtr % generate(randomNumber)
    MGCollisionDataPtr % MT = macroXSs % invert(randomNumber)
    MGCollisionDataPtr % sigma_fission = macroXSs % fission
    MGCollisionDataPtr % sigma_nuFiss = macroXSs % nuFission
    MGCollisionDataPtr % sigma_tot = macroXSs % total

  end subroutine sampleCollision

end module neutronMGCollisionProcessor_inter