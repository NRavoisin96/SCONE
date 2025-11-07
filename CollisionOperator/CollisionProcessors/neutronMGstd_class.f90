module neutronMGstd_class

  use collisionProcessor_inter, only : collisionProcessor, collisionData ,init_super => init
  use dictionary_class,         only : dictionary
  use endfConstants
  use errors_mod,               only : fatalError
  use fissionMG_class,          only : fissionMG, fissionMG_TptrCast
  use genericProcedures,        only : numToChar, rotateVector
  use MGNeutron_class,          only : castMGNeutronPtr, MGNeutron
  use mgNeutronDatabase_inter,  only : mgNeutronDatabase
  use mgNeutronMaterial_inter,  only : mgNeutronMaterial, mgNeutronMaterial_CptrCast
  use MGParticleState_class,    only : castMGParticleStatePtr, MGParticleState
  use multiScatterMG_class,     only : multiScatterMG, multiScatterMG_CptrCast
  use neutronXsPackages_class,  only : neutronMacroXSs
  use nuclearDatabase_inter,    only : nuclearDatabase
  use nuclearDataReg_mod,       only : ndReg_getNeutronMG => getNeutronMG
  use numPrecision
  use particleDungeon_class,    only : particleDungeon
  use physicalParticle_inter,   only : physicalParticle
  use RNG_class,                only : RNG
  use tallyAdmin_class,         only : tallyAdmin

  implicit none
  private

  !!
  !! Standard (default) scalar collision processor for MG neutrons
  !!   -> Preforms implicit fission site generation
  !!   -> Preforms analog capture
  !!   -> Treats fission as capture (only implicit generation of 2nd-ary neutrons)
  !!   -> Does not create secondary non-neutron projectiles
  !!
  !! Settings:
  !!  NONE
  !!
  !! Sample dictionary input:
  !!   collProcName {
  !!   type            neutronMGstd;
  !!   }
  !!
  type, public, extends(collisionProcessor) :: neutronMGstd
    private
    class(mgNeutronDatabase), pointer, public :: xsData => null()
    class(mgNeutronMaterial), pointer, public :: mat    => null()
  contains
    ! Initialisation procedure
    procedure :: init

    ! Implementation of customisable procedures
    procedure :: sampleCollision
    procedure :: implicit
    procedure :: elastic
    procedure :: inelastic
    procedure :: capture
    procedure :: fission
    procedure :: cutoffs
  end type neutronMGstd

contains

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self, dict)
    class(neutronMGstd), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    character(*), parameter :: Here = 'init (neutronMGstd_class.f90)'

    ! Call superclass
    call init_super(self, dict)

  end subroutine init

  !!
  !! Samples collision without any implicit treatment
  !!
  subroutine sampleCollision(self, p, collDat)
    class(neutronMGstd), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    type(collisionData), intent(inout)     :: collDat
    real(defReal)                          :: randomNumber
    type(MGNeutron), pointer               :: MGNeutronPtr
    type(neutronMacroXSs)                  :: macroXSs
    type(RNG), pointer                     :: RNGPtr
    character(*), parameter                :: Here = 'sampleCollision (neutronMGstd_class.f90)'

    ! Verify that particle is MG neutron
    MGNeutronPtr => castMGNeutronPtr(p, .true.)

    ! Verify and load nuclear data pointer
    self % xsData => ndReg_getNeutronMG()
    if (.not. associated(self % xsData)) call fatalError(Here, "Failed to get active database for MG Neutron")

    ! Get and verify material pointer
    self % mat => mgNeutronMaterial_CptrCast( self % xsData % getMaterial(MGNeutronPtr % getMaterialIdx()))
    if (.not. associated(self % mat)) call fatalError(Here, "Failed to get MG Neutron Material")

    ! Select Main reaction channel
    RNGPtr => MGNeutronPtr % getRNGPtr()
    call self % mat % getMacroXSs(MGNeutronPtr % getEnergyGroup(), macroXSs, RNGPtr)
    call RNGPtr % generate(randomNumber)
    collDat % MT = macroXSs % invert(randomNumber)

  end subroutine sampleCollision

  !!
  !! Preform implicit treatment
  !!
  subroutine implicit(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronMGstd), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon), intent(inout)  :: thisCycle, nextCycle
    integer(shortInt)                      :: G, G_out, i, n
    real(defReal)                          :: k_eff, mu, phi, randomNumber, sig_nufiss, sig_tot, wgt, w0
    real(defReal), dimension(3)            :: rGlobal, u, uGlobal
    type(fissionMG), pointer               :: fission
    type(MGNeutron), pointer               :: MGNeutronPtr
    type(MGParticleState), pointer         :: MGParticleStatePtr, preHistoryStatePtr
    type(neutronMacroXSs)                  :: macroXSs
    type(RNG), pointer                     :: RNGPtr
    character(*), parameter                :: Here = 'implicit (neutronMGstd_class.f90)'

    MGNeutronPtr => castMGNeutronPtr(p, .true.)
    if (self % mat % isFissile()) then
      ! Obtain required data
      wgt = MGNeutronPtr % getWeight()                ! Current weight
      preHistoryStatePtr => castMGParticleStatePtr(MGNeutronPtr % getPreHistoryStatePtr(), .true.)
      w0 = preHistoryStatePtr % getWeight() ! Starting weight
      k_eff = MGNeutronPtr % getKEff() ! k_eff for normalisation
      RNGPtr => MGNeutronPtr % getRNGPtr()
      call RNGptr % generate(randomNumber)    ! Random number to sample sites

      G = MGNeutronPtr % getEnergyGroup()
      call self % mat % getMacroXSs(G, macroXSs, RNGPtr)
      sig_tot = macroXSs % total
      sig_nuFiss = macroXSs % nuFission

      ! Sample number of fission sites generated
      !n = int(wgt * sig_nuFiss/(sig_tot*k_eff) + r1, shortInt)
      n = int(abs((wgt * sig_nuFiss) / (w0 * sig_tot * k_eff)) + randomNumber, shortInt)

      ! Shortcut if no particles were samples
      if (n < 1) return

      ! Get Fission reaction object
      fission => fissionMG_TptrCast(self % xsData % getReaction(macroFission, collDat % matIdx))
      if (.not. associated(fission)) call fatalError(Here, 'Failed to getrive fissionMG reaction object')

      ! Store new sites in the next cycle dungeon
      wgt = sign(w0, wgt)
      rGlobal = MGNeutronPtr % getGlobalPosition()
      uGlobal = MGNeutronPtr % getGlobalDirection()
      do i= 1, n
        call fission % sampleOut(mu, phi, G_out, G, RNGPtr)
        u = rotateVector(uGlobal, mu, phi)

        ! Copy extra detail from parent particle (i.e. time, flags ect.)
        MGParticleStatePtr => castMGParticleStatePtr(MGNeutronPtr % updateAndGetCurrentStatePtr(), .true.)

        ! Overwrite position, direction, energy group and weight
        call MGParticleStatePtr % setGlobalPosition(rGlobal)
        call MGParticleStatePtr % setGlobalDirection(u)
        call MGParticleStatePtr % setEnergyGroup(G_out)
        call MGParticleStatePtr % setWeight(wgt)
        call MGParticleStatePtr % setCollisionsNumber(0)

        call nextCycle % detain(MGParticleStatePtr)

        ! Report birth of new particle
        call tally % reportSpawn(N_FISSION, MGNeutronPtr, MGParticleStatePtr)

      end do

    end if

  end subroutine implicit

  !!
  !! Elastic Scattering
  !!
  subroutine elastic(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronMGstd), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon), intent(inout)  :: thisCycle, nextCycle

    ! Do nothing. Should not be called

  end subroutine elastic

  !!
  !! Preform scattering
  !!
  subroutine inelastic(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronMGstd), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon), intent(inout)  :: thisCycle, nextCycle
    class(multiScatterMG), pointer         :: scatter
    integer(shortInt)                      :: G, G_out   ! Post-collision energy group
    real(defReal)                          :: phi, w_mul
    type(MGNeutron), pointer               :: MGNeutronPtr
    character(*), parameter                :: Here = "inelastic (neutronMGstd_class.f90)"

    MGNeutronPtr => castMGNeutronPtr(p, .true.)

    ! Assign MT number
    collDat % MT = macroIEscatter

    ! Get Scatter object
    scatter => multiScatterMG_CptrCast(self % xsData % getReaction(macroIEscatter, collDat % matIdx))
    if (.not. associated(scatter)) call fatalError(Here, "Failed to get scattering reaction object for MG neutron")

    ! Sample Mu and G_out
    G = MGNeutronPtr % getEnergyGroup()
    call scatter % sampleOut(collDat % muL, phi, G_out, G, p % getRNGPtr())

    ! Read scattering multiplicity
    w_mul = scatter % production(G, G_out)

    ! Update neutron state
    call MGNeutronPtr % setEnergyGroup(G_out)
    call MGNeutronPtr % setWeight(MGNeutronPtr % getWeight() * w_mul)
    call p % rotate(collDat % muL, phi)

  end subroutine inelastic

  !!
  !! Preform capture
  !!
  subroutine capture(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronMGstd), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon), intent(inout)  :: thisCycle, nextCycle

    call p % setIsDead(.true.)

  end subroutine capture

  !!
  !! Preform fission
  !!
  subroutine fission(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronMGstd), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon), intent(inout)  :: thisCycle, nextCycle

    call p % setIsDead(.true.)

  end subroutine fission

  !!
  !! Applay cutoffs or post-collision implicit treatment
  !!
  subroutine cutoffs(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronMGstd), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon), intent(inout)  :: thisCycle, nextCycle

    ! Do nothing

  end subroutine cutoffs

end module neutronMGstd_class
