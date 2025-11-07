module neutronCEstd_class

  use CENeutron_class,                   only : castCENeutronPtr, CENeutron
  use CENeutronState_class,              only : castCENeutronStatePtr, CENeutronState, newCENeutronState
  use CEParticleState_class,             only : buildCEParticleStatePayload
  use collisionProcessor_inter,          only : collisionData
  use dictionary_class,                  only : dictionary
  use endfConstants
  use errors_mod,                        only : fatalError
  use fissionCE_class,                   only : fissionCE, fissionCE_TptrCast
  use genericProcedures,                 only : numToChar, rotateVector
  use neutronCECollisionProcessor_inter, only : neutronCECollisionProcessor, init_super => init
  use neutronXsPackages_class,           only : neutronMicroXSs
  use numPrecision
  use particleDungeon_class,             only : particleDungeon
  use physicalParticle_inter,            only : physicalParticle
  use RNG_class,                         only : RNG
  use tallyAdmin_class,                  only : tallyAdmin

  implicit none
  private

  !!
  !! Standard (default) scalar collision processor for CE neutrons
  !!   -> Preforms implicit fission site generation
  !!   -> Preforms analog capture
  !!   -> Treats fission as capture (only implicit generation of 2nd-ary neutrons)
  !!   -> Does not create secondary non-neutron projectiles
  !!
  !! Settings:
  !!  minE    -> minimum energy cut-off [MeV] (default = 1.0E-11)
  !!  maxE    -> maximum energy. Higher energies are set to maximum (not re-rolled) [MeV]
  !!             (default = 20.0)
  !!  threshE -> Energy threshold for explicit treatment of target nuclide movement [-].
  !!             Target movement is sampled if neutron energy E < kT * threshE where
  !!             kT is target material temperature in [MeV]. (default = 400.0)
  !!  threshA -> Mass threshold for explicit treatment of target nuclide movement [Mn].
  !!             Target movement is sampled if target mass A < threshA. (default = 1.0)
  !!  DBRCeMin -> Minimum energy to which DBRC is applied
  !!  DBRCeMax -> Maximum energy to which DBRC is applied
  !!
  !! Sample dictionary input:
  !!   collProcName {
  !!   type             neutronCEstd;
  !!   #minEnergy       <real>;#
  !!   #maxEnergy       <real>;#
  !!   #energyThreshold <real>;#
  !!   #massThreshold   <real>;#
  !!   }
  !!
  type, public, extends(neutronCECollisionProcessor) :: neutronCEstd
    private
  contains
    ! Initialisation procedure
    procedure :: init

    ! Implementation of customisable procedures
    procedure :: implicit
    procedure :: fission
    procedure :: cutoffs
  end type neutronCEstd

contains

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self, dict)
    class(neutronCEstd), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    character(*), parameter :: Here = 'init (neutronCEstd_class.f90)'

    ! Initialise superclass.
    call init_super(self, dict)

  end subroutine init

  !!
  !! Perform implicit treatment
  !!
  subroutine implicit(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon), intent(inout)  :: thisCycle, nextCycle
    type(buildCEParticleStatePayload)      :: payload
    type(CENeutron), pointer               :: CENeutronPtr
    type(CENeutronState)                   :: newState
    type(CENeutronState), pointer          :: preHistoryStatePtr
    type(fissionCE), pointer               :: fission
    type(neutronMicroXSs)                  :: microXSs
    type(RNG), pointer                     :: RNGPtr
    real(defReal)                          :: E, E_max, E_out, k_eff, mu, phi, randomNumber, &
                                              sig_nufiss, sig_tot, wgt, w0
    real(defReal), dimension(3)            :: uGlobal
    integer(shortInt)                      :: i, n
    character(*), parameter                :: Here = 'implicit (neutronCEstd_class.f90)'

    CENeutronPtr => castCENeutronPtr(p, .true.)

    ! Generate fission sites if nuclide is fissile
    if (self % getNuclideIsFissile()) then
      ! Obtain required data
      wgt = CENeutronPtr % getWeight()                ! Current weight
      preHistoryStatePtr => castCENeutronStatePtr(CENeutronPtr % getPreHistoryStatePtr(), .true.)
      w0 = preHistoryStatePtr % getWeight() ! Starting weight
      k_eff = CENeutronPtr % getKEff()            ! k_eff for normalisation
      RNGPtr => CENeutronPtr % getRNGPtr()
      call RNGPtr % generate(randomNumber)     ! Random number to sample sites

      ! Retrieve cross section at the energy used for reaction sampling
      call self % getNuclideMicroXS(collDat % E, collDat % kT, RNGPtr, microXSs)
      sig_nufiss = microXSs % nuFission
      sig_tot = microXSs % total

      ! Sample number of fission sites generated
      ! Support -ve weight particles
      n = int(abs((wgt * sig_nufiss) / (w0 * sig_tot * k_eff)) + randomNumber, shortInt)

      ! Shortcut particle generation if no particles were sampled
      if (n < 1) return

      ! Get fission Reaction
      fission => fissionCE_TptrCast(self % getReaction(N_FISSION, collDat % nucIdx))
      if (.not. associated(fission)) call fatalError(Here, 'Failed to retrieve fissionCE')

      ! Store new sites in the next cycle dungeon.
      call CENeutronPtr % preparePayload(payload)
      uGlobal = payload % uGlobal
      payload % weight = sign(w0, wgt)
      E = CENeutronPtr % getEnergy()
      E_max = self % getMaximumEnergy()
      do i = 1, n
        call fission % sampleOut(mu, phi, E_out, E, RNGPtr)
        payload % uGlobal = rotateVector(uGlobal, mu, phi)
        payload % energy = min(E_out, E_max)
        newState = newCENeutronState(payload)
        call nextCycle % detain(newState)

        ! Report birth of new particle
        call tally % reportSpawn(N_FISSION, CENeutronPtr, newState)

      end do
      
    end if

  end subroutine implicit

  !!
  !! Process fission reaction
  !!
  subroutine fission(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon), intent(inout)  :: thisCycle, nextCycle

    call p % setIsDead(.true.)

  end subroutine fission

  !!
  !! Apply cutoffs
  !!
  subroutine cutoffs(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon), intent(inout)  :: thisCycle, nextCycle
    type(CENeutron), pointer               :: CENeutronPtr

    CENeutronPtr => castCENeutronPtr(p, .true.)
    if (CENeutronPtr % getEnergy() < self % getMinimumEnergy()) call CENeutronPtr % setIsDead(.true.)

  end subroutine cutoffs

end module neutronCEstd_class