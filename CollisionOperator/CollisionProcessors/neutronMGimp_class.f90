module neutronMGimp_class

  use collisionProcessor_inter,    only : collisionProcessor, collisionData ,init_super => init
  use dictionary_class,            only : dictionary
  use endfConstants
  use errors_mod,                  only : fatalError
  use fissionMG_class,             only : fissionMG, fissionMG_TptrCast
  use genericProcedures,           only : numToChar, rotateVector
  use geometryReg_mod,             only : gr_fieldIdx => fieldIdx, gr_fieldPtr => fieldPtr
  use MGNeutron_class,             only : castMGNeutronPtr, MGNeutron
  use mgNeutronDatabase_inter,     only : mgNeutronDatabase
  use mgNeutronMaterial_inter,     only : mgNeutronMaterial, mgNeutronMaterial_CptrCast
  use MGParticleState_class,       only : castMGParticleStatePtr, MGParticleState
  use multiScatterMG_class,        only : multiScatterMG, multiScatterMG_CptrCast
  use neutronXsPackages_class,     only : neutronMacroXSs
  use nuclearDatabase_inter,       only : nuclearDatabase
  use nuclearDataReg_mod,          only : ndReg_getNeutronMG => getNeutronMG
  use numPrecision
  use particleDungeon_class,       only : particleDungeon
  use physicalParticle_inter,      only : physicalParticle
  use physicalParticleState_class, only : castPhysicalParticleStatePtr, physicalParticleState
  use reactionHandle_inter,        only : reactionHandle
  use RNG_class,                   only : RNG
  use tallyAdmin_class,            only : tallyAdmin
  use universalVariables,          only : nameWW
  use weightWindowsField_class,    only : weightWindowsField, weightWindowsField_TptrCast

  implicit none
  private

  !!
  !! Scalar collision processor for MG neutrons
  !!   -> Preforms implicit fission site generation
  !!   -> Preforms analog capture
  !!   -> Treats fission as capture (only implicit generation of 2nd-ary neutrons)
  !!   -> Does not create secondary non-neutron projectiles
  !!   -> Supports the use of weight windows
  !!
  !! Settings:
  !!  weightWindows -> uses a weight windows field (off by default)
  !!  maxSplit -> maximum number of splits allowed per particle (default = 1000)
  !!
  !! Sample dictionary input:
  !!   collProcName {
  !!   type            neutronMGimp;
  !!   #weightWindows  <logical>;#
  !!   #maxSplit       <integer>;#
  !!   }
  !!
  type, public, extends(collisionProcessor) :: neutronMGimp
    private
    class(mgNeutronDatabase), pointer, public :: xsData => null()
    class(mgNeutronMaterial), pointer, public :: mat    => null()

    !! Settings - private
    integer(shortInt) :: maxSplit

    ! Variance reduction options
    logical(defBool)  :: weightWindows
    type(weightWindowsField), pointer :: weightWindowsMap

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

    ! Variance reduction procedures
    procedure, private :: split
    procedure, private :: russianRoulette

  end type neutronMGimp

contains

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self, dict)
    class(neutronMGimp), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    integer(shortInt)                  :: idx
    character(*), parameter :: Here = 'init (neutronMGimp_class.f90)'

    ! Call superclass
    call init_super(self, dict)

    ! Obtain settings for variance reduction
    call dict % getOrDefault(self % maxSplit,'maxSplit', 1000)
    call dict % getOrDefault(self % weightWindows,'weightWindows', .false.)

    ! Sets up the weight windows field
    if (self % weightWindows) then
      idx = gr_fieldIdx(nameWW)
      self % weightWindowsMap => weightWindowsField_TptrCast(gr_fieldPtr(idx))
    end if

  end subroutine init

  !!
  !! Samples collision without any implicit treatment
  !!
  subroutine sampleCollision(self, p, collDat)
    class(neutronMGimp), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    type(collisionData), intent(inout)     :: collDat
    real(defReal)                          :: randomNumber
    type(MGNeutron), pointer               :: MGNeutronPtr
    type(neutronMacroXSs)                  :: macroXSs
    type(RNG), pointer                     :: RNGPtr
    character(*), parameter                :: Here =' sampleCollision (neutronMGimp_class.f90)'

    ! Verify that particle is MG neutron
    MGNeutronPtr => castMGNeutronPtr(p, .true.)

    ! Verify and load nuclear data pointer
    self % xsData => ndReg_getNeutronMG()
    if (.not. associated(self % xsData)) call fatalError(Here, "Failed to get active database for MG Neutron")

    ! Get and verify material pointer
    self % mat => mgNeutronMaterial_CptrCast(self % xsData % getMaterial(MGNeutronPtr % getMaterialIdx()))
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
    class(neutronMGimp), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon), intent(inout)  :: thisCycle, nextCycle
    type(fissionMG), pointer               :: fission
    type(MGNeutron), pointer               :: MGNeutronPtr
    type(MGParticleState), pointer         :: MGParticleStatePtr, preHistoryStatePtr
    type(neutronMacroXSs)                  :: macroXSs
    type(RNG), pointer                     :: RNGPtr
    integer(shortInt)                      :: G, G_out, n, i
    real(defReal)                          :: k_eff, mu, phi, randomNumber, sig_nufiss, sig_tot, wgt, w0
    real(defReal), dimension(3)            :: rGlobal, u, uGlobal
    character(*), parameter                :: Here = 'implicit (neutronMGimp_class.f90)'

    MGNeutronPtr => castMGNeutronPtr(p, .true.)
    if (self % mat % isFissile()) then
      ! Obtain required data
      wgt = MGNeutronPtr % getWeight()               ! Current weight
      preHistoryStatePtr => castMGParticleStatePtr(MGNeutronPtr % getPreHistoryStatePtr(), .true.)
      w0 = preHistoryStatePtr % getWeight() ! Starting weight
      k_eff = MGNeutronPtr % getKEff()            ! k_eff for normalisation
      RNGPtr => MGNeutronPtr % getRNGPtr()
      call RNGPtr % generate(randomNumber)     ! Random number to sample sites

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
      if (.not. associated(fission)) call fatalError(Here, 'Failed to retrieve fissionMG.')

      ! Store new sites in the next cycle dungeon
      wgt = sign(w0, wgt)
      rGlobal = MGNeutronPtr % getGlobalPosition()
      uGlobal = MGNeutronPtr % getGlobalDirection()
      do i = 1, n
        call fission % sampleOut(mu, phi, G_out, G, RNGPtr)
        u = rotateVector(uGlobal, mu, phi)

        ! Copy extra detail from parent particle (i.e. time, flags ect.)
        MGParticleStatePtr => castMGParticleStatePtr(MGNeutronPtr % updateAndGetCurrentStatePtr(), .true.)

        ! Overwrite position, direction, energy group and weight
        call MGParticleStatePtr % setGlobalPosition(rGlobal)
        call MGParticleStatePtr % setGlobalDirection(u)
        call MGParticleStatePtr % setEnergyGroup(G_out)
        call MGParticleStatePtr % setWeight(wgt)

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
    class(neutronMGimp), intent(inout)     :: self
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
    class(neutronMGimp), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon), intent(inout)  :: thisCycle
    class(particleDungeon), intent(inout)  :: nextCycle
    class(multiScatterMG), pointer         :: scatter
    integer(shortInt)                      :: G, G_out   ! Post-collision energy group
    real(defReal)                          :: phi, w_mul
    type(MGNeutron), pointer               :: MGNeutronPtr
    character(*), parameter                :: Here = "inelastic (neutronMGimp_class.f90)"

    MGNeutronPtr => castMGNeutronPtr(p, .true.)

    ! Assign MT number
    collDat % MT = macroIEscatter

    ! Get Scatter object
    scatter => multiScatterMG_CptrCast(self % xsData % getReaction(macroIEscatter, collDat % matIdx))
    if (.not. associated(scatter)) call fatalError(Here, "Failed to get scattering reaction object for MG neutron")

    ! Sample Mu and G_out
    G = MGNeutronPtr % getEnergyGroup()
    call scatter % sampleOut(collDat % muL, phi, G_out, G, MGNeutronPtr % getRNGPtr())

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
    class(neutronMGimp), intent(inout)     :: self
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
    class(neutronMGimp), intent(inout)     :: self
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
    class(neutronMGimp), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon), intent(inout)  :: thisCycle, nextCycle
    real(defReal)                          :: avWgt, maxWgt, minWgt, weight
    real(defReal), dimension(3)            :: val

    weight = p % getWeight()

    if (p % getIsDead()) then
      ! Do nothing !

    ! Weight Windows treatment
    elseif (self % weightWindows) then
      val = self % weightWindowsMap % at(p)
      minWgt = val(1)
      maxWgt = val(2)
      avWgt  = val(3)

      ! If a particle is outside the WW map and all the weight limits
      ! are zero nothing happens. NOTE: this holds for positive weights only
      if (maxWgt < weight .and. maxWgt /= ZERO .and. p % getSplitsNumber() < self % maxSplit) then
        call self % split(p, tally, thisCycle, maxWgt)

      elseif (weight < minWgt) then
        call self % russianRoulette(p, avWgt)

      end if

    end if

  end subroutine cutoffs

  !!
  !! Perform Russian roulette on a particle
  !!
  subroutine russianRoulette(self, p, avWgt)
    class(neutronMGimp), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    real(defReal), intent(in)              :: avWgt
    real(defReal)                          :: randomNumber
    type(RNG), pointer                     :: RNGPtr

    RNGPtr => p % getRNGPtr()
    call RNGPtr % generate(randomNumber)
    if (randomNumber < ONE - p % getWeight() / avWgt) then
      call p % setIsDead(.true.)

    else
      call p % setWeight(avWgt)

    end if

  end subroutine russianRoulette

  !!
  !! Split particle which has too large a weight
  !!
  subroutine split(self, p, tally, thisCycle, maxWgt)
    class(neutronMGimp), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    type(tallyAdmin), intent(inout)        :: tally
    class(particleDungeon), intent(inout)  :: thisCycle
    real(defReal), intent(in)              :: maxWgt
    class(physicalParticleState), pointer  :: physicalParticleStatePtr
    integer(shortInt)                      :: i, mult, nSplits
    real(defReal)                          :: newWeight, weight

    ! This value must be at least 2
    weight = p % getWeight()
    mult = ceiling(weight / maxWgt)

    ! Limit maximum split
    nSplits = p % getSplitsNumber()
    if (self % maxSplit < mult + nSplits) mult = self % maxSplit - nSplits + 1

    ! Copy particle to a particle state
    ! Note that particleState doesn't have property splitCount, so it is reset
    ! to 0 for the new particle
    newWeight = weight / mult
    physicalParticleStatePtr => castPhysicalParticleStatePtr(p % updateAndGetCurrentStatePtr(), .true.)
    call physicalParticleStatePtr % setWeight(newWeight)

    ! Add split particle's to the dungeon
    do i = 1, mult - 1
      call thisCycle % detain(physicalParticleStatePtr)
      call tally % reportSpawn(N_N_SPLIT, p, physicalParticleStatePtr)

    end do

    ! Update particle split coun and decrease original particle weight.
    call p % setSplitsNumber(nSplits + mult)
    call p % setWeight(newWeight)

  end subroutine split

end module neutronMGimp_class
