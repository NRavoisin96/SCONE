module neutronCEimp_class

  use CENeutron_class,                   only : castCENeutronPtr, CENeutron
  use CEParticleState_class,             only : castCEParticleStatePtr, CEParticleState
  use collisionProcessor_inter,          only : collisionData
  use dictionary_class,                  only : dictionary
  use endfConstants
  use errors_mod,                        only : fatalError
  use fissionCE_class,                   only : fissionCE, fissionCE_TptrCast
  use genericProcedures,                 only : numToChar, rotateVector
  use geometryReg_mod,                   only : gr_fieldIdx => fieldIdx, gr_fieldPtr => fieldPtr
  use neutronCECollisionProcessor_inter, only : init_super => init, neutronCECollisionProcessor
  use neutronXsPackages_class,           only : neutronMicroXSs
  use numPrecision
  use particleDungeon_class,             only : particleDungeon
  use physicalParticle_inter,            only : physicalParticle
  use RNG_class,                         only : RNG
  use tallyAdmin_class,                  only : tallyAdmin
  use uniFissSitesField_class,           only : uniFissSitesField, uniFissSitesField_TptrCast
  use universalVariables,                only : nameUFS, nameWW
  use weightWindowsField_class,          only : weightWindowsField, weightWindowsField_TptrCast

  implicit none
  private

  !!
  !! Standard scalar collision processor for CE neutrons
  !!   -> Preforms implicit or analog fission site generation
  !!   -> Preforms implicit or analog capture
  !!   -> Treats fission as capture (only implicit generation of 2nd-ary neutrons)
  !!   -> Does not create secondary non-neutron projectiles
  !!
  !! Settings:
  !!  minE    -> minimum energy cut-off [MeV] (default = 1.0E-11)
  !!  maxE    -> maximum energy. Higher energies are set to maximum (not re-rolled) [MeV]
  !!             (default = 20.0)
  !!  minWgt  -> minimum particle weight for rouletting (optional)
  !!  maxWgt  -> maximum particle weight for splitting (optional)
  !!  avgWgt  -> weight of a particle on surviving rouletting (optional)
  !!  impAbs  -> is implicit capture performed? (off by default)
  !!  impGen  -> are fission sites generated implicitly? (on by default)
  !!  UFS     -> uniform fission sites variance reduction
  !!  maxSplit -> maximum number of splits allowed per particle (default = 1000)
  !!  threshE  -> Energy threshold for explicit treatment of target nuclide movement [-].
  !!              Target movement is sampled if neutron energy E < kT * threshE where
  !!              kT is target material temperature in [MeV]. (default = 400.0)
  !!  threshA  -> Mass threshold for explicit treatment of target nuclide movement [Mn].
  !!              Target movment is sampled if target mass A < threshA. (default = 1.0)
  !!  DBRCeMin -> Minimum energy to which DBRC is applied
  !!  DBRCeMax -> Maximum energy to which DBRC is applied
  !!  splitting -> splits particles above certain weight (on by default)
  !!  roulette  -> roulettes particles below certain weight (off by defautl)
  !!  weightWindows -> uses a weight windows field (off by default)
  !!
  !! Sample dictionary input:
  !!   collProcName {
  !!   type            neutronCEimp;
  !!   #minEnergy      <real>;#
  !!   #maxEnergy      <real>;#
  !!   #energyTreshold <real>;#
  !!   #massTreshold   <real>;#
  !!   #splitting      <logical>;#
  !!   #roulette       <logical>;#
  !!   #minWgt         <real>;#
  !!   #maxWgt         <real>;#
  !!   #avgWgt         <real>;#
  !!   #impAbs         <logical>;#
  !!   #impGen         <logical>;#
  !!   #UFS            <logical>;#
  !!   #weightWindows  <logical>;#
  !!   #maxSplit       <integer>;#
  !!   }
  !!
  type, public, extends(neutronCECollisionProcessor) :: neutronCEimp
    private
    !! Nuclear Data block pointer -> public so it can be used by subclasses (protected member)
    class(uniFissSitesField), pointer :: ufsField => null()

    !! Settings - private
    real(defReal)     :: minWgt = ZERO, maxWgt = ZERO, avWgt = ZERO
    integer(shortInt) :: maxSplit = 0

    ! Variance reduction options
    logical(defBool)  :: weightWindows = .false., splitting = .false., roulette = .false., implicitAbsorption = .false., &
                         implicitSites = .false., uniFissSites = .false.

    ! Variance reduction requirements
    type(weightWindowsField), pointer :: weightWindowsMap
  contains
    ! Initialisation procedure
    procedure :: init

    ! Implementation of customisable procedures
    procedure :: implicit
    procedure :: fission
    procedure :: cutoffs

    ! Variance reduction procedures
    procedure, private :: split
    procedure, private :: russianRoulette
  end type neutronCEimp

contains

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self, dict)
    class(neutronCEimp), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    integer(shortInt)                  :: idx
    character(*), parameter :: Here = 'init (neutronCEimp_class.f90)'

    ! Call superclass
    call init_super(self, dict)

    ! Obtain settings for variance reduction
    call dict % getOrDefault(self % weightWindows, 'weightWindows', .false.)
    call dict % getOrDefault(self % maxSplit, 'maxSplit', 1000)
    call dict % getOrDefault(self % splitting, 'split', .false.)
    call dict % getOrDefault(self % roulette, 'roulette', .false.)
    call dict % getOrDefault(self % minWgt, 'minWgt', 0.25_defReal)
    call dict % getOrDefault(self % maxWgt, 'maxWgt', 1.25_defReal)
    call dict % getOrDefault(self % avWgt, 'avWgt', 0.5_defReal)
    call dict % getOrDefault(self % implicitAbsorption, 'impAbs', .false.)
    call dict % getOrDefault(self % implicitSites, 'impGen', .true.)
    call dict % getOrDefault(self % uniFissSites, 'UFS', .false.)

    if (self % splitting) then
      if (self % maxWgt < 2 * self % minWgt) call fatalError(Here,&
              'Upper weight bound must be at least twice the lower weight bound')
    end if

    if (self % implicitAbsorption) then
      if (.not.self % roulette .and. .not. self % weightWindows) call fatalError(Here,&
         'Must use Russian roulette or weight windows when using implicit absorption')
      if (.not.self % implicitSites) call fatalError(Here,&
         'Must generate fission sites implicitly when using implicit absorption')
    end if

    ! Sets up the uniform fission sites field
    if (self % uniFissSites) then
      idx = gr_fieldIdx(nameUFS)
      self % ufsField => uniFissSitesField_TptrCast(gr_fieldPtr(idx))
    end if

    ! Sets up the weight windows field
    if (self % weightWindows) then
      idx = gr_fieldIdx(nameWW)
      self % weightWindowsMap => weightWindowsField_TptrCast(gr_fieldPtr(idx))
    end if

  end subroutine init

  !!
  !! Perform implicit treatment
  !!
  subroutine implicit(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEimp), intent(inout)       :: self
    class(physicalParticle), intent(inout)   :: p
    type(tallyAdmin), intent(inout)          :: tally
    type(collisionData), intent(inout)       :: collDat
    class(particleDungeon), intent(inout)    :: thisCycle, nextCycle
    integer(shortInt)                        :: i, n
    logical(defBool)                         :: fiss_and_implicit
    real(defReal)                            :: E, E_max, E_out, mu, k_eff, phi, randomNumber, &
                                                sig_scatter, sig_nufiss, sig_tot, totalElastic, wgt
    real(defReal), dimension(3)              :: r, u, uGlobal, val
    type(CENeutron), pointer                 :: CENeutronPtr
    type(CEParticleState), pointer           :: CEParticleStatePtr
    type(fissionCE), pointer                 :: fission
    type(neutronMicroXSs)                    :: microXSs
    type(RNG), pointer                       :: RNGPtr
    character(*), parameter                  :: Here = 'implicit (neutronCEimp_class.f90)'

    CENeutronPtr => castCENeutronPtr(p, .true.)

    ! Generate fission sites if nuclide is fissile
    fiss_and_implicit = self % getNuclideIsFissile() .and. self % implicitSites

    ! Get pointer to RNG.
    RNGPtr => CENeutronPtr % getRNGPtr()

    if (fiss_and_implicit) then
      ! Obtain required data
      wgt = CENeutronPtr % getWeight()                ! Current weight
      k_eff = CENeutronPtr % getKEff()            ! k_eff for normalisation
      call RNGPtr % generate(randomNumber)     ! Random number to sample sites

      ! Retrieve cross section at the energy used for reaction sampling
      call self % getNuclideMicroXS(collDat % E, collDat % kT, RNGPtr, microXSs)
      sig_nufiss = microXSs % nuFission
      sig_tot    = microXSs % total

      ! Sample number of fission sites generated
      ! Support -ve weight particles
      if (self % uniFissSites) then
        val = self % ufsField % at(p)
        n = int(abs((wgt * sig_nufiss) / (sig_tot * k_eff)) * val(1) / val(2) + randomNumber, shortInt)
        wgt =  val(2) / val(1)

      else
        n = int(abs((wgt * sig_nufiss) / (sig_tot * k_eff)) + randomNumber, shortInt)
        wgt =  sign(ONE, wgt)

      end if

      ! Shortcut particle generation if no particles were sampled
      if (n < 1) return

      ! Get fission Reaction
      fission => fissionCE_TptrCast(self % getReaction(N_FISSION, collDat % nucIdx))
      if (.not. associated(fission)) call fatalError(Here, 'Failed to retrieve fissionCE.')

      ! Store new sites in the next cycle dungeon
      E = CENeutronPtr % getEnergy()
      E_max = self % getMaximumEnergy()
      r = CENeutronPtr % getGlobalPosition()
      uGlobal = CENeutronPtr % getGlobalDirection()
      do i = 1, n
        call fission % sampleOut(mu, phi, E_out, E, RNGPtr)
        u = rotateVector(uGlobal, mu, phi)
        E_out = min(E_out, E_max)

        ! Copy extra detail from parent particle (i.e. time, flags ect.)
        CEParticleStatePtr => castCEParticleStatePtr(CENeutronPtr % updateAndGetCurrentStatePtr(), .true.)

        ! Overwrite position, direction, energy and weight
        call CEParticleStatePtr % setGlobalPosition(r)
        call CEParticleStatePtr % setGlobalDirection(u)
        call CEParticleStatePtr % setEnergy(E_out)
        call CEParticleStatePtr % setWeight(wgt)
        call CEParticleStatePtr % setCollisionsNumber(0)

        call nextCycle % detain(CEParticleStatePtr)
        if (self % uniFissSites) call self % ufsField % storeFS(CEParticleStatePtr)

        ! Report birth of new particle
        call tally % reportSpawn(N_FISSION, CENeutronPtr, CEParticleStatePtr)

      end do

    end if

    ! Perform implicit absorption
    if (self % implicitAbsorption) then
      if (.not. fiss_and_implicit) then
        call self % getNuclideMicroXS(collDat % E, collDat % kT, RNGPtr, microXSs)

      end if

      sig_scatter = microXSs % elasticScatter + microXSs % inelasticScatter
      sig_tot = microXSs % total
      call CENeutronPtr % setWeight(CENeutronPtr % getWeight() * sig_scatter / sig_tot)
      
      ! Sample between elastic and inelastic
      totalElastic = microXSs % elasticScatter + microXSs % inelasticScatter
      call RNGPtr % generate(randomNumber)
      collDat % MT = merge(N_N_elastic, N_N_inelastic, randomNumber < microXSs % elasticScatter / totalElastic)

    end if

  end subroutine implicit

  !!
  !! Process fission reaction
  !!
  subroutine fission(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEimp), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon), intent(inout)  :: thisCycle, nextCycle
    integer(shortInt)                      :: n, i
    real(defReal)                          :: E, E_max, E_out, k_eff, mu, phi, randomNumber, &
                                              sig_nufiss, sig_fiss, wgt
    real(defReal), dimension(3)            :: r, u, uGlobal, val
    type(CENeutron), pointer               :: CENeutronPtr
    type(CEParticleState), pointer         :: CEParticleStatePtr
    type(fissionCE), pointer               :: fiss
    type(neutronMicroXSs)                  :: microXSs
    type(RNG), pointer                     :: RNGPtr
    character(*), parameter                :: Here = 'fission (neutronCEimp_class.f90)'

    CENeutronPtr => castCENeutronPtr(p, .true.)
    call CENeutronPtr % setIsDead(.true.)
    RNGPtr => CENeutronPtr % getRNGPtr()
    if (.not. self % implicitSites) then
      ! Obtain required data
      wgt = CENeutronPtr % getWeight()                ! Current weight
      k_eff = CENeutronPtr % getKEff()            ! k_eff for normalisation
      call RNGPtr % generate(randomNumber)     ! Random number to sample sites

      ! Retrieve cross section at the energy used for reaction sampling
      call self % getNuclideMicroXS(collDat % E, collDat % kT, RNGPtr, microXSs)
      sig_nufiss = microXSs % nuFission
      sig_fiss = microXSs % fission

      ! Sample number of fission sites generated
      ! Support -ve weight particles
      ! Note change of denominator (sig_fiss) wrt implicit generation
      if (self % uniFissSites) then
        val = self % ufsField % at(p)
        n = int(abs((wgt * sig_nufiss) / (sig_fiss * k_eff)) * val(1) / val(2) + randomNumber, shortInt)
        wgt =  val(2) / val(1)

      else
        n = int(abs((wgt * sig_nufiss) / (sig_fiss * k_eff)) + randomNumber, shortInt)
        wgt = sign(ONE, wgt)

      end if

      ! Shortcut particle generation if no particles were sampled
      if (n < 1) return

      ! Get fission Reaction
      fiss => fissionCE_TptrCast(self % getReaction(N_FISSION, collDat % nucIdx))
      if (.not. associated(fiss)) call fatalError(Here, 'Failed to retrieve fissionCE.')

      ! Store new sites in the next cycle dungeon
      E = CENeutronPtr % getEnergy()
      E_max = self % getMaximumEnergy()
      r = CENeutronPtr % getGlobalPosition()
      uGlobal = CENeutronPtr % getGlobalDirection()
      do i = 1, n
        call fiss % sampleOut(mu, phi, E_out, E, RNGPtr)
        u = rotateVector(uGlobal, mu, phi)
        E_out = min(E_out, E_max)

        ! Copy extra detail from parent particle (i.e. time, flags ect.)
        CEParticleStatePtr => castCEParticleStatePtr(CENeutronPtr % updateAndGetCurrentStatePtr(), .true.)

        ! Overwrite position, direction, energy and weight
        call CEParticleStatePtr % setGlobalPosition(r)
        call CEParticleStatePtr % setGlobalDirection(u)
        call CEParticleStatePtr % setEnergy(E)
        call CEParticleStatePtr % setWeight(wgt)
        call CEParticleStatePtr % setCollisionsNumber(0)

        call nextCycle % detain(CEParticleStatePtr)
        if (self % uniFissSites) call self % ufsField % storeFS(CEParticleStatePtr)

        ! Report birth of new particle
        call tally % reportSpawn(N_FISSION, CENeutronPtr, CEParticleStatePtr)

      end do

    end if

  end subroutine fission

  !!
  !! Apply cutoffs
  !!
  subroutine cutoffs(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEimp), intent(inout)     :: self
    class(physicalParticle), intent(inout) :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon), intent(inout)  :: thisCycle, nextCycle
    real(defReal)                          :: avWgt, maxWgt, minWgt, weight
    real(defReal), dimension(3)            :: val
    type(CENeutron), pointer               :: CENeutronPtr

    CENeutronPtr => castCENeutronPtr(p, .true.)
    weight = CENeutronPtr % getWeight()
    if (CENeutronPtr % getIsDead()) then
      ! Do nothing !

    elseif (CENeutronPtr % getEnergy() < self % getMinimumEnergy()) then
      call CENeutronPtr % setIsDead(.true.)

    ! Weight Windows treatment
    elseif (self % weightWindows) then
      val = self % weightWindowsMap % at(CENeutronPtr)
      minWgt = val(1)
      maxWgt = val(2)
      avWgt  = val(3)

      ! If a particle is outside the WW map and all the weight limits
      ! are zero nothing happens. NOTE: this holds for positive weights only
      if (maxWgt < weight .and. maxWgt /= ZERO .and. CENeutronPtr % getSplitsNumber() < self % maxSplit) then
        call self % split(CENeutronPtr, tally, thisCycle, maxWgt)

      elseif (weight < minWgt) then
        call self % russianRoulette(CENeutronPtr, avWgt)

      end if

    ! Splitting with fixed threshold
    elseif (self % splitting .and. self % maxWgt < weight) then
      call self % split(CENeutronPtr, tally, thisCycle, self % maxWgt)

    ! Roulette with fixed threshold and survival weight
    elseif (self % roulette .and. weight < self % minWgt) then
      call self % russianRoulette(CENeutronPtr, self % avWgt)

    end if

  end subroutine cutoffs

  !!
  !! Perform Russian roulette on a particle
  !!
  subroutine russianRoulette(self, n_CE, avWgt)
    class(neutronCEimp), intent(inout) :: self
    class(CENeutron), intent(inout)    :: n_CE
    real(defReal), intent(in)          :: avWgt
    real(defReal)                      :: randomNumber
    type(RNG), pointer                 :: RNGPtr

    RNGPtr => n_CE % getRNGPtr()
    call RNGPtr % generate(randomNumber)
    if (randomNumber < (ONE - n_CE % getWeight() / avWgt)) then
      call n_CE % setIsDead(.true.)

    else
      call n_CE % setWeight(avWgt)

    end if

  end subroutine russianRoulette

  !!
  !! Split particle which has too large a weight
  !!
  subroutine split(self, n_CE, tally, thisCycle, maxWgt)
    class(neutronCEimp), intent(inout)    :: self
    class(CENeutron), intent(inout)       :: n_CE
    type(tallyAdmin), intent(inout)       :: tally
    class(particleDungeon), intent(inout) :: thisCycle
    real(defReal), intent(in)             :: maxWgt
    integer(shortInt)                     :: i, mult, nSplits
    real(defReal)                         :: newWeight, weight
    type(CEParticleState), pointer        :: CEParticleStatePtr

    nSplits = n_CE % getSplitsNumber()
    weight = n_CE % getWeight()
    
    ! This value must be at least 2
    mult = ceiling(weight / maxWgt)

    ! Limit maximum split and compute newWeight.
    if (self % maxSplit < mult + nSplits) mult = self % maxSplit - nSplits + 1
    newWeight = weight / mult

    ! Copy particle to a particle state
    ! Note that particleState doesn't have property splitCount, so it is reset
    ! to 0 for the new particle
    CEParticleStatePtr => castCEParticleStatePtr(n_CE % updateAndGetCurrentStatePtr(), .true.)
    call CEParticleStatePtr % setWeight(newWeight)

    ! Add split particle's to the dungeon
    do i = 1, mult - 1
      call thisCycle % detain(CEParticleStatePtr)
      call tally % reportSpawn(N_N_SPLIT, n_CE, CEParticleStatePtr)

    end do
    ! Update particle split count and weight.
    call n_CE % setSplitsNumber(nSplits + mult)
    call n_CE % setWeight(newWeight)

  end subroutine split

end module neutronCEimp_class