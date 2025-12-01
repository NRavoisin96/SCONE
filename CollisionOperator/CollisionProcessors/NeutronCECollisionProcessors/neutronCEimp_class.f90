module neutronCEimp_class

  use CECollisionData_class,             only : castCECollisionDataPtr, CECollisionData
  use CENeutron_class,                   only : castCENeutronPtr, CENeutron
  use CEParticleState_class,             only : castCEParticleStatePtr, CEParticleState
  use collisionData_class,               only : collisionData
  use collisionProcessor_inter,          only : fission_super => fission, implicit_super => implicit
  use dictionary_class,                  only : dictionary
  use endfConstants
  use errors_mod,                        only : fatalError
  use fissionCE_class,                   only : fissionCE, fissionCE_TptrCast
  use genericProcedures,                 only : numToChar, rotateVector
  use geometryReg_mod,                   only : gr_fieldIdx => fieldIdx, gr_fieldPtr => fieldPtr
  use neutronCECollisionProcessor_inter, only : cutoffs_super => cutoffs, init_super => init, &
                                                kill_super => kill, neutronCECollisionProcessor
  use neutronXsPackages_class,           only : neutronMicroXSs
  use numPrecision
  use particleDungeon_class,             only : particleDungeon
  use physicalParticle_inter,            only : physicalParticle
  use populationComber_class,            only : populationComber
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
    !! Nuclear Data block pointer
    class(uniFissSitesField), pointer :: ufsField => null()
    logical(defBool)                  :: implicitAbsorption = .false., implicitSites = .false., uniFissSites = .false.
    type(populationComber)            :: comber
  contains
    procedure :: computeNumberOfSecondaryParticles
    procedure :: cutoffs
    procedure :: fission
    procedure :: implicit
    procedure :: init
    procedure :: getImplicitCondition
    procedure :: kill
  end type neutronCEimp

contains
  !!
  !!
  !!
  subroutine computeNumberOfSecondaryParticles(self, factor, collDat, p)
    class(neutronCEimp), intent(in)        :: self
    real(defReal), intent(in)              :: factor
    class(collisionData), intent(inout)    :: collDat
    class(physicalParticle), intent(inout) :: p
    real(defReal)                          :: n, randomNumber
    real(defReal), dimension(3)            :: ufsFieldValues

    ! Sample number of fission sites generated. Supports particles with negative weights.
    n = abs(collDat % weight * collDat % sigma_nuFiss / (factor * collDat % k_eff))
    call collDat % RNGPtr % generate(randomNumber)
    if (self % uniFissSites) then
      ufsFieldValues = self % ufsField % at(p)
      collDat % n = int(n * ufsFieldValues(1) / ufsFieldValues(2) + randomNumber, shortInt)
      collDat % weight = ufsFieldValues(2) / ufsFieldValues(1)

    else
      collDat % n = int(n + randomNumber, shortInt)
      collDat % weight = sign(ONE, collDat % weight)

    end if

  end subroutine computeNumberOfSecondaryParticles

  !!
  !! Apply cutoffs
  !!
  subroutine cutoffs(self, p, dungeon, tally)
    class(neutronCEimp), intent(in)        :: self
    class(physicalParticle), intent(inout) :: p
    type(particleDungeon), intent(inout)   :: dungeon
    type(tallyAdmin), intent(inout)        :: tally

    ! Call superclass procedure and return if particle was killed.
    call cutoffs_super(self, p, dungeon, tally)
    if (p % getIsDead()) return

    ! Call population comber.
    call self % comber % cutoffs(p, dungeon, tally)

  end subroutine cutoffs

  !!
  !! Process fission reaction
  !!
  subroutine fission(self, collDat, p, dungeon, tally)
    class(neutronCEimp), intent(in)        :: self
    class(collisionData), intent(inout)    :: collDat
    class(physicalParticle), intent(inout) :: p
    type(particleDungeon), intent(inout)   :: dungeon
    type(tallyAdmin), intent(inout)        :: tally

    ! Call superclass procedure and return immediately if using implicit sites.
    call fission_super(self, collDat, p, dungeon, tally)
    call implicit_super(self, .not. self % implicitSites, collDat % sigma_fission, collDat, p, dungeon, tally)

  end subroutine fission

  !!
  !! Perform implicit treatment
  !!
  subroutine implicit(self, implicitCondition, factor, collDat, p, dungeon, tally)
    class(neutronCEimp), intent(in)        :: self
    logical(defBool), intent(in)           :: implicitCondition
    real(defReal), intent(in)              :: factor
    class(collisionData), intent(inout)    :: collDat
    class(physicalParticle), intent(inout) :: p
    type(particleDungeon), intent(inout)   :: dungeon
    type(tallyAdmin), intent(inout)        :: tally
    real(defReal)                          :: randomNumber, sigma_totalScatter
    type(CECollisionData), pointer         :: CECollisionDataPtr
    type(CENeutron), pointer               :: CENeutronPtr

    ! Generate fission sites if nuclide is fissile
    call implicit_super(self, self % getImplicitCondition(), collDat % sigma_tot, collDat, p, dungeon, tally)

    ! Perform implicit absorption
    if (self % implicitAbsorption) then
      CECollisionDataPtr => castCECollisionDataPtr(collDat)
      CENeutronPtr => castCENeutronPtr(p)
      
      sigma_totalScatter = collDat % sigma_elasticScatter + collDat % sigma_inelasticScatter
      call CENeutronPtr % setWeight(CENeutronPtr % getWeight() * sigma_totalScatter / collDat % sigma_tot)
      
      ! Sample between elastic and inelastic
      call collDat % RNGPtr % generate(randomNumber)
      CECollisionDataPtr % MT = merge(N_N_elastic, N_N_inelastic, &
                                      randomNumber < collDat % sigma_elasticScatter / sigma_totalScatter)

    end if

  end subroutine implicit

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self, dict)
    class(neutronCEimp), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    character(*), parameter            :: HERE = 'init (neutronCEimp_class.f90)'

    ! Initialise superclass and population comber.
    call init_super(self, dict)
    call self % comber % init(dict)

    ! Obtain settings for variance reduction
    call dict % getOrDefault(self % implicitAbsorption, 'impAbs', .false.)
    call dict % getOrDefault(self % implicitSites, 'impGen', .true.)
    call dict % getOrDefault(self % uniFissSites, 'UFS', .false.)

    if (self % implicitAbsorption) then
      if (.not. (self % comber % getUsesRussianRoulette() .or. self % comber % getUsesWeightWindows())) &
      call fatalError(HERE, 'Must use Russian roulette or weight windows when using implicit absorption.')
      
      if (.not. self % implicitSites) &
      call fatalError(HERE, 'Must generate fission sites implicitly when using implicit absorption.')

    end if

    ! Sets up the uniform fission sites field
    if (self % uniFissSites) self % ufsField => uniFissSitesField_TptrCast(gr_fieldPtr(gr_fieldIdx(nameUFS)))

  end subroutine init

  !!
  !!
  !!
  elemental function getImplicitCondition(self) result(isIt)
    class(neutronCEimp), intent(in) :: self
    logical(defBool)                :: isIt

    isIt = self % getNuclideIsFissile() .and. self % implicitSites

  end function getImplicitCondition

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(neutronCEimp), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % ufsField => null()
    self % implicitAbsorption = .false.
    self % implicitSites = .false.
    self % uniFissSites = .false.
    call self % comber % kill()

  end subroutine kill

end module neutronCEimp_class