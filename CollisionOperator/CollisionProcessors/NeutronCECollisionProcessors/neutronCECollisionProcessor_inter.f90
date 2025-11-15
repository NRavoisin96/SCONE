module neutronCECollisionProcessor_inter

  use CECollisionData_class,        only : castCECollisionDataPtr, CECollisionData
  use CENeutron_class,              only : castCENeutronPtr, CENeutron
  use ceNeutronDatabase_inter,      only : ceNeutronDatabase, ceNeutronDatabase_CptrCast
  use ceNeutronMaterial_class,      only : ceNeutronMaterial, ceNeutronMaterial_CptrCast
  use ceNeutronNuclide_inter,       only : ceNeutronNuclide, ceNeutronNuclide_CptrCast
  use CENeutronState_class,         only : CENeutronState
  use CEParticleState_class,        only : buildCEParticleStatePayload, castBuildCEParticleStatePayloadPtr
  use collisionData_class,          only : collisionData
  use collisionProcessor_inter,     only : collisionProcessor, init_super => init, kill_super => kill
  use dictionary_class,             only : dictionary
  use endfConstants,                only : N_FISSION, N_N_ELASTIC, N_N_ThermEL
  use errors_mod,                   only : fatalError
  use fissionCE_class,              only : fissionCE, fissionCE_TptrCast
  use genericProcedures,            only : rotateVector
  use neutronXsPackages_class,      only : neutronMicroXSs
  use nuclearDataReg_mod,           only : getNeutronCE
  use numPrecision
  use particleDungeon_class,        only : particleDungeon
  use physicalParticle_inter,       only : physicalParticle
  use physicalParticleState_class,  only : physicalParticleState
  use reactionHandle_inter,         only : reactionHandle
  use RNG_class,                    only : RNG
  use scalarField_inter,            only : getTemperatureFieldPtr, scalarField
  use scatteringKernels_func,       only : asymptoticInelasticScatter, asymptoticScatter, targetVelocity_constXS, &
                                           targetVelocity_DBRCXS
  use tallyAdmin_class,             only : tallyAdmin
  use transportObjectState_class,   only : buildTransportObjectStatePayload
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE, uncorrelatedReactionCE_CptrCast
  use universalVariables

  implicit none
  private

  ! Public procedures.
  public :: cutoffs, init, kill

  !!
  !!
  !!
  type, public, abstract, extends(collisionProcessor) :: neutronCECollisionProcessor
    private
    class(ceNeutronNuclide), pointer :: nuc => null()
    ! Settings.
    real(defReal) :: DBRCeMax = ZERO, DBRCeMin = ZERO, maxE = ZERO, minE = ZERO, threshA = ZERO, threshE = ZERO
  contains
    procedure          :: allocateCollisionData
    procedure          :: cutoffs
    procedure          :: getFissionHandlePtr
    procedure          :: getMaximumEnergy
    procedure          :: getMinimumEnergy
    procedure          :: getNuclideIsFissile
    procedure          :: getNuclideMicroXS
    procedure          :: elastic
    procedure          :: inelastic
    procedure          :: init
    procedure          :: kill
    procedure          :: newPhysicalParticleState
    procedure          :: sampleCollision
    procedure, private :: scatterFromFixed
    procedure, private :: scatterFromMoving
    procedure, private :: scatterInLAB
  end type neutronCECollisionProcessor

contains
  !!
  !!
  !!
  subroutine allocateCollisionData(self, collDat)
    class(neutronCECollisionProcessor), intent(in) :: self
    class(collisionData), allocatable, intent(out) :: collDat

    allocate(CECollisionData :: collDat)

  end subroutine allocateCollisionData

  !!
  !!
  !!
  subroutine cutoffs(self, p, dungeon, tally)
    class(neutronCECollisionProcessor), intent(in) :: self
    class(physicalParticle), intent(inout)         :: p
    type(particleDungeon), intent(inout)           :: dungeon
    type(tallyAdmin), intent(inout)                :: tally
    type(CENeutron), pointer                       :: CENeutronPtr

    ! Kill continuous energy neutron if its energy is below threshold.
    CENeutronPtr => castCENeutronPtr(p)
    if (CENeutronPtr % getEnergy() < self % minE) call CENeutronPtr % setIsDead(.true.)

  end subroutine cutoffs

  !!
  !!
  !!
  function getFissionHandlePtr(self, collDat) result(fissionHandlePtr)
    class(neutronCECollisionProcessor), intent(in) :: self
    class(collisionData), intent(in)               :: collDat
    class(reactionHandle), pointer                 :: fissionHandlePtr

    fissionHandlePtr => fissionCE_TptrCast(self % getReaction(N_FISSION, collDat % nucIdx))

  end function getFissionHandlePtr

  !!
  !!
  !!
  elemental function getMaximumEnergy(self) result(maximumEnergy)
    class(neutronCECollisionProcessor), intent(in) :: self
    real(defReal)                                  :: maximumEnergy

    maximumEnergy = self % maxE

  end function getMaximumEnergy

  !!
  !!
  !!
  elemental function getMinimumEnergy(self) result(minimumEnergy)
    class(neutronCECollisionProcessor), intent(in) :: self
    real(defReal)                                  :: minimumEnergy

    minimumEnergy = self % minE

  end function getMinimumEnergy

  !!
  !!
  !!
  elemental function getNuclideIsFissile(self) result(isFissile)
    class(neutronCECollisionProcessor), intent(in) :: self
    logical(defBool)                               :: isFissile

    isFissile = self % nuc % isFissile()

  end function getNuclideIsFissile

  !!
  !!
  !!
  subroutine getNuclideMicroXS(self, E, kT, rand, microXS)
    class(neutronCECollisionProcessor), intent(in) :: self
    real(defReal), intent(in)                      :: E, kT
    class(RNG), intent(inout)                      :: rand
    type(neutronMicroXSs), intent(out)             :: microXS

    call self % nuc % getMicroXSs(E, kT, microXS, rand)

  end subroutine getNuclideMicroXS

  !!
  !! Process elastic scattering
  !!
  !! All CE elastic scattering happens in the CM frame
  !!
  subroutine elastic(self, collDat, p)
    class(neutronCECollisionProcessor), intent(in) :: self
    class(collisionData), intent(inout)            :: collDat
    class(physicalParticle), intent(inout)         :: p
    class(ceNeutronMaterial), pointer              :: CENeutronMaterialPtr
    class(uncorrelatedReactionCE), pointer         :: reac
    logical(defBool)                               :: isFixed, hasDBRC
    type(CECollisionData), pointer                 :: CECollisionDataPtr
    type(CENeutron), pointer                       :: CENeutronPtr
    character(*), parameter                        :: HERE = 'elastic (neutronCECollisionProcessor_inter.f90)'

    ! Downcast physicalParticle to CENeutron and retrieve neutron energy.
    CECollisionDataPtr => castCECollisionDataPtr(collDat)
    CENeutronPtr => castCENeutronPtr(p)
    
    ! Assess if thermal scattering data is needed or not
    if (self % nuc % needsSabEl(CECollisionDataPtr % initialEnergy)) CECollisionDataPtr % MT = N_N_ThermEL

    ! Get reaction
    reac => uncorrelatedReactionCE_CptrCast(self % getReaction(CECollisionDataPtr % MT, CECollisionDataPtr % nucIdx))
    if (.not. associated(reac)) call fatalError(HERE, 'Failed to get elastic neutron scatter.')

    ! Scatter particle
    CECollisionDataPtr % A =  self % nuc % getMass()

    ! Retrieve kT from either material or nuclide
    CENeutronMaterialPtr => ceNeutronMaterial_CptrCast(self % getMaterialPtr(CECollisionDataPtr % matIdx))
    if (.not. associated(CENeutronMaterialPtr)) call fatalError(HERE, 'Failed to retrieve CE neutron material.')
    if (.not. CENeutronMaterialPtr % useTMS(CECollisionDataPtr % initialEnergy)) CECollisionDataPtr % kT = self % nuc % getkT()

    ! Check if DBRC is on.
    hasDBRC = self % nuc % hasDBRC()
    isFixed = CECollisionDataPtr % kT * self % threshE < CECollisionDataPtr % initialEnergy .and. &
              self % threshA < CECollisionDataPtr % A .and. .not. hasDBRC

    ! Apply criterion for Free-Gas vs Fixed Target scattering
    if (.not. reac % inCMFrame()) then
      call self % scatterInLAB(CENeutronPtr, CECollisionDataPtr, reac)

    elseif (isFixed) then
      call self % scatterFromFixed(CENeutronPtr, CECollisionDataPtr, reac)

    else
      call self % scatterFromMoving(CENeutronPtr, CECollisionDataPtr, reac)

    end if

  end subroutine elastic

  !!
  !! Process inelastic scattering
  !!
  subroutine inelastic(self, collDat, p)
    class(neutronCECollisionProcessor), intent(in) :: self
    class(collisionData), intent(inout)            :: collDat
    class(physicalParticle), intent(inout)         :: p
    class(uncorrelatedReactionCE), pointer         :: reac
    type(CECollisionData), pointer                 :: CECollisionDataPtr
    type(CENeutron), pointer                       :: CENeutronPtr
    character(*), parameter                        :: Here = 'inelastic (neutronCECollisionProcessor_inter.f90)'

    ! Invert inelastic scattering and get reaction
    CECollisionDataPtr => castCECollisionDataPtr(collDat)
    CENeutronPtr => castCENeutronPtr(p)
    CECollisionDataPtr % MT = self % nuc % invertInelastic(CECollisionDataPtr % reactionEnergy, CECollisionDataPtr % RNGPtr)
    reac => uncorrelatedReactionCE_CptrCast(self % getReaction(CECollisionDataPtr % MT, CECollisionDataPtr % nucIdx))
    if (.not. associated(reac)) call fatalError(Here, 'Failed to retrieve scattering reaction.')

    ! Scatter particle
    if (reac % inCMFrame()) then
      CECollisionDataPtr % A =  self % nuc % getMass()
      call self % scatterFromFixed(CENeutronPtr, CECollisionDataPtr, reac)

    else
      call self % scatterInLAB(CENeutronPtr, CECollisionDataPtr, reac)

    end if

    ! Apply weigth change
    call CENeutronPtr % setWeight(CENeutronPtr % getWeight() * reac % release(CENeutronPtr % getEnergy()))

  end subroutine inelastic

  !!
  !!
  !!
  subroutine init(self, dict)
    class(neutronCECollisionProcessor), intent(inout) :: self
    class(dictionary), intent(in)                     :: dict
    character(*), parameter                           :: HERE = 'init (neutronCECollisionProcessor_inter.f90)'

    ! Initialise superclass.
    call init_super(self, dict)

    ! Load nuclear data pointer.
    call self % setNuclearDatabasePtr(getNeutronCE())

    ! Get settings from dictionary.
    call dict % getOrDefault(self % minE, 'minEnergy', 1.0E-11_defReal)
    call dict % getOrDefault(self % maxE, 'maxEnergy', 20.0_defReal)
    call dict % getOrDefault(self % threshE, 'energyThreshold', 400.0_defReal)
    call dict % getOrDefault(self % threshA, 'massThreshold', 1.0_defReal)
    call dict % getOrDefault(self % DBRCeMin, 'DBRCeMin', 1.0E-8_defReal)
    call dict % getOrDefault(self % DBRCeMax, 'DBRCeMax', 200.0E-6_defReal)

    ! Verify settings.
    if (self % minE < ZERO) call fatalError(HERE, 'Minimum energy must be positive.')
    if (self % maxE < ZERO) call fatalError(HERE, 'Maximum energy must be positive.')
    if (self % maxE <= self % minE) call fatalError(HERE, 'maxEnergy <= minEnergy.')
    if (self % threshE < ZERO) call fatalError(HERE,' Energy threshold must be positive.')
    if (self % threshA < ZERO) call fatalError(HERE, 'Mass threshold must be positive.')

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(neutronCECollisionProcessor), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % nuc => null()
    self % DBRCeMax = ZERO
    self % DBRCeMin = ZERO
    self % maxE = ZERO
    self % minE = ZERO
    self % threshA = ZERO
    self % threshE = ZERO

  end subroutine kill

  !!
  !!
  !!
  subroutine newPhysicalParticleState(self, collDat, fissionHandlePtr, payload, newState)
    class(neutronCECollisionProcessor), intent(in)         :: self
    class(collisionData), intent(in)                       :: collDat
    class(reactionHandle), pointer, intent(in)             :: fissionHandlePtr
    class(buildTransportObjectStatePayload), intent(inout) :: payload
    class(physicalParticleState), allocatable, intent(out) :: newState
    type(buildCEParticleStatePayload), pointer             :: payloadPtr
    type(CECollisionData), pointer                         :: CECollisionDataPtr
    type(fissionCE), pointer                               :: fissionCEPtr
    real(defReal)                                          :: mu, phi

    ! Downcast classes to correct types.
    payloadPtr => castBuildCEParticleStatePayloadPtr(payload)
    CECollisionDataPtr => castCECollisionDataPtr(collDat)
    fissionCEPtr => fissionCE_TptrCast(fissionHandlePtr)

    ! Allocate new state.
    allocate(CENeutronState :: newState)

    ! Finalise payload then initialise state.
    call fissionCEPtr % sampleOut(mu, phi, payloadPtr % energy, CECollisionDataPtr % initialEnergy, CECollisionDataPtr % RNGPtr)
    payloadPtr % uGlobal = rotateVector(CECollisionDataPtr % u, mu, phi)
    payloadPtr % energy = min(payloadPtr % energy, self % maxE)
    call newState % init(payload)

  end subroutine newPhysicalParticleState

  !!
  !! Samples collision without any implicit treatment
  !!
  subroutine sampleCollision(self, p, collDat)
    class(neutronCECollisionProcessor), intent(inout) :: self
    class(physicalParticle), intent(in)               :: p
    class(collisionData), intent(inout)               :: collDat
    class(ceNeutronMaterial), pointer                 :: CENeutronMaterialPtr
    class(scalarField), pointer                       :: temperatureFieldPtr
    real(defReal)                                     :: randomNumber, temperature
    type(CECollisionData), pointer                    :: CECollisionDataPtr
    type(neutronMicroXSs)                             :: microXSs
    character(*), parameter                           :: Here = 'sampleCollision (neutronCECollisionProcessor_inter.f90)'

    ! Downcast collisionData and physicalParticle to correct types.
    CECollisionDataPtr => castCECollisionDataPtr(collDat)

    ! Load material pointer.
    CENeutronMaterialPtr => ceNeutronMaterial_CptrCast(self % getMaterialPtr(CECollisionDataPtr % matIdx))
    call self % setCurrentMaterialPtr(CENeutronMaterialPtr)

    ! Retrieve material temperature from temperature field.
    CECollisionDataPtr % kT = CENeutronMaterialPtr % kT
    temperatureFieldPtr => getTemperatureFieldPtr()
    if (associated (temperatureFieldPtr)) then
      temperature = temperatureFieldPtr % at(p % getCoordsPtr())
      if (ZERO < temperature) CECollisionDataPtr % kT = kBoltzmann * temperature / joulesPerMeV

    end if

    ! Select collision nuclide.
    call CENeutronMaterialPtr % sampleNuclide(CECollisionDataPtr % initialEnergy, CECollisionDataPtr % kT, &
                                              CECollisionDataPtr % RNGPtr, CECollisionDataPtr % nucIdx, &
                                              CECollisionDataPtr % reactionEnergy)

    ! If nuclide was rejected in TMS loop return to tracking
    if (CECollisionDataPtr % nucIdx == REJECTED) return

    self % nuc => ceNeutronNuclide_CptrCast(self % getNuclidePtr(CECollisionDataPtr % nucIdx))
    if (.not. associated(self % nuc)) call fatalError(Here, 'Failed to retrieve CE Neutron Nuclide')

    ! Select Main reaction channel
    call self % nuc % getMicroXSs(CECollisionDataPtr % reactionEnergy, CECollisionDataPtr % kT, &
                                  microXSs, CECollisionDataPtr % RNGPtr)
    call CECollisionDataPtr % RNGPtr % generate(randomNumber)
    CECollisionDataPtr % MT = microXss % invert(randomNumber)
    CECollisionDataPtr % sigma_elasticScatter = microXSs % elasticScatter
    CECollisionDataPtr % sigma_fission = microXSs % fission
    CECollisionDataPtr % sigma_inelasticScatter = microXSs % inelasticScatter
    CECollisionDataPtr % sigma_nuFiss = microXSs % nuFission
    CECollisionDataPtr % sigma_tot = microXSs % total

  end subroutine sampleCollision

  !!
  !! Subroutine to perform scattering from stationary target.
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterFromFixed(self, n_CE, CECollDat, reac)
    class(neutronCECollisionProcessor), intent(in) :: self
    type(CENeutron), intent(inout)                 :: n_CE
    type(CECollisionData), intent(inout)           :: CECollDat
    class(uncorrelatedReactionCE), intent(in)      :: reac
    real(defReal)                                  :: E_outCM, phi

    ! Sample mu, phi and outgoing energy
    call reac % sampleOut(CECollDat % muL, phi, E_outCM, CECollDat % initialEnergy, CECollDat % RNGPtr)

    ! Save incident energy
    CECollDat % finalEnergy = CECollDat % initialEnergy

    if (CECollDat % MT == N_N_ELASTIC) then
      call asymptoticScatter(CECollDat % finalEnergy, CECollDat % muL, CECollDat % A)

    else
      call asymptoticInelasticScatter(CECollDat % finalEnergy, CECollDat % muL, E_outCM, CECollDat % A)

    end if

    ! Update particle state
    call n_CE % rotate(CECollDat % muL, phi)
    call n_CE % setEnergy(CECollDat % finalEnergy)

  end subroutine scatterFromFixed

  !!
  !! Subroutine to perform scattering from moving target
  !! Supports only elastic collisions
  !!
  subroutine scatterFromMoving(self, n_CE, CECollDat, reac)
    class(neutronCECollisionProcessor), intent(in) :: self
    type(CENeutron), intent(inout)                 :: n_CE
    type(CECollisionData), intent(inout)           :: CECollDat
    class(uncorrelatedReactionCE), intent(in)      :: reac
    class(ceNeutronDatabase), pointer              :: CENeutronDatabasePtr
    class(ceNeutronNuclide), pointer               :: ceNuc0K
    logical(defBool)                               :: inEnergyRange, hasDBRC
    real(defReal)                                  :: dummy, maj, mu, phi, speed
    real(defReal), dimension(3)                    :: dir_post, dir_pre, V_cm, V_n, v_t
    character(*), parameter                        :: HERE = 'scatterFromMoving (neutronCECollisionProcessor_inter.f90)'

    ! Get neutron direction and velocity
    dir_pre = n_CE % getGlobalDirection()
    V_n = dir_pre * sqrt(CECollDat % initialEnergy)

    ! Sample target velocity with constant XS or with DBRC
    ! Check energy range
    inEnergyRange = CECollDat % initialEnergy <= self % DBRCeMax .and. self % DBRCeMin <= CECollDat % initialEnergy
    
    ! Check if DBRC is on for this target nuclide
    hasDBRC = self % nuc % hasDBRC()
    if (inEnergyRange .and. hasDBRC) then
      ! Retrieve and downcast nuclear database pointer.
      CENeutronDatabasePtr => ceNeutronDatabase_CptrCast(self % getNuclearDatabasePtr())
      if (.not. associated(CENeutronDatabasePtr)) &
      call fatalError(HERE, 'Failed to retrieve active CE neutron database.')

      ! Retrieve 0K nuclide index from DBRC nuclide map
      CECollDat % nucIdx = CENeutronDatabasePtr % mapDBRCnuc % get(CECollDat % nucIdx)

      ! Assign pointer for the 0K nuclide
      ceNuc0K => ceNeutronNuclide_CptrCast(CENeutronDatabasePtr % getNuclide(CECollDat % nucIdx))
      if (.not. associated(ceNuc0K)) call fatalError(HERE, 'Failed to retrieve CE neutron nuclide.')

      ! Get elastic scattering 0K majorant
      maj = CENeutronDatabasePtr % getScattMicroMajXS(CECollDat % initialEnergy, CECollDat % kT, &
                                                      CECollDat % A, CECollDat % nucIdx)

      ! Use DBRC to sample target velocity
      V_t = targetVelocity_DBRCXS(ceNuc0K, CECollDat % initialEnergy, dir_pre, CECollDat % A, &
                                  CECollDat % kT, CECollDat % RNGPtr, maj)

    else
      ! Constant cross section approximation
      V_t = targetVelocity_constXS(CECollDat % initialEnergy, dir_pre, CECollDat % A, CECollDat % kT, CECollDat % RNGPtr)

    end if

    ! Calculate Centre-of-Mass velocity
    V_cm = (V_n + V_t * CECollDat % A) / (CECollDat % A + ONE)

    ! Move Neutron velocity to CM frame, store speed and calculate new normalised direction
    V_n = V_n - V_cm
    speed = norm2(V_n)
    V_n = V_n / speed

    ! Sample mu and phi in CM frame
    call reac % sampleOut(mu, phi, dummy, CECollDat % initialEnergy, CECollDat % RNGPtr)

    ! Obtain post collision speed
    V_n = rotateVector(V_n, mu, phi) * speed

    ! Return to LAB frame
    V_n = V_n + V_cm

    ! Calculate new neutron speed and direction
    speed = norm2(V_n)
    dir_post = V_n / speed

    ! Update particle state and calculate mu in LAB frame
    call n_CE % setEnergy(speed * speed)
    call n_CE % point(dir_post)
    CECollDat % muL = dot_product(dir_pre, dir_post)

  end subroutine scatterFromMoving

  !!
  !! Subroutine to perform scattering in LAB frame
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterInLAB(self, n_CE, CECollDat, reac)
    class(neutronCECollisionProcessor), intent(in) :: self
    type(CENeutron), intent(inout)                 :: n_CE
    type(CECollisionData), intent(inout)           :: CECollDat
    class(uncorrelatedReactionCE), intent(in)      :: reac
    real(defReal)                                  :: phi

    ! Sample scattering angles and post-collision energy
    call reac % sampleOut(CECollDat % muL, phi, CECollDat % finalEnergy, CECollDat % initialEnergy, CECollDat % RNGPtr)

    ! Update neutron state
    call n_CE % setEnergy(CECollDat % finalEnergy)
    call n_CE % rotate(CECollDat % muL, phi)

  end subroutine scatterInLAB

end module neutronCECollisionProcessor_inter