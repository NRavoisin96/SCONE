module neutronCECollisionProcessor_inter

  use CENeutron_class,              only : castCENeutronPtr, CENeutron
  use ceNeutronDatabase_inter,      only : ceNeutronDatabase
  use ceNeutronMaterial_class,      only : ceNeutronMaterial, ceNeutronMaterial_CptrCast
  use ceNeutronNuclide_inter,       only : ceNeutronNuclide, ceNeutronNuclide_CptrCast
  use collisionProcessor_inter,     only : collisionData, collisionProcessor, init_super => init
  use dictionary_class,             only : dictionary
  use endfConstants
  use genericProcedures,            only : fatalError, rotateVector
  use neutronXsPackages_class,      only : neutronMicroXSs
  use nuclearDataReg_mod,           only : ndReg_getNeutronCE => getNeutronCE
  use numPrecision
  use particleDungeon_class,        only : particleDungeon
  use physicalParticle_inter,       only : physicalParticle
  use reactionHandle_inter,         only : reactionHandle
  use RNG_class,                    only : RNG
  use scalarField_inter,            only : getTemperatureFieldPtr, scalarField
  use scatteringKernels_func,       only : asymptoticInelasticScatter, asymptoticScatter, targetVelocity_constXS, &
                                           targetVelocity_DBRCXS
  use tallyAdmin_class,             only : tallyAdmin
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE, uncorrelatedReactionCE_CptrCast
  use universalVariables

  implicit none
  private

  ! Public procedures.
  public :: init

  !!
  !!
  !!
  type, public, abstract, extends(collisionProcessor) :: neutronCECollisionProcessor
    private
    ! Nuclear Data block pointers.
    class(ceNeutronDatabase), pointer :: xsData => null()
    class(ceNeutronMaterial), pointer :: mat => null()
    class(ceNeutronNuclide), pointer :: nuc => null()

    ! Settings.
    real(defReal) :: DBRCeMax = ZERO, DBRCeMin = ZERO, maxE = ZERO, minE = ZERO, threshA = ZERO, threshE = ZERO
  contains
    procedure          :: capture
    procedure          :: getMaximumEnergy
    procedure          :: getMinimumEnergy
    procedure          :: getNuclideIsFissile
    procedure          :: getNuclideMicroXS
    procedure          :: getReaction
    procedure          :: elastic
    procedure          :: inelastic
    procedure          :: init
    procedure          :: sampleCollision
    procedure, private :: scatterFromFixed
    procedure, private :: scatterFromMoving
    procedure, private :: scatterInLAB
  end type neutronCECollisionProcessor

contains
  !!
  !! Process capture reaction
  !!
  subroutine capture(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCECollisionProcessor), intent(inout) :: self
    class(physicalParticle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)                   :: tally
    type(collisionData), intent(inout)                :: collDat
    class(particleDungeon), intent(inout)             :: thisCycle, nextCycle

    call p % setIsDead(.true.)

  end subroutine capture

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
  function getNuclideIsFissile(self) result(isFissile)
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
    type(neutronMicroXSs)                          :: microXS

    call self % nuc % getMicroXSs(E, kT, microXS, rand)

  end subroutine getNuclideMicroXS

  !!
  !!
  !!
  function getReaction(self, reactionChannel, nuclideIdx) result(reaction)
    class(neutronCECollisionProcessor), intent(in) :: self
    integer(shortInt), intent(in)                  :: reactionChannel, nuclideIdx
    class(reactionHandle), pointer                 :: reaction

    reaction => self % xsData % getReaction(reactionChannel, nuclideIdx)

  end function getReaction

  !!
  !! Process elastic scattering
  !!
  !! All CE elastic scattering happens in the CM frame
  !!
  subroutine elastic(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCECollisionProcessor), intent(inout) :: self
    class(physicalParticle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)                   :: tally
    type(collisionData), intent(inout)                :: collDat
    class(particleDungeon), intent(inout)             :: thisCycle, nextCycle
    class(uncorrelatedReactionCE), pointer            :: reac
    logical(defBool)                                  :: isFixed, hasDBRC
    real(defReal)                                     :: energy
    type(CENeutron), pointer                          :: CENeutronPtr
    character(*), parameter                           :: Here = 'elastic (neutronCECollisionProcessor_inter.f90)'

    ! Downcast physicalParticle to CENeutron and retrieve neutron energy.
    CENeutronPtr => castCENeutronPtr(p, .true.)
    energy = CENeutronPtr % getEnergy()
    
    ! Assess if thermal scattering data is needed or not
    if (self % nuc % needsSabEl(energy)) collDat % MT = N_N_ThermEL

    ! Get reaction
    reac => uncorrelatedReactionCE_CptrCast(self % xsData % getReaction(collDat % MT, collDat % nucIdx))
    if (.not. associated(reac)) call fatalError(Here, 'Failed to get elastic neutron scatter.')

    ! Scatter particle
    collDat % A =  self % nuc % getMass()

    ! Retrieve kT from either material or nuclide
    if (.not. self % mat % useTMS(energy)) collDat % kT = self % nuc % getkT()

    ! Check is DBRC is on
    hasDBRC = self % nuc % hasDBRC()
    isFixed = collDat % kT * self % threshE < energy .and. self % threshA < collDat % A .and. .not. hasDBRC

    ! Apply criterion for Free-Gas vs Fixed Target scattering
    if (.not. reac % inCMFrame()) then
      call self % scatterInLAB(CENeutronPtr, collDat, reac)

    elseif (isFixed) then
      call self % scatterFromFixed(CENeutronPtr, collDat, reac)

    else
      call self % scatterFromMoving(CENeutronPtr, collDat, reac)

    end if

  end subroutine elastic

  !!
  !!
  !!
  !!
  !! Process inelastic scattering
  !!
  subroutine inelastic(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCECollisionProcessor), intent(inout) :: self
    class(physicalParticle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)                   :: tally
    type(collisionData), intent(inout)                :: collDat
    class(particleDungeon), intent(inout)             :: thisCycle, nextCycle
    class(uncorrelatedReactionCE), pointer            :: reac
    type(CENeutron), pointer                          :: CENeutronPtr
    character(*), parameter                           :: Here = 'inelastic (neutronCECollisionProcessor_inter.f90)'

    ! Invert inelastic scattering and get reaction
    CENeutronPtr => castCENeutronPtr(p, .true.)
    collDat % MT = self % nuc % invertInelastic(collDat % E, CENeutronPtr % getRNGPtr())
    reac => uncorrelatedReactionCE_CptrCast(self % xsData % getReaction(collDat % MT, collDat % nucIdx))
    if (.not. associated(reac)) call fatalError(Here, 'Failed to retrieve scattering reaction.')

    ! Scatter particle
    if (reac % inCMFrame()) then
      collDat % A =  self % nuc % getMass()
      call self % scatterFromFixed(CENeutronPtr, collDat, reac)

    else
      call self % scatterInLAB(CENeutronPtr, collDat, reac)

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
    character(*), parameter                           :: here = 'init (neutronCECollisionProcessor_inter.f90)'

    ! Initialise superclass.
    call init_super(self, dict)

    ! Get settings from dictionary.
    call dict % getOrDefault(self % minE, 'minEnergy', 1.0E-11_defReal)
    call dict % getOrDefault(self % maxE, 'maxEnergy', 20.0_defReal)
    call dict % getOrDefault(self % threshE, 'energyThreshold', 400.0_defReal)
    call dict % getOrDefault(self % threshA, 'massThreshold', 1.0_defReal)
    call dict % getOrDefault(self % DBRCeMin, 'DBRCeMin', 1.0E-8_defReal)
    call dict % getOrDefault(self % DBRCeMax, 'DBRCeMax', 200.0E-6_defReal)

    ! Verify settings.
    if (self % minE < ZERO) call fatalError(here, 'Minimum energy must be positive.')
    if (self % maxE < ZERO) call fatalError(here, 'Maximum energy must be positive.')
    if (self % maxE <= self % minE) call fatalError(here, 'maxEnergy <= minEnergy.')
    if (self % threshE < ZERO) call fatalError(here,' Energy threshold must be positive.')
    if (self % threshA < ZERO) call fatalError(here, 'Mass threshold must be positive.')

  end subroutine init

  !!
  !! Samples collision without any implicit treatment
  !!
  subroutine sampleCollision(self, p, collDat)
    class(neutronCECollisionProcessor), intent(inout) :: self
    class(physicalParticle), intent(inout)            :: p
    type(collisionData), intent(inout)                :: collDat
    class(scalarField), pointer                       :: temperatureFieldPtr
    real(defReal)                                     :: randomNumber, temperature
    type(CENeutron), pointer                          :: CENeutronPtr
    type(neutronMicroXSs)                             :: microXSs
    type(RNG), pointer                                :: RNGPtr
    character(*), parameter                           :: Here = 'sampleCollision (neutronCECollisionProcessor_inter.f90)'

    ! Verify that particle is CE neutron
    CENeutronPtr => castCENeutronPtr(p, .true.)

    ! Verify and load nuclear data pointer
    self % xsData => ndReg_getNeutronCE()
    if (.not. associated(self % xsData)) call fatalError(Here, 'There is no active Neutron CE data!')

    ! Verify and load material pointer
    self % mat => ceNeutronMaterial_CptrCast(self % xsData % getMaterial(collDat % matIdx))
    if (.not. associated(self % mat)) call fatalError(Here, 'Material is not ceNeutronMaterial')

    ! Retrieve material temperature from temperature field.
    collDat % kT = self % mat % kT
    temperatureFieldPtr => getTemperatureFieldPtr()
    if (associated (temperatureFieldPtr)) then
      temperature = temperatureFieldPtr % at(CENeutronPtr % getCoordsPtr())
      if (ZERO < temperature) collDat % kT = kBoltzmann * temperature / joulesPerMeV

    end if

    ! Select collision nuclide.
    RNGPtr => CENeutronPtr % getRNGPtr()
    call self % mat % sampleNuclide(CENeutronPtr % getEnergy(), collDat % kT, RNGPtr, collDat % nucIdx, collDat % E)

    ! If nuclide was rejected in TMS loop return to tracking
    if (collDat % nucIdx == REJECTED) then
      collDat % MT = noInteraction
      return
      
    end if

    self % nuc => ceNeutronNuclide_CptrCast(self % xsData % getNuclide(collDat % nucIdx))
    if (.not. associated(self % nuc)) call fatalError(Here, 'Failed to retrieve CE Neutron Nuclide')

    ! Select Main reaction channel
    call self % nuc % getMicroXSs(collDat % E, collDat % kT, microXSs, RNGPtr)
    call RNGPtr % generate(randomNumber)
    collDat % MT = microXss % invert(randomNumber)

  end subroutine sampleCollision

  !!
  !! Subroutine to perform scattering from stationary target.
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterFromFixed(self, n_CE, collDat, reac)
    class(neutronCECollisionProcessor), intent(inout) :: self
    type(CENeutron), intent(inout)                    :: n_CE
    type(collisionData), intent(inout)                :: collDat
    class(uncorrelatedReactionCE), intent(in)         :: reac
    real(defReal)                                     :: E_in, E_out, E_outCM, mu, phi
    integer(shortInt)                                 :: MT

    ! Read data
    MT = collDat % MT

    ! Sample mu, phi and outgoing energy
    E_in = n_CE % getEnergy()
    call reac % sampleOut(mu, phi, E_outCM, E_in, n_CE % getRNGPtr())

    ! Save incident energy
    E_out = E_in

    if (MT == N_N_elastic) then
      call asymptoticScatter(E_out, mu, collDat % A)

    else
      call asymptoticInelasticScatter(E_out, mu, E_outCM, collDat % A)

    end if

    ! Update particle state
    call n_CE % rotate(mu, phi)
    call n_CE % setEnergy(E_out)
    collDat % muL = mu

  end subroutine scatterFromFixed

  !!
  !! Subroutine to perform scattering from moving target
  !! Supports only elastic collisions
  !!
  subroutine scatterFromMoving(self, n_CE, collDat, reac)
    class(neutronCECollisionProcessor), intent(inout) :: self
    type(CENeutron), intent(inout)                    :: n_CE
    type(collisionData), intent(inout)                :: collDat
    class(uncorrelatedReactionCE), intent(in)         :: reac
    class(ceNeutronNuclide), pointer                  :: ceNuc0K
    integer(shortInt)                                 :: nucIdx
    logical(defBool)                                  :: inEnergyRange, hasDBRC
    real(defReal)                                     :: A, dummy, energy, kT, maj, mu, phi, speed
    real(defReal), dimension(3)                       :: dir_post, dir_pre, V_cm, V_n, v_t
    type(RNG), pointer                                :: RNGPtr
    character(*), parameter                           :: here = 'scatterFromMoving (neutronCECollisionProcessor_inter.f90)'

    ! Read collision data
    A = collDat % A
    kT = collDat % kT
    nucIdx = collDat % nucIdx

    ! Get neutron direction and velocity
    dir_pre = n_CE % getGlobalDirection()
    energy = n_CE % getEnergy()
    V_n = dir_pre * sqrt(energy)

    ! Sample target velocity with constant XS or with DBRC
    ! Check energy range
    inEnergyRange = energy <= self % DBRCeMax .and. self % DBRCeMin <= energy
    
    ! Check if DBRC is on for this target nuclide
    RNGPtr => n_CE % getRNGPtr()
    hasDBRC = self % nuc % hasDBRC()
    if (inEnergyRange .and. hasDBRC) then
      ! Retrieve 0K nuclide index from DBRC nuclide map
      nucIdx = self % xsData % mapDBRCnuc % get(nucIdx)

      ! Assign pointer for the 0K nuclide
      ceNuc0K => ceNeutronNuclide_CptrCast(self % xsData % getNuclide(nucIdx))
      if (.not. associated(ceNuc0K)) call fatalError(here, 'Failed to retrieve CE neutron nuclide.')

      ! Get elastic scattering 0K majorant
      maj = self % xsData % getScattMicroMajXS(energy, kT, A, nucIdx)

      ! Use DBRC to sample target velocity
      V_t = targetVelocity_DBRCXS(ceNuc0K, energy, dir_pre, A, kT, RNGPtr, maj)

    else
      ! Constant cross section approximation
      V_t = targetVelocity_constXS(energy, dir_pre, A, kT, RNGPtr)

    end if

    ! Calculate Centre-of-Mass velocity
    V_cm = (V_n + V_t * A) / (A + 1)

    ! Move Neutron velocity to CM frame, store speed and calculate new normalised direction
    V_n = V_n - V_cm
    speed = norm2(V_n)
    V_n = V_n / speed

    ! Sample mu and phi in CM frame
    call reac % sampleOut(mu, phi, dummy, energy, RNGPtr)

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
    collDat % muL = dot_product(dir_pre, dir_post)

  end subroutine scatterFromMoving

  !!
  !! Subroutine to perform scattering in LAB frame
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterInLAB(self, n_CE, collDat, reac)
    class(neutronCECollisionProcessor), intent(inout) :: self
    type(CENeutron), intent(inout)                    :: n_CE
    type(collisionData), intent(inout)                :: collDat
    class(uncorrelatedReactionCE), intent(in)         :: reac
    real(defReal)                                     :: E_out, mu, phi ! Azimuthal scatter angle

    ! Sample scattering angles and post-collision energy
    call reac % sampleOut(mu, phi, E_out, n_CE % getEnergy(), n_CE % getRNGPtr())

    ! Update neutron state
    call n_CE % setEnergy(E_out)
    call n_CE % rotate(mu, phi)
    collDat % muL = mu

  end subroutine scatterInLAB

end module neutronCECollisionProcessor_inter