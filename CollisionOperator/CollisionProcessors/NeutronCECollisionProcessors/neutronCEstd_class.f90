module neutronCEstd_class

  use numPrecision
  use endfConstants
  use genericProcedures,                 only : fatalError, rotateVector, numToChar
  use dictionary_class,                  only : dictionary
  use RNG_class,                         only : RNG

  ! Particle types
  use particle_class,                    only : particle, particleState
  use particleDungeon_class,             only : particleDungeon

  ! Abstarct interface
  use collisionProcessor_inter,          only : collisionData
  use neutronCECollisionProcessor_inter, only : neutronCECollisionProcessor, init_super => init

  ! Nuclear Data Interfaces

  ! Nuclear reactions
  use fissionCE_class,                   only : fissionCE, fissionCE_TptrCast

  ! Cross-Section Packages
  use neutronXsPackages_class,           only : neutronMicroXSs

  ! Tally interfaces
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
    class(neutronCEstd), intent(inout)    :: self
    class(particle), intent(inout)        :: p
    type(tallyAdmin), intent(inout)       :: tally
    type(collisionData), intent(inout)    :: collDat
    class(particleDungeon), intent(inout) :: thisCycle, nextCycle
    type(fissionCE), pointer              :: fission
    type(neutronMicroXSs)                 :: microXSs
    type(particleState)                   :: pTemp
    real(defReal), dimension(3)           :: r, dir
    integer(shortInt)                     :: n, i
    real(defReal)                         :: wgt, w0, randomNumber, E_out, mu, phi, &
                                             sig_nufiss, sig_tot, k_eff
    character(*), parameter               :: Here = 'implicit (neutronCEstd_class.f90)'

    ! Generate fission sites if nuclide is fissile
    if (self % getNuclideIsFissile()) then
      ! Obtain required data
      wgt = p % w                ! Current weight
      w0 = p % preHistory % wgt ! Starting weight
      k_eff = p % k_eff            ! k_eff for normalisation
      call p % pRNG % generate(randomNumber)     ! Random number to sample sites

      ! Retrieve cross section at the energy used for reaction sampling
      call self % getNuclideMicroXS(collDat % E, collDat % kT, p % pRNG, microXSs)

      sig_nufiss = microXSs % nuFission
      sig_tot = microXSs % total

      ! Sample number of fission sites generated
      ! Support -ve weight particles
      n = int(abs((wgt * sig_nufiss) / (w0 * sig_tot * k_eff)) + randomNumber, shortInt)

      ! Shortcut particle generation if no particles were sampled
      if (n < 1) return

      ! Get fission Reaction
      fission => fissionCE_TptrCast(self % getReaction(N_FISSION, collDat % nucIdx))
      if (.not.associated(fission)) call fatalError(Here, "Failed to get fissionCE")

      ! Store new sites in the next cycle dungeon
      wgt =  sign(w0, wgt)
      r = p % rGlobal()

      do i = 1, n
        call fission % sampleOut(mu, phi, E_out, p % E, p % pRNG)
        dir = rotateVector(p % dirGlobal(), mu, phi)
        E_out = min(E_out, self % getMaximumEnergy())

        ! Copy extra detail from parent particle (i.e. time, flags ect.)
        pTemp = p

        ! Overwrite position, direction, energy and weight
        pTemp % r = r
        pTemp % dir = dir
        pTemp % E = E_out
        pTemp % wgt = wgt
        pTemp % collisionN = 0

        call nextCycle % detain(pTemp)

        ! Report birth of new particle
        call tally % reportSpawn(N_FISSION, p, pTemp)

      end do
      
    end if

  end subroutine implicit

  !!
  !! Process fission reaction
  !!
  subroutine fission(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(tallyAdmin), intent(inout)      :: tally
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon), intent(inout) :: thisCycle
    class(particleDungeon), intent(inout) :: nextCycle

    p % isDead =.true.

  end subroutine fission

  !!
  !! Apply cutoffs
  !!
  subroutine cutoffs(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)    :: self
    class(particle), intent(inout)        :: p
    type(tallyAdmin), intent(inout)       :: tally
    type(collisionData), intent(inout)    :: collDat
    class(particleDungeon), intent(inout) :: thisCycle, nextCycle

    if (p % E < self % getMinimumEnergy()) p % isDead = .true.

  end subroutine cutoffs

end module neutronCEstd_class