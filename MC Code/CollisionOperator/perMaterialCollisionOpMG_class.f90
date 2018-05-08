module perMaterialCollisionOpMG_class

  use numPrecision
  use endfConstants

  use genericProcedures,              only : fatalError, rotateVector
  use RNG_class,                      only : RNG
  use particle_class,                 only : particle
  use particleDungeon_class,          only : particleDungeon
  use perMaterialNuclearDataMG_inter, only : perMaterialNuclearDataMG
  use collisionOperatorBase_inter,    only : collisionOperatorBase

  ! Cross-section packages to interface with nuclear data
  use xsMacroSet_class,       only : xsMacroSet, xsMacroSet_ptr

  implicit none
  private

  type, public,extends(collisionOperatorBase) :: perMaterialCollisionOpMG
    class(perMaterialNuclearDataMG), pointer :: xsData => null()
    integer(shortInt)                        :: G

  contains
    ! collisionOperatorBase interface
    procedure :: sampleCollision
    procedure :: implicit
    procedure :: scatter
    procedure :: capture
    procedure :: fission
    procedure :: cutoffs
  end type perMaterialCollisionOpMG

contains

  subroutine sampleCollision(self,p,thisCycle,nextCycle)
    class(perMaterialCollisionOpMG), intent(inout) :: self
    class(particle), intent(inout)                 :: p
    class(particleDungeon),intent(inout)           :: thisCycle
    class(particleDungeon),intent(inout)           :: nextCycle
    real(defReal)                                  :: r
    type(xsMacroSet_ptr)                           :: materialXS
    character(100),parameter :: Here =' sampleCollision (perMaterialCollisionOpMG_class.f90)'

    ! Check if particle is multigroup
    if( .not.(p % isMG) ) then
      call fatalError(Here,'CE neutron given to MG operator')
    end if

    ! Retrive Energy
    self % G = p % G

    ! Select Main reaction channel
    call self % xsData % getMatMacroXS(materialXS, self % G, self % matIdx)

    r = self % pRNG % get()

    self % MT = materialXS % invert(r)


  end subroutine sampleCollision

  subroutine implicit(self,p,thisCycle,nextCycle)
    class(perMaterialCollisionOpMG), intent(inout) :: self
    class(particle), intent(inout)                 :: p
    class(particleDungeon),intent(inout)           :: thisCycle
    class(particleDungeon),intent(inout)           :: nextCycle
    type(xsMacroSet_ptr)                       :: materialXS
    type(particle)                             :: pTemp
    real(defReal),dimension(3)                 :: r, dir
    integer(shortInt)                          :: G, G_out, n, i, matIdx
    real(defReal)                              :: nu, wgt, r1, mu, phi
    real(defReal)                              :: sig_fiss, sig_tot, k_eff

    if ( self % xsData % isFissileMat(self % matIdx)) then
      ! Obtain required data
      G      = self % G
      wgt    = p % w
      matIdx = self % matIdx

      nu    = self % xsData % releaseAt(G,G, macroFission, matIdx)
      k_eff = nextCycle % k_eff
      r1    = self % pRNG % get()

      call self % xsData % getMatMacroXS(materialXS, self % G, matIdx)

      sig_fiss = materialXS % fissionXS()
      sig_tot  = materialXS % totalXS()

      r   = p % rGlobal()
      dir = p % dirGlobal()

      ! Sample number of fission sites generated
      n = int(wgt * nu * sig_fiss/(sig_tot*k_eff) + r1, shortInt)

      ! Throw new sites to the next cycle dungeon
      wgt = 1.0
      do i=1,n
        call self % xsData % sampleMuGout(mu, G_out, G, self % pRNG, macroFission, matIdx)
        phi = 2*PI * self % pRNG % get()
        dir = rotateVector(dir,mu,phi)

        call pTemp % build(r,dir,G_out,wgt)

        call nextCycle % throw(pTemp)
      end do
    end if

  end subroutine implicit

  subroutine scatter(self,p,thisCycle,nextCycle)
    class(perMaterialCollisionOpMG), intent(inout) :: self
    class(particle), intent(inout)                 :: p
    class(particleDungeon),intent(inout)           :: thisCycle
    class(particleDungeon),intent(inout)           :: nextCycle
    integer(shortInt)                         :: G_out     ! Post-collision energy group
    real(defReal)                             :: muL       ! Cosine of scattering in LAB frame
    real(defReal)                             :: phi       ! Azimuthal scatter angle
    real(defReal)                             :: w_mul   ! Weight multiplier

    ! Assign MT number
    self % MT = macroAllScatter

    ! Sample Mu and G_out
    call self % xsData % sampleMuGout(muL, G_out, self % G, self % pRNG, self % MT, self % matIdx)
    phi = TWO*PI * self % pRNG % get()

    ! Read scattering multiplicity
    w_mul = self % xsData % releaseAt(self % G, G_out, self % MT, self % matIdx)

    ! Update neutron state
    p % G = G_out
    p % w = p % w * w_mul
    call p % rotate(muL,phi)

  end subroutine scatter

  subroutine capture(self,p,thisCycle,nextCycle)
    class(perMaterialCollisionOpMG), intent(inout) :: self
    class(particle), intent(inout)                 :: p
    class(particleDungeon),intent(inout)           :: thisCycle
    class(particleDungeon),intent(inout)           :: nextCycle

    p % isDead = .true.

  end subroutine capture

  subroutine fission(self,p,thisCycle,nextCycle)
    class(perMaterialCollisionOpMG), intent(inout) :: self
    class(particle), intent(inout)                 :: p
    class(particleDungeon),intent(inout)           :: thisCycle
    class(particleDungeon),intent(inout)           :: nextCycle

    p % isDead = .true.

  end subroutine fission

  subroutine cutoffs(self,p,thisCycle,nextCycle)
    class(perMaterialCollisionOpMG), intent(inout) :: self
    class(particle), intent(inout)                 :: p
    class(particleDungeon),intent(inout)           :: thisCycle
    class(particleDungeon),intent(inout)           :: nextCycle

  end subroutine cutoffs

end module perMaterialCollisionOpMG_class