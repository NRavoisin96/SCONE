module neutronMGimp_class

  use collisionProcessor_inter,          only : kill_super => kill
  use dictionary_class,                  only : dictionary
  use neutronMGCollisionProcessor_inter, only : init_super => init, neutronMGCollisionProcessor
  use numPrecision
  use particleDungeon_class,             only : particleDungeon
  use physicalParticle_inter,            only : physicalParticle
  use populationComber_class,            only : populationComber
  use tallyAdmin_class,                  only : tallyAdmin

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
  type, public, extends(neutronMGCollisionProcessor) :: neutronMGimp
    private
    type(populationComber) :: comber
  contains
    procedure :: cutoffs
    procedure :: init
    procedure :: kill
  end type neutronMGimp

contains
  !!
  !! Applay cutoffs or post-collision implicit treatment
  !!
  subroutine cutoffs(self, p, dungeon, tally)
    class(neutronMGimp), intent(in)        :: self
    class(physicalParticle), intent(inout) :: p
    type(particleDungeon), intent(inout)   :: dungeon
    type(tallyAdmin), intent(inout)        :: tally

    call self % comber % cutoffs(p, dungeon, tally)

  end subroutine cutoffs

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self, dict)
    class(neutronMGimp), intent(inout) :: self
    class(dictionary), intent(in)      :: dict

    ! Initialise superclass and population comber.
    call init_super(self, dict)
    call self % comber % init(dict)

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(neutronMGimp), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    call self % comber % kill()

  end subroutine kill

end module neutronMGimp_class