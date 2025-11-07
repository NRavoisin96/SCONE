module distributedCESource_inter

  use ceNeutronDatabase_inter,    only : ceNeutronDatabase, ceNeutronDatabase_CptrCast
  use ceNeutronMaterial_class,    only : ceNeutronMaterial, ceNeutronMaterial_CptrCast
  use CEParticleState_class,      only : castCEParticleStatePtr, CEParticleState
  use dictionary_class,           only : dictionary
  use distributedSource_inter,    only : distributedSource, init_super => init, kill_super => kill
  use errors_mod,                 only : fatalError
  use genericProcedures,          only : numToChar
  use geometry_inter,             only : geometry
  use neutronMaterial_inter,      only : neutronMaterial
  use nuclearDatabase_inter,      only : nuclearDatabase
  use numPrecision
  use RNG_class,                  only : RNG
  use transportObjectState_class, only : transportObjectState

  implicit none
  private

  ! Parameters.
  real(defReal), parameter :: DEFAULT_ENERGY = 1.0e-6_defReal

  ! Public procedures.
  public :: init, kill

  !!
  !!
  !!
  type, public, abstract, extends(distributedSource) :: distributedCESource
    private
    real(defReal) :: energy = ZERO
  contains
    procedure                             :: finaliseState
    procedure                             :: getEnergy
    procedure                             :: init
    procedure                             :: kill
    procedure(sampleFinalState), deferred :: sampleFinalState
  end type distributedCESource

  abstract interface
    !!
    !!
    !!
    subroutine sampleFinalState(self, CEMaterial, CEDatabase, temperature, rand, energy, mu, phi)
      import :: ceNeutronDatabase, ceNeutronMaterial, defReal, distributedCESource, RNG
      class(distributedCESource), intent(in) :: self
      class(ceNeutronMaterial), intent(in)   :: CEMaterial
      class(ceNeutronDatabase), intent(in)   :: CEDatabase
      real(defReal), intent(in)              :: temperature
      type(RNG), intent(inout)               :: rand
      real(defReal), intent(out)             :: energy, mu, phi
    end subroutine sampleFinalState

  end interface

contains
  !!
  !!
  !!
  subroutine finaliseState(self, mat, database, temperature, rand, state, mu, phi)
    class(distributedCESource), intent(in)      :: self
    class(neutronMaterial), pointer, intent(in) :: mat
    class(nuclearDatabase), pointer, intent(in) :: database
    real(defReal), intent(in)                   :: temperature
    type(RNG), intent(inout)                    :: rand
    class(transportObjectState), intent(inout)  :: state
    real(defReal), intent(out)                  :: mu, phi
    class(ceNeutronDatabase), pointer           :: CENeutronDatabasePtr
    class(ceNeutronMaterial), pointer           :: CENeutronMaterialPtr
    class(CEParticleState), pointer             :: CEParticleStatePtr
    real(defReal)                               :: energy

    ! Downcast arguments to correct types.
    CENeutronDatabasePtr => ceNeutronDatabase_CptrCast(database)
    CENeutronMaterialPtr => ceNeutronMaterial_CptrCast(mat)
    CEParticleStatePtr => castCEParticleStatePtr(state, .true.)

    ! Sample final state then set energy.
    call self % sampleFinalState(CENeutronMaterialPtr, CENeutronDatabasePtr, temperature, rand, energy, mu, phi)
    call CEParticleStatePtr % setEnergy(energy)

  end subroutine finaliseState

  !!
  !!
  !!
  elemental function getEnergy(self) result(energy)
    class(distributedCESource), intent(in) :: self
    real(defReal)                          :: energy

    energy = self % energy

  end function getEnergy

  !!
  !!
  !!
  subroutine init(self, dict, geom)
    class(distributedCESource), intent(inout) :: self
    class(dictionary), intent(in)             :: dict
    class(geometry), pointer, intent(in)      :: geom

    ! Initialise superclass.
    call init_super(self, dict, geom)

    ! Retrieve energy.
    call dict % getOrDefault(self % energy, 'E', DEFAULT_ENERGY)

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(distributedCESource), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % energy = ZERO

  end subroutine kill

end module distributedCESource_inter