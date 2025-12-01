module distributedCESource_inter

  use ceNeutronDatabase_inter,    only : ceNeutronDatabase, ceNeutronDatabase_CptrCast
  use ceNeutronMaterial_class,    only : ceNeutronMaterial, ceNeutronMaterial_CptrCast
  use CENeutronState_class,       only : CENeutronState
  use CEParticleState_class,      only : buildCEParticleStatePayload, castBuildCEParticleStatePayloadPtr
  use dictionary_class,           only : dictionary
  use distributedSource_inter,    only : distributedSource, init_super => init, kill_super => kill
  use errors_mod,                 only : fatalError
  use genericProcedures,          only : numToChar
  use geometry_inter,             only : geometry
  use neutronMaterial_inter,      only : neutronMaterial
  use nuclearDatabase_inter,      only : nuclearDatabase
  use numPrecision
  use RNG_class,                  only : RNG
  use transportObjectState_class, only : buildTransportObjectStatePayload, transportObjectState

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
    procedure                                         :: allocatePayloadAndState
    procedure                                         :: finalisePayload
    procedure                                         :: getEnergy
    procedure                                         :: init
    procedure                                         :: kill
    procedure(sampleFinalPayloadComponents), deferred :: sampleFinalPayloadComponents
  end type distributedCESource

  abstract interface
    !!
    !!
    !!
    subroutine sampleFinalPayloadComponents(self, CEMaterial, CEDatabase, densityFactor, temperature, rand, payload, mu, phi)
      import :: buildCEParticleStatePayload, ceNeutronDatabase, ceNeutronMaterial, defReal, distributedCESource, RNG
      class(distributedCESource), intent(in)           :: self
      class(ceNeutronMaterial), intent(in)             :: CEMaterial
      class(ceNeutronDatabase), intent(in)             :: CEDatabase
      real(defReal), intent(in)                        :: densityFactor, temperature
      type(RNG), intent(inout)                         :: rand
      type(buildCEParticleStatePayload), intent(inout) :: payload
      real(defReal), intent(out)                       :: mu, phi
    end subroutine sampleFinalPayloadComponents

  end interface

contains
  !!
  !!
  !!
  subroutine allocatePayloadAndState(self, payload, state)
    class(distributedCESource), intent(in)                            :: self
    class(buildTransportObjectStatePayload), allocatable, intent(out) :: payload
    class(transportObjectState), allocatable, intent(out)             :: state

    ! Allocate payload and state. For now only allocate CENeutronState. 
    ! Can be extended in the future to include more CE particle types.
    allocate(buildCEParticleStatePayload :: payload)
    allocate(CENeutronState :: state)

  end subroutine allocatePayloadAndState

  !!
  !!
  !!
  subroutine finalisePayload(self, mat, database, densityFactor, temperature, rand, payload, mu, phi)
    class(distributedCESource), intent(in)                 :: self
    class(neutronMaterial), pointer, intent(in)            :: mat
    class(nuclearDatabase), pointer, intent(in)            :: database
    real(defReal), intent(in)                              :: densityFactor, temperature
    type(RNG), intent(inout)                               :: rand
    class(buildTransportObjectStatePayload), intent(inout) :: payload
    real(defReal), intent(out)                             :: mu, phi
    class(ceNeutronDatabase), pointer                      :: CENeutronDatabasePtr
    class(ceNeutronMaterial), pointer                      :: CENeutronMaterialPtr
    type(buildCEParticleStatePayload), pointer             :: CEPayloadPtr

    ! Downcast arguments to correct types.
    CENeutronDatabasePtr => ceNeutronDatabase_CptrCast(database)
    CENeutronMaterialPtr => ceNeutronMaterial_CptrCast(mat)
    CEPayloadPtr => castBuildCEParticleStatePayloadPtr(payload, .true.)

    ! Sample final state components.
    call self % sampleFinalPayloadComponents(CENeutronMaterialPtr, CENeutronDatabasePtr, densityFactor, temperature, rand, &
                                             CEPayloadPtr, mu, phi)

  end subroutine finalisePayload

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