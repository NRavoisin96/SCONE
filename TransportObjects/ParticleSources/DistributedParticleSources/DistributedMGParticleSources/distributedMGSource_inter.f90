module distributedMGSource_inter

  use dictionary_class,           only : dictionary
  use distributedSource_inter,    only : distributedSource, init_super => init, kill_super => kill
  use geometry_inter,             only : geometry
  use mgNeutronDatabase_inter,    only : mgNeutronDatabase, mgNeutronDatabase_CptrCast
  use mgNeutronMaterial_inter,    only : mgNeutronMaterial, mgNeutronMaterial_CptrCast
  use MGNeutronState_class,       only : MGNeutronState
  use MGParticleState_class,      only : buildMGParticleStatePayload, castBuildMGParticleStatePayloadPtr
  use neutronMaterial_inter,      only : neutronMaterial
  use nuclearDatabase_inter,      only : nuclearDatabase
  use numPrecision
  use RNG_class,                  only : RNG
  use transportObjectState_class, only : buildTransportObjectStatePayload, transportObjectState

  implicit none
  private

  ! Parameters.
  integer(shortInt), parameter :: DEFAULT_ENERGY_GROUP = 1

  !!
  !!
  !!
  type, public, abstract, extends(distributedSource) :: distributedMGSource
    private
    integer(shortInt) :: energyGroup = 0
  contains
    procedure                                         :: allocatePayloadAndState
    procedure                                         :: finalisePayload
    procedure                                         :: getEnergyGroup
    procedure                                         :: init
    procedure                                         :: kill
    procedure(sampleFinalPayloadComponents), deferred :: sampleFinalPayloadComponents
  end type distributedMGSource

  abstract interface
    !!
    !!
    !!
    subroutine sampleFinalPayloadComponents(self, MGMaterial, MGDatabase, temperature, rand, payload, mu, phi)
      import :: buildMGParticleStatePayload, mgNeutronDatabase, mgNeutronMaterial, defReal, distributedMGSource, RNG
      class(distributedMGSource), intent(in)           :: self
      class(mgNeutronMaterial), intent(in)             :: MGMaterial
      class(mgNeutronDatabase), intent(in)             :: MGDatabase
      real(defReal), intent(in)                        :: temperature
      type(RNG), intent(inout)                         :: rand
      type(buildMGParticleStatePayload), intent(inout) :: payload
      real(defReal), intent(out)                       :: mu, phi
    end subroutine sampleFinalPayloadComponents

  end interface

contains
  !!
  !!
  !!
  subroutine allocatePayloadAndState(self, payload, state)
    class(distributedMGSource), intent(in)                            :: self
    class(buildTransportObjectStatePayload), allocatable, intent(out) :: payload
    class(transportObjectState), allocatable, intent(out)             :: state

    ! Allocate payload and state. For now only allocate CENeutronState. 
    ! Can be extended in the future to include more CE particle types.
    allocate(buildMGParticleStatePayload :: payload)
    allocate(MGNeutronState :: state)

  end subroutine allocatePayloadAndState

  !!
  !!
  !!
  subroutine finalisePayload(self, mat, database, temperature, rand, payload, mu, phi)
    class(distributedMGSource), intent(in)                 :: self
    class(neutronMaterial), pointer, intent(in)            :: mat
    class(nuclearDatabase), pointer, intent(in)            :: database
    real(defReal), intent(in)                              :: temperature
    type(RNG), intent(inout)                               :: rand
    class(buildTransportObjectStatePayload), intent(inout) :: payload
    real(defReal), intent(out)                             :: mu, phi
    class(mgNeutronDatabase), pointer                      :: MGNeutronDatabasePtr
    class(mgNeutronMaterial), pointer                      :: MGNeutronMaterialPtr
    type(buildMGParticleStatePayload), pointer             :: MGPayloadPtr

    ! Downcast arguments to correct types.
    MGNeutronDatabasePtr => mgNeutronDatabase_CptrCast(database)
    MGNeutronMaterialPtr => mgNeutronMaterial_CptrCast(mat)
    MGPayloadPtr => castBuildMGParticleStatePayloadPtr(payload)

    ! Sample final state components.
    call self % sampleFinalPayloadComponents(MGNeutronMaterialPtr, MGNeutronDatabasePtr, temperature, rand, &
                                             MGPayloadPtr, mu, phi)

  end subroutine finalisePayload

  !!
  !!
  !!
  elemental function getEnergyGroup(self) result(energyGroup)
    class(distributedMGSource), intent(in) :: self
    integer(shortInt)                      :: energyGroup

    energyGroup = self % energyGroup

  end function getEnergyGroup

  !!
  !!
  !!
  subroutine init(self, dict, geom)
    class(distributedMGSource), intent(inout) :: self
    class(dictionary), intent(in)             :: dict
    class(geometry), pointer, intent(in)      :: geom

    ! Initialise superclass.
    call init_super(self, dict, geom)

    ! Retrieve energy group.
    call dict % getOrDefault(self % energyGroup, 'G', DEFAULT_ENERGY_GROUP)

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(distributedMGSource), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % energyGroup = 0

  end subroutine kill

end module distributedMGSource_inter