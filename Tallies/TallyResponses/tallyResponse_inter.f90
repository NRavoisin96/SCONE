module tallyResponse_inter

  use dictionary_class,        only : dictionary
  use errors_mod,              only : fatalError
  use neutronMaterial_inter,   only : neutronMaterial, neutronMaterial_CptrCast
  use neutronXsPackages_class, only : neutronMacroXSs
  use nuclearDatabase_inter,   only : nuclearDatabase
  use numPrecision
  use physicalParticle_inter,  only : castPhysicalParticlePtr, physicalParticle
  use transportObject_inter,   only : transportObject
  use universalVariables,      only : VOID_MAT

  implicit none
  private


  !!
  !! Abstract interface for all tallyResponses
  !!
  !! Very simple class, which given a particle returns a real number
  !! Real number is used to weight FLUX sample to score reaction rates etc.
  !! Returns only scalar to move all logic for dealing with multiple responces to tallyClerks.
  !! Thus tallyResponses should be quick and easy to write
  !!
  !! Interface:
  !!   init -> Initialise
  !!   get  -> Get velue of the response
  !!   kill -> Return to uninitialised state
  !!
  type, public,abstract :: tallyResponse
    private
  contains
    procedure(get), deferred  :: get
    procedure                 :: getNeutronMacroXS
    procedure(init), deferred :: init
    procedure(kill), deferred :: kill
  end type tallyResponse

  abstract interface

    !!
    !! Initialise Response from dictionary
    !!
    !! Args:
    !!   dict [in] -> DIctionary with the data
    !!
    !! Errors:
    !!   Depend on specific implementation.
    !!   fatalError if there is a mistake in definition
    !!
    subroutine init(self, dict)
      import                              :: dictionary, tallyResponse
      class(tallyResponse), intent(inout) :: self
      class(dictionary), intent(in)       :: dict
    end subroutine init

    !!
    !! Get value of response
    !!
    !! Args:
    !!   p [in]         -> Particle to provide state
    !!   xsData [inout] -> Nuclear Database used by the particle
    !!
    !! Result:
    !!   Value of the response for particle p
    !!
    !! Errors:
    !!   Depend on specific implementation
    !!
    subroutine get(self, object, val, xsData)
      import :: defReal, nuclearDatabase, tallyResponse, transportObject
      class(tallyResponse), intent(in)                :: self
      class(transportObject), intent(in)              :: object
      real(defReal), intent(out)                      :: val
      class(nuclearDatabase), intent(inout), optional :: xsData
    end subroutine get

    !!
    !! Return to uninitialised state
    !!
    !! Args:
    !!   None
    !!
    !! Errors:
    !!   None
    !!
    elemental subroutine kill(self)
      import :: tallyResponse
      class(tallyResponse), intent(inout) :: self
    end subroutine kill

  end interface

contains
  !!
  !!
  !!
  subroutine getNeutronMacroXS(self, object, channel, val, materialIdx, xsData)
    class(tallyResponse), intent(in)                :: self
    class(transportObject), intent(in)              :: object
    integer(shortInt), intent(in)                   :: channel
    real(defReal), intent(out)                      :: val
    integer(shortInt), intent(in), optional         :: materialIdx
    class(nuclearDatabase), intent(inout), optional :: xsData
    class(neutronMaterial), pointer                 :: mat
    class(physicalParticle), pointer                :: p
    integer(shortInt)                               :: particleMaterialIdx, searchMaterialIdx
    type(neutronMacroXSs)                           :: xss
    character(*), parameter                         :: here = 'getNeutronMacroXS (tallyResponse_inter.f90)'

    ! Initialise val = ZERO
    val = ZERO

    ! Downcast transport object to a physical particle. Return immediately if the transport object is not a physical particle.
    p => castPhysicalParticlePtr(object)
    if (.not. associated(p)) return

    ! Get material occupied by the particle. Return if the particle is in the void.
    particleMaterialIdx = p % getMaterialIdx()
    if (particleMaterialIdx == VOID_MAT) return
    searchMaterialIdx = particleMaterialIdx
    if (present(materialIdx)) searchMaterialIdx = materialIdx

    ! Get pointer to active material data. Return if material is not a neutron material.
    if (.not. present(xsData)) call fatalError(here, 'Nuclear database was not provided.')
    mat => neutronMaterial_CptrCast(xsData % getMaterial(searchMaterialIdx))
    if (.not. associated(mat)) return

    ! Retrieve macroscopic cross section and the value corresponding to the specific channel.
    call mat % getMacroXSs(p, xss)
    val = xss % get(channel)

  end subroutine getNeutronMacroXS

end module tallyResponse_inter