module fissionPowerResponse_class

  use dictionary_class,        only : dictionary
  use endfConstants
  use genericProcedures,       only : fatalError
  use neutronMaterial_inter,   only : neutronMaterial, neutronMaterial_CptrCast
  use neutronXsPackages_class, only : neutronMacroXSs
  use nuclearDatabase_inter,   only : nuclearDatabase
  use numPrecision
  use particle_class,          only : particle, P_NEUTRON
  use tallyResponse_inter,     only : tallyResponse
  use universalVariables,      only : energyPerFission, joulesPerMeV, VOID_MAT, ZERO

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(tallyResponse) :: fissionPowerResponse
    private
  contains
    procedure :: get
    procedure :: init
    procedure :: kill
  end type fissionPowerResponse

contains
  !!
  !!
  !!
  subroutine get(self, p, val, xsData)
    class(fissionPowerResponse), intent(in)         :: self
    class(particle), intent(in)                     :: p
    real(defReal), intent(out)                      :: val
    class(nuclearDatabase), intent(inout), optional :: xsData
    integer(shortInt)                               :: matIdx
    type(neutronMacroXSs)                           :: xss
    class(neutronMaterial), pointer                 :: mat
    character(*), parameter                         :: here = 'get (fissionPowerResponse_class.f90)'

    val = ZERO

    ! Return zero if particle is not neutron or if the particle is in void
    if (p % type /= P_NEUTRON) return

    matIdx = p % getMatIdx()
    if (matIdx == VOID_MAT) return

    ! Get pointer to active material data
    if (.not. present(xsData)) call fatalError(here, 'Nuclear database was not provided.')
    mat => neutronMaterial_CptrCast(xsData % getMaterial(matIdx))

    ! Return if material is not a neutronMaterial
    if (.not. associated(mat)) return

    call mat % getMacroXSs(p, xss)
    val = xss % get(macroFission)

    ! Multiply by energy per fission and convert from MeV to J.
    val = val * energyPerFission * joulesPerMeV

  end subroutine get

  !!
  !!
  !!
  subroutine init(self, dict)
    class(fissionPowerResponse), intent(inout) :: self
    class(dictionary), intent(in)              :: dict

    ! Do nothing.

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(fissionPowerResponse), intent(inout) :: self

    ! Local.

  end subroutine kill

end module fissionPowerResponse_class