module fissionPowerResponse_class

  use dictionary_class,        only : dictionary
  use endfConstants,           only : macroFission
  use errors_mod,              only : fatalError
  use nuclearDatabase_inter,   only : nuclearDatabase
  use numPrecision
  use tallyResponse_inter,     only : tallyResponse
  use transportObject_inter,   only : transportObject
  use universalVariables,      only : energyPerFission, joulesPerMeV, ZERO

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
  subroutine get(self, object, val, xsData)
    class(fissionPowerResponse), intent(in)         :: self
    class(transportObject), intent(in)              :: object
    real(defReal), intent(out)                      :: val
    class(nuclearDatabase), intent(inout), optional :: xsData

    call self % getNeutronMacroXS(object, macroFission, val, xsData = xsData)

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