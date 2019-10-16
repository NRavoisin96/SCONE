module macroResponse_class

  use numPrecision
  use endfConstants
  use genericProcedures,          only : fatalError, numToChar
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, P_NEUTRON
  use tallyResponse_inter,        only : tallyResponse

  ! Nuclear Data interfaces
  use nuclearDataReg_mod,         only : ndReg_get => get
  use nuclearDatabase_inter,      only : nuclearDatabase
  use ceNeutronMaterial_class,    only : ceNeutronMaterial
  use mgNeutronMaterial_inter,    only : mgNeutronMaterial
  use neutronXsPackages_class,    only : neutronMacroXSs

!  use transportNuclearData_inter, only : transportNuclearData
!  use xsMacroSet_class,           only : xsMacroSet_ptr


  implicit none
  private


  !!
  !! tallyResponse for scoring a single macroscopicXSs
  !!  Currently supports neutrons only
  !!
  !! Sample dictionary input
  !!  name {
  !!     type macroResponse;
  !!     MT   <int>;
  !!  }
  !!
  type, public,extends(tallyResponse) :: macroResponse
    private
    !! Response MT number
    integer(shortInt) :: MT = 0
  contains
    procedure  :: init
    procedure  :: build
    procedure  :: get
    procedure  :: kill
  end type macroResponse

contains

  !!
  !! Initialise Response from dictionary
  !!
  subroutine init(self, dict)
    class(macroResponse), intent(inout) :: self
    class(dictionary), intent(in)       :: dict
    integer(shortInt)                   :: MT
    character(100), parameter :: Here = 'init ( macroResponse_class.f90)'

    ! Load MT number
    call dict % get(MT, 'MT')

    ! Build response
    call self % build(MT)

  end subroutine init

  !!
  !! Build macroResponse from MT number
  !!
  subroutine build(self, MT)
    class(macroResponse), intent(inout) :: self
    integer(shortInt), intent(in)       :: MT
    character(100), parameter :: Here = 'build ( macroResponse_class.f90)'

    ! Check that MT number is valid
    select case(MT)
      case(macroTotal, macroCapture, macroFission, macroNuFission, macroAbsorbtion)
        ! Do nothing. MT is Valid

      case(macroEscatter)
        call fatalError(Here,'Macroscopic Elastic scattering is not implemented yet')

      case default
        call fatalError(Here,'Unrecognised MT number: '// numToChar(self % MT))
    end select

    ! Load MT
    self % MT = MT

  end subroutine build

  !!
  !! Return response value
  !!  Returns 0.0 if the xs type is invalid
  !!
  function get(self, p) result(val)
    class(macroResponse), intent(in) :: self
    class(particle), intent(in)      :: p
    real(defReal)                    :: val
    type(neutronMacroXSs)            :: xss
    class(nuclearDatabase), pointer  :: data

    val = ZERO

    ! Return 0.0 if particle is not neutron
    if(p % type /= P_NEUTRON) return

    ! Get pointer to active data
    data => ndReg_get(p % getType())

    ! Return 0.0 if there is no active data
    if(.not.associated(data)) return

    ! Get XSs
    select type (mat => data % getMaterial(p % matIdx()))
      class is(mgNeutronMaterial)
        call mat % getMacroXSs(xss, p % G, p % pRNG)

      class is(ceNeutronMaterial)
        call mat % getMacroXSs(xss, p % E, p % pRNG)

      class default ! Return if data is uknown
        return

    end select

    val = xss % get(self % MT)

  end function get

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(macroResponse), intent(inout) :: self

    self % MT = 0

  end subroutine kill

end module macroResponse_class
