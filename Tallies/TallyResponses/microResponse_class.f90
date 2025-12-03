module microResponse_class

  use dictionary_class,        only : dictionary
  use endfConstants
  use errors_mod,              only : fatalError
  use genericProcedures,       only : numToChar
  use materialMenu_mod,        only : materialItem, matName, nMat, getMatPtr
  use nuclearDatabase_inter,   only : nuclearDatabase
  use numPrecision
  use tallyResponse_inter,     only : tallyResponse
  use transportObject_inter,   only : transportObject

  implicit none
  private

  !!
  !! tallyResponse for scoring a single microscopicXSs
  !!   Currently supports neutrons only
  !!
  !! Private Members:
  !!   MT     -> MT number of the microscopic reaction for weighting
  !!   matIdx -> index of the material that contains (only) the nuclide wanted
  !!   dens   -> atomic density of the nuclide
  !!
  !! Interface:
  !!   tallyResponse interface
  !!   build -> Initialise directly from MT number
  !!
  !! Sample dictionary input:
  !!  name {
  !!     type microResponse;
  !!     MT   <int>;
  !!     material <matName>;
  !!  }
  !!
  !! Note:
  !!   The material <matName> must include only one nuclide. The final estimate
  !!   is independent of the nuclide atomic density, which can be any value but zero.
  !!
  type, public, extends(tallyResponse) :: microResponse
    private
    !! Response MT number
    integer(shortInt) :: matIdx = 0, MT = 0
    real(defReal)     :: dens = ZERO
  contains
    ! Superclass Procedures
    procedure  :: init
    procedure  :: get
    procedure  :: kill

    ! Local Procedures
    procedure  :: build

  end type microResponse

contains

  !!
  !! Initialise Response from dictionary
  !!
  !! See tallyResponse_inter for details
  !!
  !! Errors:
  !!   fatalError if the material contains more than one nuclide
  !!   fatalError if the nuclide has density 0.0
  !!
  subroutine init(self, dict)
    class(microResponse), intent(inout)      :: self
    class(dictionary), intent(in)            :: dict
    character(15)                            :: mName
    integer(shortInt)                        :: MT, i
    real(defReal), dimension(:), allocatable :: atomicDensities
    type(materialItem), pointer              :: mat
    character(*), parameter                  :: HERE = 'init ( microResponse_class.f90)'

    ! Load MT number and material name
    call dict % get(MT, 'MT')
    call dict % get(mName, 'material')

    ! Find corresponding material index
    do i = 1, nMat()
      if (mName == matName(i)) self % matIdx = i

    end do

    ! Get pointer to the material
    mat => getMatPtr(self % matIdx)

    atomicDensities = mat % getAtomicDensities()
    if (1 < size(atomicDensities)) call fatalError(HERE, 'Material: '//trim(mName)//' has more than one nuclide.')
    self % dens = atomicDensities(1)
    if (self % dens == ZERO) call fatalError(HERE, 'Density of material: '//trim(mName)//' cannot be ZERO.')
    
    ! Build response.
    call self % build(MT)

  end subroutine init

  !!
  !! Build microResponse from MT number
  !!
  !! Args:
  !!   MT [in] -> MT number for weighting
  !!
  !! Errors:
  !!   fatalError if MT is invalid
  !!
  subroutine build(self, MT)
    class(microResponse), intent(inout) :: self
    integer(shortInt), intent(in)       :: MT
    character(*), parameter             :: Here = 'build (microResponse_class.f90)'

    ! Check that MT number is valid and load MT
    select case(MT)
      case(N_ABSORPTION)
        self % MT = macroAbsorption

      case(N_FISSION)
        self % MT = macroFission

      case(N_GAMMA)
        self % MT = macroCapture

      case(N_heating)
        self % MT = macroHeating

      case(N_N_ELASTIC)
        self % MT = macroEScatter

      case(N_TOTAL)
        self % MT = macroTotal

      case default
        call fatalError(Here,'Unrecognised MT number: '//numToChar(MT)//'.')

    end select

  end subroutine build

  !!
  !! Return response value
  !!
  !! See tallyResponse_inter for details
  !!
  !! Errors:
  !!   Return ZERO if particle is not a Neutron
  !!
  subroutine get(self, object, val, xsData)
    class(microResponse), intent(in)                :: self
    class(transportObject), intent(in)              :: object
    real(defReal), intent(out)                      :: val
    class(nuclearDatabase), intent(inout), optional :: xsData

    call self % getNeutronMacroXS(object, self % MT, val, materialIdx = self % matIdx, xsData = xsData)

    ! Normalise the macroscopic cross section with the atomic density
    val = val / self % dens

  end subroutine get

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(microResponse), intent(inout) :: self

    self % matIdx = 0
    self % MT = 0
    self % dens = ZERO

  end subroutine kill

end module microResponse_class
