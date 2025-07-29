!!
!! Transport operator for hybrid tracking
!!
module transportOperatorHT_class

  use dictionary_class,          only : dictionary
  use errors_mod,                only : fatalError
  use genericProcedures,         only : numToChar
  use geometry_inter,            only : geometry
  use nuclearDatabase_inter,     only : nuclearDatabase
  use nuclearDataReg_mod,        only : ndReg_get => get
  use numPrecision
  use particle_class,            only : particle
  use tallyCodes
  use tallyAdmin_class,          only : tallyAdmin
  use transportOperatorDT_class, only : transportOperatorDT
  use transportOperatorST_class, only : transportOperatorST
  use transportOperator_inter,   only : transportOperator, init_super => init, kill_super => kill
  use universalVariables

  implicit none
  private

  !!
  !! Transport operator that moves a particle with hybrid tracking
  !!
  type, public, extends(transportOperator) :: transportOperatorHT
    private
    type(transportOperatorDT)              :: deltaTracking
    type(transportOperatorST)              :: surfaceTracking
    real(defReal)                          :: cutoff = ZERO  ! Cutoff threshold between ST and DT
  contains
    procedure :: transit => tracking_selection
    ! Override procedure
    procedure :: init
    procedure :: kill

  end type transportOperatorHT

contains

  subroutine tracking_selection(self, p, tally)
    class(transportOperatorHT), intent(inout)              :: self
    class(particle), intent(inout)                         :: p
    type(tallyAdmin), intent(inout)                        :: tally
    real(defReal)                                          :: majorant_inv, sigmaT
    integer(shortInt)                                      :: matIdx

    ! Get majornat XS inverse: 1/Sigma_majorant
    matIdx = p % getMatIdx()
    majorant_inv = ONE / self % xsData % getTrackingXS(p, matIdx, MAJORANT_XS)

    ! Obtain the local cross-section
    sigmaT = self % xsData % getTrackMatXS(p, matIdx)

    ! Cut-off criterion to decide on tracking method
    if (ONE - self % cutoff < sigmaT * majorant_inv) then
      call self % deltaTracking % transit(p, tally)

    else
      call self % surfaceTracking % transit(p, tally)

    end if

  end subroutine tracking_selection

  !!
  !! Initialise HT operator from a dictionary
  !!
  !! See transportOperator_inter for more details
  !!
  subroutine init(self, dict)
    class(transportOperatorHT), intent(inout) :: self
    class(dictionary), intent(in)             :: dict

    ! Initialise superclass
    call init_super(self, dict)

    ! Retrieve DT-ST probability cutoff
    call dict % getOrDefault(self % cutoff, 'cutoff', 0.9_defReal)

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(transportOperatorHT), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % cutoff = ZERO

  end subroutine kill

end module transportOperatorHT_class