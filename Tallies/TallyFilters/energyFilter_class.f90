module energyFilter_class

  use CEParticleState_class,      only : CEParticleState
  use dictionary_class,           only : dictionary
  use errors_mod,                 only : fatalError
  use genericProcedures,          only : numToChar
  use MGParticleState_class,      only : MGParticleState
  use numPrecision
  use tallyFilter_inter,          only : tallyFilter
  use transportObjectState_class, only : transportObjectState

  implicit none
  private

  !!
  !! Filter that tests if energy of a state is in a closed, energy (or group) interval.
  !!
  !! Private members:
  !!   Emin -> minimum value of energy
  !!   Emax -> maximum value of energy
  !!   Gtop -> Top energy group
  !!   Glow -> Bottom of allowable energy range
  !!
  !! Interface:
  !!   tallyFilter Interface
  !!   build -> build filter from components
  !!
  !! NOTE: Energy is not enforced to be +ve. Could be useful in debugging
  !!
  !! Sample Dictionary Input:
  !!
  !!   CEfilter {
  !!     type energyFilter;
  !!     Emin 1.0;
  !!     Emax 2.0;
  !!   }
  !!
  !!   MGfilter {
  !!     type energyFilter;
  !!     Gtop 3;
  !!     Glow 17;
  !!  }
  !!
  !!  filter {
  !!     type energyFilter;
  !!     Emin 1.0;
  !!     Emax 2.0;
  !!     Gtop 3;
  !!     Glow 17;
  !! }
  !!
  type, public, extends(tallyFilter) :: energyFilter
    private
    real(defReal)     :: Emax = ZERO, Emin = ZERO
    integer(shortInt) :: Glow = 0, Gtop = 0
  contains
    procedure :: init
    procedure :: isPass

    !! Instance specific procedures
    generic   :: build => build_CE, build_MG, build_CEMG
    procedure :: build_CE
    procedure :: build_MG
    procedure :: build_CEMG

  end type energyFilter

contains

  !!
  !! Initialise energyFilter from dictionary
  !!
  subroutine init(self,dict)
    class(energyFilter), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    real(defReal)                      :: E1, E2
    integer(shortInt)                  :: Gtop, Glow
    logical(defBool)                   :: hasMG, hasCE

    ! Detect which case
    hasMG = dict % isPresent('Glow') .or. dict % isPresent('Gtop')
    hasCE = dict % isPresent('Emin') .or. dict % isPresent('Emax')

    if (hasMG .and. hasCE) then
      ! CE-MG case
      call dict % get(E1, 'Emin')
      call dict % get(E2, 'Emax')
      call dict % get(Gtop, 'Gtop')
      call dict % get(Glow, 'Glow')
      call self % build(E1, E2, Gtop, Glow)

    else if (hasMG) then
      ! MG Case
      call dict % get(Gtop, 'Gtop')
      call dict % get(Glow, 'Glow')
      call self % build(Gtop, Glow)

    else
      ! CE Case
      call dict % get(E1, 'Emin')
      call dict % get(E2, 'Emax')
      call self % build(E1, E2)

    end if

  end subroutine init

  !!
  !! Returns true if energy value is between specified bounds
  !!
  function isPass(self, state) result(passed)
    class(energyFilter), intent(in)         :: self
    class(transportObjectState), intent(in) :: state
    logical(defBool)                        :: passed
    real(defReal)                           :: E
    integer(shortInt)                       :: G
    character(*), parameter                 :: here = 'isPass (energyFilter_class.f90)'

    ! Determine if test is passed depending on the type of particle fed.
    select type(ptr => state)
      class is(CEParticleState)
        E = ptr % getEnergy()
        passed = self % Emin <= E .and. E <= self % Emax

      class is(MGParticleState)
        G = ptr % getEnergyGroup()
        passed = self % Gtop <= G .and. G <= self % Glow

      class default
        call fatalError(here, 'Invalid state type.')
        passed = .false.

    end select

  end function isPass

  !!
  !! Build energyFilter for CE particles only from components
  !!
  !! Args:
  !!   Emin [in] -> minimum energy [MeV]
  !!   Emax [in] -> maximum energy [MeV]
  !!
  !! Errors:
  !!   fatalError if Emin > Emax
  !!
  subroutine build_CE(self, Emin, Emax)
    class(energyFilter), intent(inout) :: self
    real(defReal), intent(in)          :: Emin, Emax
    character(*), parameter            :: here = 'build_CE (energyFilter_class.f90)'

    ! Verify bounds then set values.
    if (Emax <= Emin) &
    call fatalError(here, 'Emin = '//numToChar(Emin)//' is larger than or equal to Emax = '//numToChar(Emax)//'.')
    self % Emin = Emin
    self % Emax = Emax

  end subroutine build_CE

  !!
  !! Build energyFilter for MG particles only from components
  !!
  !! Args:
  !!   Gtop [in] -> maximum energy (lowest index) energy group
  !!   Glow [in] -> minimum energy (highest index) enegy group
  !!
  !! Errors:
  !!   fatalError if Gtop > Elow
  !!
  subroutine build_MG(self, Gtop, Glow)
    class(energyFilter), intent(inout) :: self
    integer(shortInt), intent(in)      :: Gtop, Glow
    character(*), parameter            :: here = 'build_MG (energyFilter_class.f90)'

    ! Verify bounds then set values.
    if (Glow < Gtop) &
    call fatalError(here, 'Gtop = '//numToChar(Gtop)//' is larger than Glow = '//numToChar(Glow)//'.')
    self % Gtop = Gtop
    self % Glow = Glow

  end subroutine build_MG

  !!
  !! Build energyFilter for MG and CE particles from components
  !!
  !! Args:
  !!   Emin [in] -> minimum energy [MeV]
  !!   Emax [in] -> maximum energy [MeV]
  !!   Gtop [in] -> maximum energy (lowest index) energy group
  !!   Glow [in] -> minimum energy (highest index) enegy group
  !!
  !! Errors:
  !!   fatalError if Emin > Emax
  !!   fatalError if Gtop > Elow
  !!
  subroutine build_CEMG(self, Emin, Emax, Gtop, Glow)
    class(energyFilter), intent(inout) :: self
    real(defReal), intent(in)          :: Emin, Emax
    integer(shortInt), intent(in)      :: Gtop, Glow

    call self % build_CE(Emin, Emax)
    call self % build_MG(Gtop, Glow)

  end subroutine build_CEMG

end module energyFilter_class