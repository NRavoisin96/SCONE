module nuclideInfo_class

  use dictionary_class,  only : dictionary
  use errors_mod,        only : fatalError
  use genericProcedures, only : charToInt
  use numPrecision

  implicit none
  private

  !!
  !!
  !!
  type, public :: buildNuclideInfoPayload
    character(nameLen), dimension(:), allocatable :: sabFiles
    integer(shortInt)                             :: atomicNumber = 0, evaluationNumber = 0, massNumber = 0
    logical(defBool)                              :: hasSab = .false., sabMix = .false.
  end type buildNuclideInfoPayload

  !!
  !! Information about a single nuclide
  !!
  !! Based somewhat on MCNP conventions.
  !! Atomic and Mass number identify clearly a nuclide species
  !! Evaluation number T allows to refer to multiple states/evaluations of the same nuclide species
  !! E.G. at a Different temperature as in MCNP Library.
  !!
  !! Public members:
  !!   Z -> Atomic number
  !!   A -> Mass number
  !!   T -> Evaluation number
  !!   hasSab -> Does the nuclide have S(a,b) data?
  !!   sabMix -> Does the nuclide mix S(a,b) data?
  !!   file_Sab1 -> First (and maybe only) S(a,b) file
  !!   file_Sab2 -> Second S(a,b) file
  !!
  !! Interface:
  !!   init -> build from a string
  !!
  type, public :: nuclideInfo
    private
    character(nameLen), dimension(:), allocatable :: sabFiles
    integer(shortInt)                             :: A = 0, T = 0, Z = 0
    logical(defBool)                              :: hasSab = .false., sabMix = .false.
    real(defReal)                                 :: density = ZERO ! Atomic density in b⁻¹ cm⁻¹.
  contains
    procedure :: display
    procedure :: getAtomicNumber
    procedure :: getDensity
    procedure :: getEvaluationNumber
    procedure :: getHasSab
    procedure :: getMassNumber
    procedure :: getSabFiles
    procedure :: getSabMix
    procedure :: init
    procedure :: kill
    procedure :: setDensity
    procedure :: toChar
  end type nuclideInfo

contains
  !!
  !!
  !!
  subroutine display(self, format)
    class(nuclideInfo), intent(in) :: self
    character(*), intent(in)       :: format

    print format, self % Z, self % A, self % T, self % density

  end subroutine display

  !!
  !!
  !!
  elemental function getAtomicNumber(self) result(atomicNumber)
    class(nuclideInfo), intent(in) :: self
    integer(shortInt)              :: atomicNumber

    atomicNumber = self % Z

  end function getAtomicNumber

  !!
  !!
  !!
  elemental function getDensity(self) result(density)
    class(nuclideInfo), intent(in) :: self
    real(defReal)                  :: density

    density = self % density

  end function getDensity

  !!
  !!
  !!
  elemental function getEvaluationNumber(self) result(evaluationNumber)
    class(nuclideInfo), intent(in) :: self
    integer(shortInt)              :: evaluationNumber

    evaluationNumber = self % T

  end function getEvaluationNumber

  !!
  !!
  !!
  elemental function getHasSab(self) result(hasSab)
    class(nuclideInfo), intent(in) :: self
    logical(defBool)               :: hasSab

    hasSab = self % hasSab

  end function getHasSab

  !!
  !!
  !!
  elemental function getMassNumber(self) result(massNumber)
    class(nuclideInfo), intent(in) :: self
    integer(shortInt)              :: massNumber

    massNumber = self % A

  end function getMassNumber

  !!
  !!
  !!
  pure function getSabFiles(self) result(sabFiles)
    class(nuclideInfo), intent(in)                :: self
    character(nameLen), dimension(:), allocatable :: sabFiles

    sabFiles = self % sabFiles

  end function getSabFiles

  !!
  !!
  !!
  elemental function getSabMix(self) result(sabMix)
    class(nuclideInfo), intent(in) :: self
    logical(defBool)               :: sabMix

    sabMix = self % sabMix

  end function getSabMix

  !!
  !!
  !!
  elemental subroutine init(self, payload)
    class(nuclideInfo), intent(inout)         :: self
    type(buildNuclideInfoPayload), intent(in) :: payload

    ! Copy everything from payload.
    if (allocated(payload % sabFiles)) self % sabFiles = payload % sabFiles
    self % A = payload % massNumber
    self % T = payload % evaluationNumber
    self % Z = payload % atomicNumber
    self % hasSab = payload % hasSab
    self % sabMix = payload % sabMix

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(nuclideInfo), intent(inout) :: self

    if (allocated(self % sabFiles)) deallocate(self % sabFiles)
    self % A = 0
    self % T = 0
    self % Z = 0
    self % hasSab = .false.
    self % sabMix = .false.
    self % density = ZERO

  end subroutine kill

  !!
  !!
  !!
  elemental subroutine setDensity(self, density)
    class(nuclideInfo), intent(inout) :: self
    real(defReal), intent(in)         :: density

    self % density = density

  end subroutine setDensity

  !!
  !! Convert nuclide information to the definition character
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Character in format ZZAAA.TT that describes nuclide definition
  !!
  !! Errors:
  !!   None
  !!
  elemental function toChar(self) result(str)
    class(nuclideInfo), intent(in) :: self
    character(nameLen)             :: str
    character(3)                   :: ZZ
    character(3)                   :: AAA
    character(2)                   :: TT

    write(ZZ, '(I3)') self % Z
    write(AAA, '(I3.3)') self % A
    write(TT, '(I2.2)') self % T

    str = trim(adjustl(ZZ)) // AAA // "." // TT

  end function toChar

end module nuclideInfo_class