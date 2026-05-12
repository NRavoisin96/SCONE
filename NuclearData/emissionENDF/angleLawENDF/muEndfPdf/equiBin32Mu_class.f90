module equiBin32Mu_class

  use aceCard_class,      only : aceCard
  use errors_mod,         only : fatalError
  use genericProcedures,  only : linearSearchFloor
  use muEndfPdf_inter,    only : muEndfPdf
  use numPrecision
  use RNG_class,          only : RNG
  use universalVariables, only : valueOutsideArray

  implicit none
  private

  !!
  !! Constructor interface
  !!
  interface equiBin32Mu
    module procedure new_equiBin32Mu
    module procedure new_equiBin32Mu_fromACE
  end interface

  !!
  !! Class that stores PDF of mu in 32 equiprobable bins.
  !! Extends muEndfPdf abstract interface
  !!
  type, public,extends(muEndfPdf) :: equiBin32Mu
    private
    real(defReal),dimension(33) :: boundaries = ZERO
  contains
    ! Superclass procedures
    procedure :: sample
    procedure :: probabilityOf
    procedure :: kill

    ! Local procedures
    procedure :: build
  end type equiBin32Mu

contains

  !!
  !! Samples angle given random number generator
  !!
  function sample(self,rand) result (mu)
    class(equiBin32Mu), intent(in)  :: self
    class(RNG), intent(inout)       :: rand
    real(defReal)                   :: mu
    integer(shortInt)               :: bin
    real(defReal)                   :: f

    ! Sample bin
    bin = floor(32 * rand % get())

    ! Sample angle within bin
    f = rand % get()
    mu = (ONE-f)*self % boundaries(bin) + f* self% boundaries(bin+1)

  end function sample

  !!
  !! Returns probability density of mu
  !!
  function probabilityOf(self,mu) result(prob)
    class(equiBin32Mu), intent(in) :: self
    real(defReal), intent(in)      :: mu
    real(defReal)                  :: prob
    integer(shortInt)              :: idx
    real(defReal)                  :: binWidth
    character(*),parameter         :: HERE = 'probabilityOf (equiBin32Mu_class.f90)'

    ! Find bin location
    idx = linearSearchFloor(self % boundaries, mu)
    if(idx == valueOutsideArray) call fatalError(HERE, 'Requested value is outside array bounds.')

    ! Calculate bin width
    binWidth = self % boundaries(idx + 1) - self % boundaries(idx)

    ! Calculate probability density
    prob = ONE / 32.0_defReal / binWidth

  end function probabilityOf

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(equiBin32Mu), intent(inout) :: self

    self % boundaries = ZERO

  end subroutine kill


  !!
  !! Initialise equiBin32Mu from array of bin boundaries
  !!
  subroutine build(self,boundaries)
    class(equiBin32Mu), intent(inout)        :: self
    real(defReal), dimension(33), intent(in) :: boundaries
    character(100),parameter                 :: Here='init (equiBin32Mu_class.f90)'

    ! Check if the first element of the bin corresponds to -1 and last to 1
    if ( (boundaries(1) /= -1.0_defReal) .or. (boundaries(33) /= 1.0_defReal)) then

       call fatalError(Here, 'Provided bin boundaries do not begin with -1 and end with 1')

    end if

    self % boundaries = boundaries

  end subroutine build

  !!
  !! Constructor from array of bin boundaries
  !!
  function new_equiBin32Mu(boundaries) result(new)
    real(defReal),dimension(33), intent(in) :: boundaries
    type(equiBin32Mu)                       :: new

    ! Allocate space and call initialisation procedure
    call new % build(boundaries)

  end function new_equiBin32Mu

  !!
  !! Constructor form aceCard
  !! aceCard head needs to be set to beginning of data
  !!
  function new_equiBin32Mu_fromACE(ACE) result(new)
    type(aceCard), intent(inout) :: ACE
    type(equiBin32Mu)            :: new
    real(defReal),dimension(33)  :: boundaries

    ! Read Boundaries from ACE library
    boundaries = ACE % readRealArray(33)

    ! Allocate space and call initialisation procedure
    call new % build(boundaries)

  end function new_equiBin32Mu_fromACE

end module equiBin32Mu_class
