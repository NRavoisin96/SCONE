module densityResponse_class

  use dictionary_class,       only : dictionary
  use nuclearDatabase_inter,  only : nuclearDatabase
  use numPrecision
  use physicalParticle_inter, only : castPhysicalParticlePtr, physicalParticle
  use tallyResponse_inter,    only : tallyResponse
  use transportObject_inter,  only : transportObject

  implicit none
  private

  !!
  !! tallyResponse to score particle density contribution
  !!
  !! Returns the inverse of the particle speed in [cm/s]
  !!
  !! NOTE:
  !!  The speeds are computed from non-relativistic formula for massive particles.
  !!  The small error might appear in MeV range (e.g. for fusion applications)
  !!
  !! Interface:
  !!   tallyResponse Interface
  !!
  type, public, extends(tallyResponse) :: densityResponse
    private
  contains
    procedure :: get
    procedure :: init
    procedure :: kill
  end type densityResponse

contains

  !!
  !! Initialise Response from dictionary
  !!
  !! See tallyResponse_inter for details
  !!
  subroutine init(self, dict)
    class(densityResponse), intent(inout) :: self
    class(dictionary), intent(in)         :: dict

    ! Do nothing

  end subroutine init

  !!
  !! Returns the inverse of the particle speed (response to score particle density)
  !!
  !! See tallyResponse_inter for details
  !!
  subroutine get(self, object, val, xsData)
    class(densityResponse), intent(in)              :: self
    class(transportObject), intent(in)              :: object
    real(defReal), intent(out)                      :: val
    class(nuclearDatabase), intent(inout), optional :: xsData
    class(physicalParticle), pointer                :: p

    ! Gets the particle speed from the particle.
    p => castPhysicalParticlePtr(object, .true.)
    val = ONE / p % getSpeed()

  end subroutine get

  !!
  !! Return to uninitialised State
  !!
  elemental subroutine kill(self)
    class(densityResponse), intent(inout) :: self

    ! Do nothing for nothing can be done

  end subroutine kill

end module densityResponse_class
