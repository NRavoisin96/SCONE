module weightResponse_class

  use dictionary_class,       only : dictionary
  use endfConstants
  use errors_mod,             only : fatalError
  use genericProcedures,      only : numToChar
  use numPrecision
  use nuclearDatabase_inter,  only : nuclearDatabase
  use neutronMaterial_inter,  only : neutronMaterial, neutronMaterial_CptrCast
  use physicalParticle_inter, only : castPhysicalParticlePtr, physicalParticle
  use tallyResponse_inter,    only : tallyResponse
  use transportObject_inter,  only : transportObject

  implicit none
  private

  !!
  !! tallyResponse for scoring particle weights
  !!  Currently supports neutrons only
  !!
  !! Interface:
  !!   tallyResponse interface
  !!
  !! Sample dictionary input
  !!  name {
  !!     type weightResponse; moment 1;
  !!  }
  !!
  type, public, extends(tallyResponse) :: weightResponse
    private
    integer(shortInt) :: moment = 0
  contains
    ! Superclass Procedures
    procedure  :: init
    procedure  :: get
    procedure  :: kill

  end type weightResponse

contains

  !!
  !! Initialise Response from dictionary
  !!
  !! See tallyResponse_inter for details
  !!
  subroutine init(self, dict)
    class(weightResponse), intent(inout) :: self
    class(dictionary), intent(in)        :: dict
    character(*), parameter              :: here = 'init (weightResponse_class.f90)'

    ! Get response moment to be calculated
    call dict % getOrDefault(self % moment, 'moment', 1)
    if (self % moment < 0) call fatalError(here, 'Moment must be larger than or equal to 0.')

  end subroutine init

  !!
  !! Return response value
  !!
  !! See tallyResponse_inter for details
  !!
  !! Errors:
  !!   Return ZERO if particle is not a Neutron
  !!
  subroutine get(self, object, val, xsData)
    class(weightResponse), intent(in)               :: self
    class(transportObject), intent(in)              :: object
    real(defReal), intent(out)                      :: val
    class(nuclearDatabase), intent(inout), optional :: xsData
    class(neutronMaterial), pointer                 :: mat
    class(physicalParticle), pointer                :: p
    integer(shortInt)                               :: matIdx
    real(defReal)                                   :: factor, weight
    character(*), parameter                         :: here = 'get (weightResponse_class.f90)'

    val = ZERO

    ! Return if particle is not neutron.
    p => castPhysicalParticlePtr(object)
    if (.not. associated(p)) return

    ! Get pointer to active material data
    matIdx = p % getMaterialIdx()
    if (.not. present(xsData)) call fatalError(here, 'Nuclear database was not provided.')
    mat => neutronMaterial_CptrCast(xsData % getMaterial(matIdx))

    ! Return if material is not a neutronMaterial
    if (.not.associated(mat)) return

    weight = p % getWeight()
    if (self % moment == 0) then
      factor = ONE / weight

    else
      factor = weight * (self % moment - 1)

    end if
    val = xsData % getTotalMatXS(p, matIdx) * factor

  end subroutine get

  !!
  !! Return to uninitialised State
  !!
  elemental subroutine kill(self)
    class(weightResponse), intent(inout) :: self

    self % moment = 0

  end subroutine kill

end module weightResponse_class
