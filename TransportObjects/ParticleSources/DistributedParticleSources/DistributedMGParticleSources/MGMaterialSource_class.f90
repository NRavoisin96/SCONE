module MGMaterialSource_class

  use distributedMGSource_inter, only : distributedMGSource
  use errors_mod,                only : fatalError
  use genericProcedures,         only : numToChar
  use mgNeutronDatabase_inter,   only : mgNeutronDatabase
  use mgNeutronMaterial_inter,   only : mgNeutronMaterial
  use MGParticleState_class,     only : buildMGParticleStatePayload
  use numPrecision
  use RNG_class,                 only : RNG
  use universalVariables,        only : OUTSIDE_MAT

  implicit none
  private

  ! Parameters.
  integer(shortInt) :: N_MAXIMUM_ITERATIONS = 200

  !!
  !!
  !!
  type, public, extends(distributedMGSource) :: MGMaterialSource
    private
    integer(shortInt) :: materialIdx = -1
  contains
    procedure :: getDefaultMaxIterationsNumber
    procedure :: isMaterialInvalid
    procedure :: printInfiniteLoopError
    procedure :: sampleFinalPayloadComponents
  end type MGMaterialSource

contains
  !!
  !!
  !!
  elemental function getDefaultMaxIterationsNumber(self) result(nMaxIterations)
  class(MGMaterialSource), intent(in) :: self
  integer(shortInt)                   :: nMaxIterations

  nMaxIterations = N_MAXIMUM_ITERATIONS

  end function getDefaultMaxIterationsNumber

  !!
  !!
  !!
  elemental function isMaterialInvalid(self, materialIdx) result(isInvalid)
    class(MGMaterialSource), intent(in) :: self
    integer(shortInt), intent(in)       :: materialIdx
    logical(defBool)                    :: isInvalid

    isInvalid = materialIdx == OUTSIDE_MAT .or. materialIdx /= self % materialIdx

  end function isMaterialInvalid

  !!
  !!
  !!
  subroutine printInfiniteLoopError(self, nMaxIterations)
    class(MGMaterialSource), intent(in) :: self
    integer(shortInt), intent(in)       :: nMaxIterations
    character(*), parameter             :: HERE = 'printInfiniteLoopError (MGMaterialSource_class.f90)'

    call fatalError(HERE, 'Failed to find source material in: '//numToChar(nMaxIterations)//' attempts.&
                           & Please check that the volume defined contains the source material.')

  end subroutine printInfiniteLoopError

  !!
  !!
  !!
  subroutine sampleFinalPayloadComponents(self, MGMaterial, MGDatabase, temperature, rand, payload, mu, phi)
    class(MGMaterialSource), intent(in)              :: self
    class(mgNeutronMaterial), intent(in)             :: MGMaterial
    class(mgNeutronDatabase), intent(in)             :: MGDatabase
    real(defReal), intent(in)                        :: temperature
    type(RNG), intent(inout)                         :: rand
    type(buildMGParticleStatePayload), intent(inout) :: payload
    real(defReal), intent(out)                       :: mu, phi

    ! Sample energy and direction.
    payload % energyGroup = self % getEnergyGroup()
    call rand % generateMu(mu)
    call rand % generatePhi(phi)

  end subroutine sampleFinalPayloadComponents

end module MGMaterialSource_class