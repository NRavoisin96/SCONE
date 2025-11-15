module CEMaterialSource_class

  use ceNeutronDatabase_inter,   only : ceNeutronDatabase
  use ceNeutronMaterial_class,   only : ceNeutronMaterial
  use CEParticleState_class,     only : buildCEParticleStatePayload
  use distributedCESource_inter, only : distributedCESource
  use errors_mod,                only : fatalError
  use genericProcedures,         only : numToChar
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
  type, public, extends(distributedCESource) :: CEMaterialSource
    private
    integer(shortInt) :: materialIdx = -1
  contains
    procedure :: getDefaultMaxIterationsNumber
    procedure :: isMaterialInvalid
    procedure :: printInfiniteLoopError
    procedure :: sampleFinalPayloadComponents
  end type CEMaterialSource

contains
  !!
  !!
  !!
  elemental function getDefaultMaxIterationsNumber(self) result(nMaxIterations)
  class(CEMaterialSource), intent(in) :: self
  integer(shortInt)                   :: nMaxIterations

  nMaxIterations = N_MAXIMUM_ITERATIONS

  end function getDefaultMaxIterationsNumber

  !!
  !!
  !!
  elemental function isMaterialInvalid(self, materialIdx) result(isInvalid)
    class(CEMaterialSource), intent(in) :: self
    integer(shortInt), intent(in)       :: materialIdx
    logical(defBool)                    :: isInvalid

    isInvalid = materialIdx == OUTSIDE_MAT .or. materialIdx /= self % materialIdx

  end function isMaterialInvalid

  !!
  !!
  !!
  subroutine printInfiniteLoopError(self, nMaxIterations)
    class(CEMaterialSource), intent(in) :: self
    integer(shortInt), intent(in)       :: nMaxIterations
    character(*), parameter             :: HERE = 'printInfiniteLoopError (CEMaterialSource_class.f90)'

    call fatalError(HERE, 'Failed to find source material in: '//numToChar(nMaxIterations)//' attempts.&
                           & Please check that the volume defined contains the source material.')

  end subroutine printInfiniteLoopError

  !!
  !!
  !!
  subroutine sampleFinalPayloadComponents(self, CEMaterial, CEDatabase, temperature, rand, payload, mu, phi)
    class(CEMaterialSource), intent(in)              :: self
    class(ceNeutronMaterial), intent(in)             :: CEMaterial
    class(ceNeutronDatabase), intent(in)             :: CEDatabase
    real(defReal), intent(in)                        :: temperature
    type(RNG), intent(inout)                         :: rand
    type(buildCEParticleStatePayload), intent(inout) :: payload
    real(defReal), intent(out)                       :: mu, phi

    ! Sample energy and direction.
    payload % energy = self % getEnergy()
    call rand % generateMu(mu)
    call rand % generatePhi(phi)

  end subroutine sampleFinalPayloadComponents

end module CEMaterialSource_class