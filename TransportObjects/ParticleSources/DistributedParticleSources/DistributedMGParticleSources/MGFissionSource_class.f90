module MGFissionSource_class

  use distributedMGSource_inter, only : distributedMGSource
  use endfConstants,             only : macroFission
  use errors_mod,                only : fatalError
  use fissionMG_class,           only : fissionMG, fissionMG_TptrCast
  use genericProcedures,         only : numToChar
  use mgNeutronDatabase_inter,   only : mgNeutronDatabase
  use mgNeutronMaterial_inter,   only : mgNeutronMaterial
  use MGParticleState_class,     only : buildMGParticleStatePayload
  use neutronMaterial_inter,     only : neutronMaterial
  use numPrecision
  use RNG_class,                 only : RNG
  use universalVariables,        only : OUTSIDE_MAT, VOID_MAT

  implicit none
  private

  ! Parameters.
  integer(shortInt) :: N_MAXIMUM_ITERATIONS = 10000

  !!
  !!
  !!
  type, public, extends(distributedMGSource) :: MGFissionSource
    private
  contains
    procedure :: getDefaultMaxIterationsNumber
    procedure :: isMaterialInvalid
    procedure :: printInfiniteLoopError
    procedure :: rejectMaterial
    procedure :: sampleFinalPayloadComponents
  end type MGFissionSource

contains
  !!
  !!
  !!
  elemental function getDefaultMaxIterationsNumber(self) result(nMaxIterations)
    class(MGFissionSource), intent(in) :: self
    integer(shortInt)                  :: nMaxIterations

    nMaxIterations = N_MAXIMUM_ITERATIONS

  end function getDefaultMaxIterationsNumber

  !!
  !!
  !!
  elemental function isMaterialInvalid(self, materialIdx) result(isInvalid)
    class(MGFissionSource), intent(in) :: self
    integer(shortInt), intent(in)      :: materialIdx
    logical(defBool)                   :: isInvalid

    isInvalid = materialIdx == OUTSIDE_MAT .or. materialIdx == VOID_MAT

  end function isMaterialInvalid

  !!
  !!
  !!
  subroutine printInfiniteLoopError(self, nMaxIterations)
    class(MGFissionSource), intent(in) :: self
    integer(shortInt), intent(in)      :: nMaxIterations
    character(*), parameter            :: HERE = 'printInfiniteLoopError (MGFissionSource_class.f90)'

    call fatalError(HERE, 'Failed to find a fissile material in: '//numToChar(nMaxIterations)//' attempts.&
                           & Increase the number of maximum attempts or verify that fissile materials are present.')

  end subroutine printInfiniteLoopError

  !!
  !!
  !!
  elemental function rejectMaterial(self, material) result(isRejected)
    class(MGFissionSource), intent(in) :: self
    class(neutronMaterial), intent(in) :: material
    logical(defBool)                   :: isRejected

    isRejected = .not. material % isFissile()

  end function rejectMaterial

  !!
  !!
  !!
  subroutine sampleFinalPayloadComponents(self, MGMaterial, MGDatabase, temperature, rand, payload, mu, phi)
    class(MGFissionSource), intent(in)               :: self
    class(mgNeutronMaterial), intent(in)             :: MGMaterial
    class(mgNeutronDatabase), intent(in)             :: MGDatabase
    real(defReal), intent(in)                        :: temperature
    type(RNG), intent(inout)                         :: rand
    type(buildMGParticleStatePayload), intent(inout) :: payload
    real(defReal), intent(out)                       :: mu, phi
    type(fissionMG), pointer                         :: fissMG
    character(*), parameter                          :: HERE = 'sampleFinalState (MGFissionSource_class.f90)'

    ! Get reaction object
    fissMG => fissionMG_TptrCast(MGDatabase % getReaction(macroFission, payload % materialIdx))
    if (.not. associated(fissMG)) call fatalError(HERE, "Failed to get MG Fission Reaction Object")

    ! Get mu, phi, and energy group.
    call fissMG % sampleOut(mu, phi, payload % energyGroup, self % getEnergyGroup(), rand)

  end subroutine sampleFinalPayloadComponents

end module MGFissionSource_class