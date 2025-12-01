module CEFissionSource_class

  use ceNeutronDatabase_inter,    only : ceNeutronDatabase
  use ceNeutronMaterial_class,    only : ceNeutronMaterial
  use CEParticleState_class,      only : buildCEParticleStatePayload
  use distributedCESource_inter,  only : distributedCESource
  use endfConstants,              only : N_FISSION
  use errors_mod,                 only : fatalError
  use fissionCE_class,            only : fissionCE, fissionCE_TptrCast
  use genericProcedures,          only : numToChar
  use neutronMaterial_inter,      only : neutronMaterial
  use numPrecision
  use RNG_class,                  only : RNG
  use scalarField_inter,          only : getScalarFieldValue
  use universalVariables,         only : kBoltzmann_MeV, nameDensity, OUTSIDE_MAT, VOID_MAT

  implicit none
  private

  ! Parameters.
  integer(shortInt) :: N_MAXIMUM_ITERATIONS = 10000

  !!
  !!
  !!
  type, public, extends(distributedCESource) :: CEFissionSource
    private
  contains
    procedure :: getDefaultMaxIterationsNumber
    procedure :: isMaterialInvalid
    procedure :: printInfiniteLoopError
    procedure :: rejectMaterial
    procedure :: sampleFinalPayloadComponents
  end type CEFissionSource

contains
  !!
  !!
  !!
  elemental function getDefaultMaxIterationsNumber(self) result(nMaxIterations)
    class(CEFissionSource), intent(in) :: self
    integer(shortInt)                  :: nMaxIterations

    nMaxIterations = N_MAXIMUM_ITERATIONS

  end function getDefaultMaxIterationsNumber

  !!
  !!
  !!
  elemental function isMaterialInvalid(self, materialIdx) result(isInvalid)
    class(CEFissionSource), intent(in) :: self
    integer(shortInt), intent(in)      :: materialIdx
    logical(defBool)                   :: isInvalid

    isInvalid = materialIdx == OUTSIDE_MAT .or. materialIdx == VOID_MAT

  end function isMaterialInvalid

  !!
  !!
  !!
  subroutine printInfiniteLoopError(self, nMaxIterations)
    class(CEFissionSource), intent(in) :: self
    integer(shortInt), intent(in)      :: nMaxIterations
    character(*), parameter            :: HERE = 'printInfiniteLoopError (CEFissionSource_class.f90)'

    call fatalError(HERE, 'Failed to find a fissile material in: '//numToChar(nMaxIterations)//' attempts.&
                           & Increase the number of maximum attempts or verify that fissile materials are present.')

  end subroutine printInfiniteLoopError

  !!
  !!
  !!
  elemental function rejectMaterial(self, material) result(isRejected)
    class(CEFissionSource), intent(in) :: self
    class(neutronMaterial), intent(in) :: material
    logical(defBool)                   :: isRejected

    isRejected = .not. material % isFissile()

  end function rejectMaterial

  !!
  !!
  !!
  subroutine sampleFinalPayloadComponents(self, CEMaterial, CEDatabase, densityFactor, temperature, rand, payload, mu, phi)
    class(CEFissionSource), intent(in)               :: self
    class(ceNeutronMaterial), intent(in)             :: CEMaterial
    class(ceNeutronDatabase), intent(in)             :: CEDatabase
    real(defReal), intent(in)                        :: densityFactor, temperature
    type(RNG), intent(inout)                         :: rand
    type(buildCEParticleStatePayload), intent(inout) :: payload
    real(defReal), intent(out)                       :: mu, phi
    integer(shortInt)                                :: nuclideIdx
    real(defReal)                                    :: E_down, E_out, E_up, kT, sourceEnergy
    type(fissionCE), pointer                         :: fissCE
    character(*), parameter                          :: HERE = 'sampleFinalState (CEFissionSource_class.f90)'

    ! Get energy of the source.
    sourceEnergy = self % getEnergy()

    ! Get energy bounds
    call CEDatabase % energyBounds(E_down, E_up)

    ! Get Nuclide.
    kT = merge(temperature * kBoltzmann_MeV, CEMaterial % kT, ZERO < temperature)
    nuclideIdx = CEMaterial % sampleFission(densityFactor, sourceEnergy, kT, rand)

    ! Get reaction object
    fissCE => fissionCE_TptrCast(CEDatabase % getReaction(N_FISSION, nuclideIdx))
    if (.not. associated(fissCE)) call fatalError(HERE, "Failed to get CE Fission Reaction Object")

    ! Get mu, phi, E_out
    call fissCE % sampleOut(mu, phi, E_out, sourceEnergy, rand)
    payload % energy = min(E_out, E_up)

  end subroutine sampleFinalPayloadComponents

end module CEFissionSource_class