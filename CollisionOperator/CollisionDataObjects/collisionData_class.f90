module collisionData_class

  use endfConstants,        only : noInteraction
  use numPrecision
  use RNG_class,            only : RNG

  implicit none
  private

  !!
  !! Data package with all relevant data about the collision to move beetween customisable
  !! procedures
  !!
  type, public :: collisionData
    integer(shortInt)           :: matIdx = -1, MT = noInteraction, n = 0, nucIdx = -1 ! matIdx = material index at collision, nucIdx = nuclide index of target, MT = MT number of reaction.
    real(defReal)               :: A = ZERO, implicitWeight = ZERO, initialWeight = ZERO, kT = ZERO, k_eff = ZERO, &
                                   muL = ZERO, sigma_elasticScatter = ZERO, sigma_fission = ZERO, &
                                   sigma_inelasticScatter = ZERO, sigma_nuFiss = ZERO, sigma_tot = ZERO, weight = ZERO     ! A = target mass [neutron mass], E = collision energy (could be relative to target) [MeV], kT = target temperature [MeV], muL = cosine of deflection angle in LAB-frame (-).
    real(defReal), dimension(3) :: r = ZERO, u = ZERO
    type(RNG), pointer          :: RNGPtr => null()
  end type

end module collisionData_class