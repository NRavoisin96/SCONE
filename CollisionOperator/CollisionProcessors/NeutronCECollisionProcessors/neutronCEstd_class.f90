module neutronCEstd_class

  use neutronCECollisionProcessor_inter, only : neutronCECollisionProcessor
  use numPrecision

  implicit none
  private

  !!
  !! Standard (default) scalar collision processor for CE neutrons
  !!   -> Preforms implicit fission site generation
  !!   -> Preforms analog capture
  !!   -> Treats fission as capture (only implicit generation of 2nd-ary neutrons)
  !!   -> Does not create secondary non-neutron projectiles
  !!
  !! Settings:
  !!  minE    -> minimum energy cut-off [MeV] (default = 1.0E-11)
  !!  maxE    -> maximum energy. Higher energies are set to maximum (not re-rolled) [MeV]
  !!             (default = 20.0)
  !!  threshE -> Energy threshold for explicit treatment of target nuclide movement [-].
  !!             Target movement is sampled if neutron energy E < kT * threshE where
  !!             kT is target material temperature in [MeV]. (default = 400.0)
  !!  threshA -> Mass threshold for explicit treatment of target nuclide movement [Mn].
  !!             Target movement is sampled if target mass A < threshA. (default = 1.0)
  !!  DBRCeMin -> Minimum energy to which DBRC is applied
  !!  DBRCeMax -> Maximum energy to which DBRC is applied
  !!
  !! Sample dictionary input:
  !!   collProcName {
  !!   type             neutronCEstd;
  !!   #minEnergy       <real>;#
  !!   #maxEnergy       <real>;#
  !!   #energyThreshold <real>;#
  !!   #massThreshold   <real>;#
  !!   }
  !!
  type, public, extends(neutronCECollisionProcessor) :: neutronCEstd
    private
  contains
    procedure :: getImplicitCondition
  end type neutronCEstd

contains
  !!
  !!
  !!
  elemental function getImplicitCondition(self) result(isIt)
    class(neutronCEstd), intent(in) :: self
    logical(defBool)                :: isIt

    isIt = self % getNuclideIsFissile()

  end function getImplicitCondition

end module neutronCEstd_class