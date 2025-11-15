module neutronMGstd_class
  
  use neutronMGCollisionProcessor_inter, only : neutronMGCollisionProcessor

  implicit none
  private

  !!
  !! Standard (default) scalar collision processor for MG neutrons
  !!   -> Preforms implicit fission site generation
  !!   -> Preforms analog capture
  !!   -> Treats fission as capture (only implicit generation of 2nd-ary neutrons)
  !!   -> Does not create secondary non-neutron projectiles
  !!
  !! Settings:
  !!  NONE
  !!
  !! Sample dictionary input:
  !!   collProcName {
  !!   type            neutronMGstd;
  !!   }
  !!
  type, public, extends(neutronMGCollisionProcessor) :: neutronMGstd
    private
  end type neutronMGstd

end module neutronMGstd_class