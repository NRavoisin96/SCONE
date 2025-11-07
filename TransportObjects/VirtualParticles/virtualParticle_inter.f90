module virtualParticle_inter

  use transportObject_inter,      only : transportObject
  use transportObjectState_class, only : transportObjectState

  implicit none
  private

  !!
  !!
  !!
  type, public, abstract, extends(transportObject) :: virtualParticle
    private
  end type virtualParticle

end module virtualParticle_inter