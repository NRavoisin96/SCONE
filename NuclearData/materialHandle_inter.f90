module materialHandle_inter

  use errors_mod, only : fatalError
  use numPrecision

  implicit none
  private

  !!
  !! Material Handle allows to interact with different types of materials
  !!
  !! This is top abstract class for all materials.
  !! Diffrent types of materials like MG Neutron, CE Neutron etc. are its subclasses
  !!
  !! Interface:
  !!   kill -> returns to uninitialised state
  !!
  type, public, abstract :: materialHandle
    private
    real(defReal) :: inverseDensity = ZERO
  contains
    procedure                 :: getInverseDensity
    procedure(kill), deferred :: kill
    procedure                 :: setInverseDensity
  end type materialHandle

  abstract interface
    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      import :: materialHandle
      class(materialHandle), intent(inout) :: self
    end subroutine kill

  end interface

contains
  !!
  !!
  !!
  elemental function getInverseDensity(self) result(inverseDensity)
    class(materialHandle), intent(in) :: self
    real(defReal)                     :: inverseDensity

    inverseDensity = self % inverseDensity

  end function getInverseDensity

  !!
  !!
  !!
  elemental subroutine setInverseDensity(self, inverseDensity)
    class(materialHandle), intent(inout) :: self
    real(defReal), intent(in)            :: inverseDensity

    self % inverseDensity = inverseDensity

  end subroutine setInverseDensity

end module materialHandle_inter