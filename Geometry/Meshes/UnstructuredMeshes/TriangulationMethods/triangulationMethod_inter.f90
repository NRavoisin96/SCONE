module triangulationMethod_inter

  use edgeShelf_class,    only : edgeShelf
  use elementShelf_class, only : elementShelf
  use faceShelf_class,    only : faceShelf
  use numPrecision
  use vertexShelf_class,  only : vertexShelf

  implicit none
  private

  !!
  !!
  !!
  type, public, abstract :: triangulationMethod
    private
  contains
    procedure(triangulate), deferred :: triangulate
  end type triangulationMethod

  !!
  !!
  !!
  abstract interface
    !!
    !!
    !!
    subroutine triangulate(self, edges, elements, faces, vertices)
      import                                 :: edgeShelf, elementShelf, faceShelf, triangulationMethod, vertexShelf
      class(triangulationMethod), intent(in) :: self
      type(edgeShelf), intent(inout)         :: edges
      type(elementShelf), intent(inout)      :: elements
      type(faceShelf), intent(inout)         :: faces
      type(vertexShelf), intent(inout)       :: vertices
    end subroutine triangulate

  end interface

end module triangulationMethod_inter