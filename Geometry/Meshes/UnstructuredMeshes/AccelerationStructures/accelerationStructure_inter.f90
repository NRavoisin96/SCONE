module accelerationStructure_inter

  use coord_class,        only : coord
  use elementShelf_class, only : elementShelf
  use faceShelf_class,    only : faceShelf
  use vertexShelf_class,  only : vertexShelf
  use edgeShelf_class,    only : edgeShelf

  implicit none
  private

  !!
  !!
  !!
  type, public, abstract :: accelerationStructure
    private
  contains
    procedure(findHostElement), deferred :: findHostElement
    procedure(init), deferred            :: init
    procedure(kill), deferred            :: kill
  end type accelerationStructure

  !!
  !!
  !!
  abstract interface
    !!
    !!
    !!
    subroutine findHostElement(self, vertices, edges, faces, elements, coords)
      import :: accelerationStructure, coord, faceShelf, elementShelf, vertexShelf, edgeShelf
      class(accelerationStructure), intent(in) :: self
      class(vertexShelf), intent(in)           :: vertices
      class(edgeShelf), intent(in)             :: edges
      type(faceShelf), intent(in)              :: faces
      type(elementShelf), intent(in)           :: elements
      type(coord), intent(inout)               :: coords
    end subroutine findHostElement

    !!
    !!
    !!
    subroutine init(self, vertices, edges, faces, elements)
      import :: accelerationStructure, elementShelf, faceShelf, vertexShelf, edgeShelf
      class(accelerationStructure), intent(inout) :: self
      type(vertexShelf), intent(in)               :: vertices
      type(faceShelf), intent(inout)              :: faces
      type(elementShelf), intent(in)              :: elements
      type(edgeShelf), intent(inout)              :: edges
    end subroutine init

    !!
    !!
    !!
    elemental subroutine kill(self)
      import :: accelerationStructure
      class(accelerationStructure), intent(inout) :: self
    end subroutine kill

  end interface

end module accelerationStructure_inter