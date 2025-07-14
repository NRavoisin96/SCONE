module accelerationStructure_inter

  use coord_class,        only : coord
  use dictionary_class,   only : dictionary
  use edgeShelf_class,    only : edgeShelf
  use elementShelf_class, only : elementShelf
  use faceShelf_class,    only : faceShelf
  use vertexShelf_class,  only : vertexShelf

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
    subroutine findHostElement(self, elements, coords)
      import :: accelerationStructure, coord, faceShelf, elementShelf
      class(accelerationStructure), intent(in) :: self
      type(elementShelf), intent(in)           :: elements
      type(coord), intent(inout)               :: coords
    end subroutine findHostElement

    !!
    !!
    !!
    subroutine init(self, dict, edges, elements, faces, vertices)
      import :: accelerationStructure, dictionary, edgeShelf, elementShelf, faceShelf, vertexShelf
      class(accelerationStructure), intent(inout) :: self
      class(dictionary), intent(in)               :: dict
      type(edgeShelf), intent(in)                 :: edges
      type(elementShelf), intent(in)              :: elements
      type(faceShelf), target, intent(in)         :: faces
      type(vertexShelf), intent(in)               :: vertices
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