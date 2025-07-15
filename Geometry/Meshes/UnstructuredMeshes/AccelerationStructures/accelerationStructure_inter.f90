module accelerationStructure_inter

  use coord_class,                  only : coord
  use dictionary_class,             only : dictionary
  use topologicalObjectShelf_class, only : topologicalObjectShelf

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
      import :: accelerationStructure, coord, topologicalObjectShelf
      class(accelerationStructure), intent(in) :: self
      type(topologicalObjectShelf), intent(in) :: elements
      type(coord), intent(inout)               :: coords
    end subroutine findHostElement

    !!
    !!
    !!
    subroutine init(self, dict, edges, elements, faces, vertices)
      import :: accelerationStructure, dictionary, topologicalObjectShelf
      class(accelerationStructure), intent(inout)      :: self
      class(dictionary), intent(in)                    :: dict
      type(topologicalObjectShelf), target, intent(in) :: edges, elements, faces, vertices
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