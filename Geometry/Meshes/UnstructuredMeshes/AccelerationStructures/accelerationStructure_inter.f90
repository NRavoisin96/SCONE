module accelerationStructure_inter

  use coord_class,                  only : coord
  use dictionary_class,             only : dictionary
  use element_class,                only : elementBox
  use face_class,                   only : faceBox
  use numPrecision
  use topologicalObjectShelf_class, only : topologicalObjectShelf

  implicit none
  private

  !!
  !!
  !!
  type, public, abstract :: accelerationStructure
    private
  contains
    procedure(findEntranceBoundaryFace), deferred :: findEntranceBoundaryFace
    procedure(findHostElement), deferred          :: findHostElement
    procedure(init), deferred                     :: init
    procedure(kill), deferred                     :: kill
  end type accelerationStructure

  abstract interface
    !!
    !!
    !!
    subroutine findEntranceBoundaryFace(self, faces, coords, d, boundaryFace)
      import :: accelerationStructure, coord, defReal, faceBox, topologicalObjectShelf
      class(accelerationStructure), intent(in) :: self
      type(topologicalObjectShelf), intent(in) :: faces
      type(coord), intent(in)                  :: coords
      real(defReal), intent(inout)             :: d
      type(faceBox), intent(out)               :: boundaryFace
    end subroutine findEntranceBoundaryFace

    !!
    !!
    !!
    subroutine findHostElement(self, elements, coords, stopSearch)
      import :: accelerationStructure, coord, defBool, topologicalObjectShelf
      class(accelerationStructure), intent(in) :: self
      type(topologicalObjectShelf), intent(in) :: elements
      type(coord), intent(inout)               :: coords
      logical(defBool), intent(out)            :: stopSearch
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