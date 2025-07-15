module patchSearchAcceleration_class

  use accelerationStructure_inter,  only : accelerationStructure
  use cartesianGrid_class,          only : cartesianGrid
  use coord_class,                  only : coord
  use dictionary_class,             only : dictionary
  use topologicalObjectShelf_class, only : topologicalObjectShelf

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(accelerationStructure) :: patchSearchAcceleration
    private
    type(cartesianGrid)                        :: grid
  contains
    procedure :: findHostElement
    procedure :: init
    procedure :: kill
  end type patchSearchAcceleration

contains
  !!
  !!
  !!
  subroutine findHostElement(self, elements, coords)
    class(patchSearchAcceleration), intent(in) :: self
    type(topologicalObjectShelf), intent(in)   :: elements
    type(coord), intent(inout)                 :: coords

  end subroutine findHostElement

  !!
  !!
  !!
  subroutine init(self, dict, edges, elements, faces, vertices)
    class(patchSearchAcceleration), intent(inout)    :: self
    class(dictionary), intent(in)                    :: dict
    type(topologicalObjectShelf), target, intent(in) :: edges, elements, faces, vertices

    ! Initialise grid.
    call self % grid % init()

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(patchSearchAcceleration), intent(inout) :: self

    ! Local.
    call self % grid % kill()

  end subroutine kill

end module patchSearchAcceleration_class