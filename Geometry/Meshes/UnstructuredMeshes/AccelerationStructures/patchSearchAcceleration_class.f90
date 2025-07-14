module patchSearchAcceleration_class

  use accelerationStructure_inter, only : accelerationStructure
  use cartesianGrid_class,         only : cartesianGrid
  use coord_class,                 only : coord
  use dictionary_class,            only : dictionary
  use edgeShelf_class,             only : edgeShelf
  use elementShelf_class,          only : elementShelf
  use faceShelf_class,             only : faceShelf
  use vertexShelf_class,           only : vertexShelf

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
    type(elementShelf), intent(in)             :: elements
    type(coord), intent(inout)                 :: coords

  end subroutine findHostElement

  !!
  !!
  !!
  subroutine init(self, dict, edges, elements, faces, vertices)
    class(patchSearchAcceleration), intent(inout) :: self
    class(dictionary), intent(in)                 :: dict
    type(edgeShelf), intent(in)                   :: edges
    type(elementShelf), intent(in)                :: elements
    type(faceShelf), target, intent(in)           :: faces
    type(vertexShelf), intent(in)                 :: vertices

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