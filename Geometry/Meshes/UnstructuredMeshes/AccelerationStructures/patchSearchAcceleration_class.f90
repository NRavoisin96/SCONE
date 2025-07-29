module patchSearchAcceleration_class

  use accelerationStructure_inter,  only : accelerationStructure, initAccelerationStructurePayload
  use cartesianGrid_class,          only : cartesianGrid
  use dictionary_class,             only : dictionary
  use face_class,                   only : faceBox
  use genericProcedures,            only : fatalError
  use numPrecision
  use publicObjects,                only : coordData
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
    procedure :: findEntranceBoundaryFace
    procedure :: findHostElement
    procedure :: init
    procedure :: kill
  end type patchSearchAcceleration

contains
  !!
  !!
  !!
  subroutine findEntranceBoundaryFace(self, faces, data, boundaryFace)
    class(patchSearchAcceleration), intent(in) :: self
    type(topologicalObjectShelf), intent(in)   :: faces
    type(coordData), intent(inout)             :: data
    type(faceBox), intent(out)                 :: boundaryFace
    character(*), parameter                    :: here = 'distanceToBoundaryFace (patchSearchAcceleration_class.f90)'

    ! Call fatalError for now.
    call fatalError(here, 'Unsupported procedure.')

  end subroutine findEntranceBoundaryFace

  !!
  !!
  !!
  subroutine findHostElement(self, elements, data, stopSearch)
    class(patchSearchAcceleration), intent(in) :: self
    type(topologicalObjectShelf), intent(in)   :: elements
    type(coordData), intent(inout)             :: data
    logical(defBool), intent(out)              :: stopSearch
    character(*), parameter                    :: here = 'findHostElement (patchSearchAcceleration_class.f90)'

    stopSearch = .true.
    call fatalError(here, 'Unsupported procedure.')

  end subroutine findHostElement

  !!
  !!
  !!
  subroutine init(self, payload)
    class(patchSearchAcceleration), intent(inout)      :: self
    type(initAccelerationStructurePayload), intent(in) :: payload

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