module accelerationStructure_inter

  use dictionary_class,             only : dictionary
  use element_class,                only : elementBox
  use face_class,                   only : faceBox
  use numPrecision
  use publicObjects,                only : coordData
  use topologicalObjectShelf_class, only : topologicalObjectShelf
  use universalVariables,           only : VALENCE

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

  !!
  !!
  !!
  type, public :: initAccelerationStructurePayload
    class(dictionary), pointer            :: dict => null()
    type(topologicalObjectShelf), pointer :: edges => null(), elements => null(), faces => null(), vertices => null()
  end type initAccelerationStructurePayload

  abstract interface
    !!
    !!
    !!
    subroutine findEntranceBoundaryFace(self, faces, data, nIntersectedFaces, intersectedFaceIdxs)
      import :: accelerationStructure, coordData, shortInt, topologicalObjectShelf, VALENCE
      class(accelerationStructure), intent(in)           :: self
      type(topologicalObjectShelf), intent(in)           :: faces
      type(coordData), intent(inout)                     :: data
      integer(shortInt), intent(out)                     :: nIntersectedFaces
      integer(shortInt), dimension(VALENCE), intent(out) :: intersectedFaceIdxs
    end subroutine findEntranceBoundaryFace

    !!
    !!
    !!
    subroutine findHostElement(self, elements, data, stopSearch)
      import :: accelerationStructure, coordData, defBool, topologicalObjectShelf
      class(accelerationStructure), intent(in) :: self
      type(topologicalObjectShelf), intent(in) :: elements
      type(coordData), intent(inout)           :: data
      logical(defBool), intent(out)            :: stopSearch
    end subroutine findHostElement

    !!
    !!
    !!
    subroutine init(self, payload)
      import :: accelerationStructure, initAccelerationStructurePayload
      class(accelerationStructure), intent(inout)        :: self
      type(initAccelerationStructurePayload), intent(in) :: payload
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