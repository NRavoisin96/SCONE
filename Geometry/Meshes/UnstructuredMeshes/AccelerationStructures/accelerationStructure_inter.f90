module accelerationStructure_inter

  use dictionary_class,   only : dictionary
  use edgeShelf_class,    only : edgeShelf
  use elementShelf_class, only : elementShelf
  use faceShelf_class,    only : faceShelf
  use numPrecision
  use vertexShelf_class,  only : vertexShelf

  implicit none
  private

  ! Public procedures.
  public :: getStorageSize

  !!
  !!
  !!
  type, public, abstract :: accelerationStructure
    private
  contains
    procedure(findHostElementIdx), deferred :: findHostElementIdx
    procedure                               :: getStorageSize
    procedure(init), deferred               :: init
    procedure(kill), deferred               :: kill
  end type accelerationStructure

  !!
  !!
  !!
  abstract interface
    !!
    !!
    !!
    subroutine findHostElementIdx(self, u, edges, elements, faces, vertices, elementIdx, r)
      import :: accelerationStructure, defReal, edgeShelf, elementShelf, faceShelf, shortInt, vertexShelf
      class(accelerationStructure), intent(in)   :: self
      real(defReal), dimension(3), intent(in)    :: u
      type(edgeShelf), intent(in)                :: edges
      type(elementShelf), intent(in)             :: elements
      type(faceShelf), intent(in)                :: faces
      type(vertexShelf), intent(in)              :: vertices
      integer(shortInt), intent(inout)           :: elementIdx
      real(defReal), dimension(3), intent(inout) :: r
    end subroutine findHostElementIdx

    !!
    !!
    !!
    subroutine init(self, dict, vertices, edges, faces, elements)
      import :: accelerationStructure, dictionary, elementShelf, faceShelf, vertexShelf, edgeShelf
      class(accelerationStructure), intent(inout) :: self
      type(dictionary), intent(in)                :: dict
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

contains
  !!
  !!
  !!
  elemental function getStorageSize(self) result(storageSize)
    class(accelerationStructure), intent(in) :: self
    integer(longInt)                         :: storageSize

    storageSize = storage_size(self) / 8

  end function getStorageSize

end module accelerationStructure_inter