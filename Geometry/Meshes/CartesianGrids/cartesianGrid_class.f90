module cartesianGrid_class

  use coord_class,      only : coord
  use dictionary_class, only : dictionary
  use faceShelf_class,  only : faceShelf
  use kdTree_class,     only : kdTree
  use mesh_inter,       only : mesh
  use numPrecision

  implicit none
  private

  type, public :: cartesianGrid
    private
    type(kdTree)              :: faceCentroidsTree
  contains
    ! Build procedures.
    procedure                 :: init
    ! Runtime procedures.
    procedure                 :: kill
  end type cartesianGrid

contains

  !!
  !!
  !!
  subroutine init(self, faces)
    class(cartesianGrid), intent(inout) :: self
    type(faceShelf), intent(in)         :: faces

    ! Initialise k-d tree for the unstructured mesh faces.
    ! call self % faceCentroidsTree % init(faces % getAllFaceCentroids(), .true.)

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(cartesianGrid), intent(inout) :: self

    ! Local.
    call self % faceCentroidsTree % kill()

  end subroutine kill

end module cartesianGrid_class