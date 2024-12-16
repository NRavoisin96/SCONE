module faceShelf_class
  
  use numPrecision
  use genericProcedures,   only : findCommon
  use face_class,          only : face
  
  implicit none
  private
  
  !!
  !! Storage space for faces of an OpenFOAM mesh.
  !!
  !! Private members:
  !!   shelf -> Array to store faces.
  !!
  type, public :: faceShelf
    private
    type(face), dimension(:), allocatable, public :: shelf
  contains
    procedure                                     :: findCommonEdgeIdx
    procedure                                     :: getSize
    procedure                                     :: kill
  end type

contains

  !! Function 'findCommonEdgeIdx'
  !!
  !! Basic description:
  !!   Returns the index of the common edge between two faces.
  !!
  !! Arguments:
  !!   firstFaceIdx [in]  -> Index of the first face.
  !!   secondFaceIdx [in] -> Index of the second face.
  !!
  !! Result:
  !!   edgeIdx            -> Index of the common edge between the two faces.
  !!
  elemental function findCommonEdgeIdx(self, firstFaceIdx, secondFaceIdx) result(edgeIdx)
    class(faceShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: firstFaceIdx, secondFaceIdx
    integer(shortInt)                            :: edgeIdx
    integer(shortInt), dimension(:), allocatable :: commonIdxs

    ! Initialise edgeIdx = 0 then find common edge indices between the two faces.
    edgeIdx = 0
    commonIdxs = findCommon(self % shelf(firstFaceIdx) % getEdgeIdxs(), self % shelf(secondFaceIdx) % getEdgeIdxs())

    ! If a common edge has been found update edgeIdx.
    if (size(commonIdxs) > 0) edgeIdx = commonIdxs(1)

  end function findCommonEdgeIdx
  
  !! Function 'getSize'
  !!
  !! Basic description:
  !!   Returns the size of the faceShelf.
  !!
  !! Result:
  !!   size -> Size of the faceShelf.
  !!
  elemental function getSize(self) result(nFaces)
    class(faceShelf), intent(in) :: self
    integer(shortInt)            :: nFaces
    
    nFaces = size(self % shelf)
    
  end function getSize
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(faceShelf), intent(inout) :: self
    
    if (allocated(self % shelf)) deallocate(self % shelf)

  end subroutine kill
  
end module faceShelf_class