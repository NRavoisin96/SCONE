module triangleShelf_class
  
  use genericProcedures, only : findCommon
  use numPrecision
  use triangle_class,    only : triangle
  
  implicit none
  private
  
  !!
  !! Storage space for triangles in a given OpenFOAM mesh.
  !!
  !! Private members:
  !!   shelf -> An array of triangles.
  !!
  type, public :: triangleShelf
    private
    type(triangle), dimension(:), allocatable, public :: shelf
  contains
    procedure                                         :: findCommonEdgeIdx
    procedure                                         :: kill
  end type triangleShelf

contains

  !! Function 'findCommonEdgeIdx'
  !!
  !! Basic description:
  !!   Returns the index of the common edge between two triangles.
  !!
  !! Arguments:
  !!   firstTriangleIdx [in]  -> Index of the first triangle.
  !!   secondTriangleIdx [in] -> Index of the second triangle.
  !!
  !! Result:
  !!   edgeIdx                -> Index of the common edge between the two triangles.
  !!
  elemental function findCommonEdgeIdx(self, firstTriangleIdx, secondTriangleIdx) result(edgeIdx)
    class(triangleShelf), intent(in) :: self
    integer(shortInt), intent(in)    :: firstTriangleIdx, secondTriangleIdx
    integer(shortInt)                :: edgeIdx
    integer(shortInt), dimension(:), allocatable :: commonIdxs

    ! Initialise edgeIdx = 0 then find common edge indices between the two triangles.
    edgeIdx = 0
    commonIdxs = findCommon(self % shelf(firstTriangleIdx) % getEdgeIdxs(), self % shelf(secondTriangleIdx) % getEdgeIdxs())
    
    ! If a common edge has been found update edgeIdx.
    if (size(commonIdxs) > 0) edgeIdx = commonIdxs(1)

  end function findCommonEdgeIdx
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(triangleShelf), intent(inout) :: self
    
    if (allocated(self % shelf)) deallocate(self % shelf)

  end subroutine kill
  
end module triangleShelf_class