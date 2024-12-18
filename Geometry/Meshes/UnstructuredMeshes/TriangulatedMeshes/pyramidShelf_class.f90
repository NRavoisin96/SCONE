module pyramidShelf_class
  
  use edgeShelf_class,        only : edgeShelf
  use faceShelf_class,        only : faceShelf
  use numPrecision
  use pyramid_class,          only : pyramid
  use tetrahedronShelf_class, only : tetrahedronShelf
  use triangleShelf_class,    only : triangleShelf
  use vertexShelf_class,      only : vertexShelf
  
  implicit none
  private
  
  !!
  !! Storage space for pyramids in a given OpenFOAM mesh.
  !!
  !! Private members:
  !!   shelf -> Array of pyramids.
  !!
  type, public                               :: pyramidShelf
    private
    type(pyramid), dimension(:), allocatable :: shelf
  contains
    procedure                                :: addTriangleIdxToPyramid
    procedure                                :: allocateShelf
    procedure                                :: getPyramidCentroid
    procedure                                :: getPyramidFaceIdx
    procedure                                :: initPyramid
    procedure                                :: kill
    procedure                                :: splitPyramid
  end type pyramidShelf

contains

  !! Subroutine 'addTriangleIdxToPyramid'
  !!
  !! Basic description:
  !!   Adds the index of a triangle to a pyramid in the shelf.
  !!
  !! Arguments:
  !!   idx [in]         -> Index of the pyramid in the shelf.
  !!   triangleIdx [in] -> Index of the triangle in the pyramid.
  !!
  elemental subroutine addTriangleIdxToPyramid(self, idx, triangleIdx)
    class(pyramidShelf), intent(inout) :: self
    integer(shortInt), intent(in)      :: idx, triangleIdx

    call self % shelf(idx) % addTriangle(triangleIdx)

  end subroutine addTriangleIdxToPyramid

  !! Subroutine 'allocateShelf'
  !!
  !! Basic description:
  !!   Allocates memory in the shelf.
  !!
  !! Arguments:
  !!   nPyramids [in] -> Number of pyramids to be included in the shelf.
  !!
  elemental subroutine allocateShelf(self, nPyramids)
    class(pyramidShelf), intent(inout) :: self
    integer(shortInt), intent(in)      :: nPyramids

    allocate(self % shelf(nPyramids))

  end subroutine allocateShelf

  !! Function 'getPyramidCentroid'
  !!
  !! Basic description:
  !!   Returns the centroid of a pyramid in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the pyramid in the shelf.
  !!
  !! Result:
  !!   centroid -> 3-D coordinates of the centroid of the pyramid.
  !!
  pure function getPyramidCentroid(self, idx) result(centroid)
    class(pyramidShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    real(defReal), dimension(3)     :: centroid

    centroid = self % shelf(idx) % getCentroid()

  end function getPyramidCentroid

  !! Function 'getPyramidFaceIdx'
  !!
  !! Basic description:
  !!   Returns the index of the base face of a pyramid in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the pyramid in the shelf.
  !!
  !! Result:
  !!   faceIdx  -> Index of the base face of the pyramid.
  !!
  pure function getPyramidFaceIdx(self, idx) result(faceIdx)
    class(pyramidShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    integer(shortInt)               :: faceIdx

    faceIdx = self % shelf(idx) % getFace()

  end function getPyramidFaceIdx

  !! Subroutine 'initPyramid'
  !!
  !! Basic description:
  !!   Initialises a pyramid in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the pyramid in the shelf.
  !!   elementIdx [in] -> Index of the element from which the pyramid originates.
  !!   faceIdx [in]    -> Index of the base face of the pyramid.
  !!   vertexIdxs [in] -> Indices of the vertices of the pyramid.
  !!   apexCoords [in] -> 3-D coordinates of the apex vertex of the pyramid.
  !!   faces [in]      -> A faceShelf.
  !!
  pure subroutine initPyramid(self, idx, elementIdx, faceIdx, vertexIdxs, apexCoords, faces)
    class(pyramidShelf), intent(inout)          :: self
    integer(shortInt), intent(in)               :: idx, elementIdx, faceIdx
    integer(shortInt), dimension(:), intent(in) :: vertexIdxs
    real(defReal), dimension(3), intent(in)     :: apexCoords
    type(faceShelf), intent(in)                 :: faces

    call self % shelf(idx) % setIdx(idx)
    call self % shelf(idx) % setElement(elementIdx)
    call self % shelf(idx) % setFace(faceIdx)
    call self % shelf(idx) % setVertices(vertexIdxs)
    call self % shelf(idx) % computeCentroid(faces, apexCoords)

  end subroutine initPyramid
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(pyramidShelf), intent(inout) :: self
    
    if (allocated(self % shelf)) deallocate(self % shelf)

  end subroutine kill

  !! Subroutine 'splitPyramid'
  !!
  !! Basic description:
  !!   Splits a pyramid in the shelf into tetrahedra.
  !!
  !! Arguments:
  !!   idx [in]                   -> Index of the pyramid in the shelf.
  !!   edges [inout]              -> An edgeShelf.
  !!   faces [in]                 -> A faceShelf.
  !!   tetrahedra [inout]         -> A tetrahedronShelf.
  !!   triangles [inout]          -> A triangleShelf.
  !!   vertices [inout]           -> A vertexShelf.
  !!   lastTetrahedronIdx [inout] -> Index of the last tetrahedron in the tetrahedronShelf.
  !!   lastTriangleIdx [inout]    -> Index of the last triangle in the triangleShelf.
  !!
  elemental subroutine splitPyramid(self, idx, edges, faces, tetrahedra, triangles, vertices, lastTetrahedronIdx, lastTriangleIdx)
    class(pyramidShelf), intent(inout)    :: self
    integer(shortInt), intent(in)         :: idx
    type(edgeShelf), intent(inout)        :: edges
    type(faceShelf), intent(in)           :: faces
    type(tetrahedronShelf), intent(inout) :: tetrahedra
    type(triangleShelf), intent(inout)    :: triangles
    type(vertexShelf), intent(inout)      :: vertices
    integer(shortInt), intent(inout)      :: lastTetrahedronIdx, lastTriangleIdx

    call self % shelf(idx) % split(faces, edges, vertices, triangles, tetrahedra, lastTriangleIdx, lastTetrahedronIdx)

  end subroutine splitPyramid
  
end module pyramidShelf_class