module edge_class
  
  use numPrecision
  use genericProcedures, only : append
  
  implicit none
  private
  
  !!
  !! Edge of a mesh linking two vertices.
  !!
  !! Private members:
  !!   idx            -> Index of the edge.
  !!   startVertexIdx -> Index of the first vertex in the edge.
  !!   endVertexIdx   -> Index of the end vertex in the edge.
  !!   edgeToFaces    -> Array that stores edge-to-faces connectivity information.
  !!   edgeToElements -> Array that stores edge-to-elements connectivity information.
  !!
  type, public                                   :: edge
    private
    integer(shortInt)                            :: idx = 0
    integer(shortInt), dimension(2)              :: vertexIdxs = 0
    integer(shortInt), dimension(:), allocatable :: faceIdxs, elementIdxs, tetrahedronIdxs, triangleIdxs
  contains
    ! Build procedures.
    procedure                                    :: addElementIdx
    procedure                                    :: addFaceIdx
    procedure                                    :: addTetrahedronIdx
    procedure                                    :: addTriangleIdx
    procedure                                    :: kill
    procedure                                    :: setIdx
    procedure                                    :: setVertexIdxs
    ! Runtime procedures.
    procedure                                    :: getElementIdxs
    procedure                                    :: getFaceIdxs
    procedure                                    :: getIdx
    procedure                                    :: getTetrahedronIdxs
    procedure                                    :: getTriangleIdxs
    procedure                                    :: getVertexIdxs
  end type edge

contains

  !! Subroutine 'addElementIdx'
  !!
  !! Basic description:
  !!   Adds the index of an element sharing the edge.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element.
  !!
  elemental subroutine addElementIdx(self, idx)
    class(edge), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx

    call append(self % elementIdxs, idx, .true.)

  end subroutine addElementIdx

  !! Subroutine 'addFaceIdx'
  !!
  !! Basic description:
  !!   Adds the index of a face sharing the edge.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the face.
  !!
  elemental subroutine addFaceIdx(self, idx)
    class(edge), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx

    call append(self % faceIdxs, idx)

  end subroutine addFaceIdx

  !! Subroutine 'addTetrahedronIdx'
  !!
  !! Basic description:
  !!   Adds the index of a tetrahedron sharing the edge.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the tetrahedron.
  !!
  elemental subroutine addTetrahedronIdx(self, idx)
    class(edge), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx

    call append(self % tetrahedronIdxs, idx, .true.)

  end subroutine addTetrahedronIdx

  !! Subroutine 'addTriangleIdx'
  !!
  !! Basic description:
  !!   Adds the index of a triangle sharing the edge.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the triangle.
  !!
  elemental subroutine addTriangleIdx(self, idx)
    class(edge), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx

    call append(self % triangleIdxs, idx)

  end subroutine addTriangleIdx

  !! Function 'getElementIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the elements sharing the edge.
  !!
  !! Result:
  !!   elementIdxs -> Indices of the elements sharing the edge.
  !!
  pure function getElementIdxs(self) result(faceIdxs)
    class(edge), intent(in)                                :: self
    integer(shortInt), dimension(size(self % elementIdxs)) :: faceIdxs

    faceIdxs = self % elementIdxs

  end function getElementIdxs

  !! Function 'getFaceIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the faces sharing the edge.
  !!
  !! Result:
  !!   faceIdxs -> Indices of the faces sharing the edge.
  !!
  pure function getFaceIdxs(self) result(faceIdxs)
    class(edge), intent(in)                             :: self
    integer(shortInt), dimension(size(self % faceIdxs)) :: faceIdxs

    faceIdxs = self % faceIdxs

  end function getFaceIdxs

  !! Function 'getIdx'
  !!
  !! Basic description:
  !!   Returns the index of the edge.
  !!
  !! Result:
  !!   idx -> Index of the edge.
  !!
  elemental function getIdx(self) result(idx)
    class(edge), intent(in) :: self
    integer(shortInt)       :: idx

    idx = self % idx

  end function getIdx

  !! Function 'getTetrahedronIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the tetrahedra sharing the edge.
  !!
  !! Result:
  !!   tetrahedronIdxs -> Indices of the tetrahedra sharing the edge.
  !!
  pure function getTetrahedronIdxs(self) result(tetrahedronIdxs)
    class(edge), intent(in)                                    :: self
    integer(shortInt), dimension(size(self % tetrahedronIdxs)) :: tetrahedronIdxs

    tetrahedronIdxs = self % tetrahedronIdxs

  end function getTetrahedronIdxs

  !! Function 'getTriangleIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the triangles sharing the edge.
  !!
  !! Result:
  !!   triangleIdxs -> Indices of the triangles sharing the edge.
  !!
  pure function getTriangleIdxs(self) result(triangleIdxs)
    class(edge), intent(in)                                 :: self
    integer(shortInt), dimension(size(self % triangleIdxs)) :: triangleIdxs

    triangleIdxs = self % triangleIdxs

  end function getTriangleIdxs

  !! Function 'getVertexIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in the edge.
  !!
  !! Result:
  !!   vertexIdxs -> Indices of the vertices in the edge.
  !!
  pure function getVertexIdxs(self) result(vertexIdxs)
    class(edge), intent(in)         :: self
    integer(shortInt), dimension(2) :: vertexIdxs

    vertexIdxs = self % vertexIdxs

  end function getVertexIdxs

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an unitialised state.
  !!
  elemental subroutine kill(self)
    class(edge), intent(inout) :: self

    self % idx = 0
    self % vertexIdxs = 0
    if (allocated(self % faceIdxs)) deallocate(self % faceIdxs)
    if (allocated(self % elementIdxs)) deallocate(self % elementIdxs)
    if (allocated(self % tetrahedronIdxs)) deallocate(self % tetrahedronIdxs)
    if (allocated(self % triangleIdxs)) deallocate(self % triangleIdxs)

  end subroutine kill

  !! Subroutine 'setIdx'
  !!
  !! Basic description:
  !!   Sets the index of the edge.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the edge.
  !!
  elemental subroutine setIdx(self, idx)
    class(edge), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx

    self % idx = idx

  end subroutine setIdx

  !! Subroutine 'setVertexIdxs'
  !!
  !! Basic description:
  !!   Sets the indices of the vertices in the edge.
  !!
  !! Arguments:
  !!   vertexIdxs [in] -> Array containing the indices of the vertices in the edge.
  !!
  pure subroutine setVertexIdxs(self, vertexIdxs)
    class(edge), intent(inout)                  :: self
    integer(shortInt), dimension(2), intent(in) :: vertexIdxs

    self % vertexIdxs = vertexIdxs

  end subroutine setVertexIdxs

end module edge_class