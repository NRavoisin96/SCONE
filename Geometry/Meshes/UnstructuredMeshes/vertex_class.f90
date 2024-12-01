module vertex_class
  
  use numPrecision
  use universalVariables
  use genericProcedures, only : append, linFind
  
  implicit none
  private
  
  !!
  !! Vertex of a given OpenFOAM mesh.
  !!
  !! Private members:
  !!   idx             -> Index of the vertex.
  !!   coordinates     -> 3-D coordinates of the vertex.
  !!   faceIdxs        -> Array of indices of the faces sharing the vertex.
  !!   elementIdxs     -> Array of indices of the elements sharing the vertex.
  !!   tetrahedronIdxs -> Array of indices of the tetrahedra sharing the vertex.
  !!   triangleIdxs    -> Array of indices of the triangles sharing the vertex.
  !!
  type, public :: vertex
    private
    integer(shortInt)                            :: idx = 0
    real(defReal), dimension(3)                  :: coordinates = ZERO
    integer(shortInt), dimension(:), allocatable :: faceIdxs, elementIdxs, &
                                                    tetrahedronIdxs, triangleIdxs
  contains
    procedure                                    :: addFaceIdx
    procedure                                    :: addElementIdx
    procedure                                    :: addTetrahedronIdx
    procedure                                    :: addTriangleIdx
    procedure                                    :: getCoordinates
    procedure                                    :: getIdx
    procedure                                    :: getVertexToElements
    procedure                                    :: getVertexToFaces
    procedure                                    :: getVertexToTetrahedra
    procedure                                    :: getVertexToTriangles
    procedure                                    :: kill
    procedure                                    :: setCoordinates
    procedure                                    :: setIdx
  end type

contains
  
  !! Subroutine 'addTetrahedronIdx'
  !!
  !! Basic description:
  !!   Adds the index of a tetrahedron sharing the vertex.
  !!
  !! Arguments:
  !!   tetrahedronIdx [in] -> Index of the tetrahedron.
  !!
  elemental subroutine addTetrahedronIdx(self, tetrahedronIdx)
    class(vertex), intent(inout)  :: self
    integer(shortInt), intent(in) :: tetrahedronIdx
    
    call append(self % tetrahedronIdxs, tetrahedronIdx)
  end subroutine addTetrahedronIdx
  
  !! Subroutine 'addTriangleIdx'
  !!
  !! Basic description:
  !!   Adds the index of a triangle sharing the vertex.
  !!
  !! Arguments:
  !!   triangleIdx [in] -> Index of the triangle.
  !!
  elemental subroutine addTriangleIdx(self, triangleIdx)
    class(vertex), intent(inout)  :: self
    integer(shortInt), intent(in) :: triangleIdx
    
    call append(self % triangleIdxs, triangleIdx)
  end subroutine addTriangleIdx
  
  !! Subroutine 'addFaceIdx'
  !!
  !! Basic description:
  !!   Adds the index of a face sharing the vertex.
  !!
  !! Arguments:
  !!   faceIdx [in] -> Index of the face.
  !!
  elemental subroutine addFaceIdx(self, faceIdx)
    class(vertex), intent(inout)  :: self
    integer(shortInt), intent(in) :: faceIdx
    
    call append(self % faceIdxs, faceIdx)
  end subroutine addFaceIdx
  
  !! Subroutine 'addElementIdx'
  !!
  !! Basic description:
  !!   Adds the index of an element sharing the vertex.
  !!
  !! Notes:
  !!   Due to the nature of the mesh importation process, here the index is only added if not
  !!   already present so as to avoid duplicates.
  !!
  !! Arguments:
  !!   elementIdx [in] -> Index of the element.
  !!
  !! TODO: Check if mesh importation process can be refactored to change this.
  elemental subroutine addElementIdx(self, elementIdx)
    class(vertex), intent(inout)  :: self
    integer(shortInt), intent(in) :: elementIdx
    
    if (.not. allocated(self % elementIdxs)) then
      allocate(self % elementIdxs(1))
      self % elementIdxs(1) = elementIdx
      return 

    end if

    if (linFind(self % elementIdxs, elementIdx) == targetNotFound) call append(self % elementIdxs, elementIdx)

  end subroutine addElementIdx
  
  !! Function 'getCoordinates'
  !!
  !! Basic description:
  !!   Returns the 3-D coordinates of the vertex.
  !!
  !! Result:
  !!   coordinates -> Array listing the x-, y- and z-coordinates of the vertex.
  !!
  pure function getCoordinates(self) result(coordinates)
    class(vertex), intent(in)   :: self
    real(defReal), dimension(3) :: coordinates
    
    coordinates = self % coordinates
  end function getCoordinates
  
  !! Function 'getIdx'
  !!
  !! Basic description:
  !!   Returns the index of the vertex.
  !!
  !! Result:
  !!   idx -> Index of the vertex.
  !!
  elemental function getIdx(self) result(idx)
    class(vertex), intent(in) :: self
    integer(shortInt)         :: idx
    
    idx = self % idx
  end function getIdx
  
  !! Function 'getVertexToElements'
  !!
  !! Basic description:
  !!   Returns the elements containing the vertex.
  !!
  !! Result:
  !!   elementIdxs -> Array listing the indices of the elements containing the vertex.
  !!
  pure function getVertexToElements(self) result(elementIdxs)
    class(vertex), intent(in)                              :: self
    integer(shortInt), dimension(size(self % elementIdxs)) :: elementIdxs
    
    elementIdxs = self % elementIdxs
  end function getVertexToElements
  
  !! Function 'getVertexToFaces'
  !!
  !! Basic description:
  !!   Returns the faces containing the vertex.
  !!
  !! Result:
  !!   faceIdxs -> Array listing the indices of the faces containing the vertex.
  !!
  pure function getVertexToFaces(self) result(faceIdxs)
    class(vertex), intent(in)                           :: self
    integer(shortInt), dimension(size(self % faceIdxs)) :: faceIdxs
    
    faceIdxs = self % faceIdxs
  end function getVertexToFaces
  !! Function 'getVertexToTetrahedra'
  !!
  !! Basic description:
  !!   Returns the tetrahedra containing the vertex.
  !!
  !! Result:
  !!   tetrahedronIdxs -> Array listing the indices of the tetrahedra containing the vertex.
  !!
  pure function getVertexToTetrahedra(self) result(tetrahedronIdxs)
    class(vertex), intent(in)                                  :: self
    integer(shortInt), dimension(size(self % tetrahedronIdxs)) :: tetrahedronIdxs
    
    tetrahedronIdxs = self % tetrahedronIdxs
  end function getVertexToTetrahedra
  
  !! Function 'getVertexToTriangles'
  !!
  !! Basic description:
  !!   Returns the triangles containing the vertex.
  !!
  !! Result:
  !!   triangleIdxs -> Array listing the indices of the triangles containing the vertex.
  !!
  pure function getVertexToTriangles(self) result(triangleIdxs)
    class(vertex), intent(in)                               :: self
    integer(shortInt), dimension(size(self % triangleIdxs)) :: triangleIdxs
    
    triangleIdxs = self % triangleIdxs
  end function getVertexToTriangles
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(vertex), intent(inout) :: self
    
    self % idx = 0
    self % coordinates = ZERO
    if (allocated(self % faceIdxs)) deallocate(self % faceIdxs)
    if (allocated(self % elementIdxs)) deallocate(self % elementIdxs)
    if (allocated(self % triangleIdxs)) deallocate(self % triangleIdxs)
    if (allocated(self % tetrahedronIdxs)) deallocate(self % tetrahedronIdxs)
  end subroutine kill
  
  !! Subroutine 'setCoordinates'
  !!
  !! Basic description:
  !!   Sets the 3-D coordinates of the vertex.
  !!
  !! Arguments:
  !!   coordinates [in] -> 3-D coordinates of the vertex.
  !!
  pure subroutine setCoordinates(self, coordinates)
    class(vertex), intent(inout)            :: self
    real(defReal), dimension(3), intent(in) :: coordinates
    
    self % coordinates = coordinates
  end subroutine setCoordinates
  
  !! Subroutine 'setIdx'
  !!
  !! Basic description:
  !!   Sets the index of the vertex.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the vertex.
  !!
  elemental subroutine setIdx(self, idx)
    class(vertex), intent(inout)  :: self
    integer(shortInt), intent(in) :: idx
    
    self % idx = idx
  end subroutine setIdx
end module vertex_class