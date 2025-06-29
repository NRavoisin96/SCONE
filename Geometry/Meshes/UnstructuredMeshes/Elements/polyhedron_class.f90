module polyhedron_class
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use edge_class,                   only : edgeBox
  use edgeShelf_class,              only : edgeShelf
  use element_inter,                only : element, elementBox, kill_super => kill
  use face_inter,                   only : faceBox
  use faceShelf_class,              only : faceShelf
  use genericProcedures,            only : append, areEqual, computePyramidCentre, computePyramidVolume, &
                                           computeTetrahedronCentre, computeTetrahedronVolume, findCommon, &
                                           fatalError, numToChar
  use numPrecision
  use tetrahedron_class,            only : tetrahedron
  use triangle_class,               only : triangle
  use universalVariables,           only : SURF_TOL, INF, ZERO
  use vertex_class,                 only : vertexBox
  use vertexShelf_class,            only : vertexShelf
  
  implicit none
  private
  
  !!
  !! Element (cell) of an OpenFOAM mesh. Consists of a list of vertices and faces indices, as well
  !! as a list of tetrahedra indices into which the element is decomposed.
  !!
  !! Private members:
  !!   idx      -> Index of the element.
  !!   vertices -> Array of vertices indices making the element up.
  !!   faces    -> Array of faces indices making the element up.
  !!   Volume   -> Volume of the element.
  !!   Centroid -> Vector pointing to the centroid of the element.
  !!
  type, public, extends(element) :: polyhedron
    private
  contains
    ! Build procedures.
    procedure                    :: computeComponents
    procedure                    :: split
    ! Runtime procedures.
    procedure                    :: kill
  end type polyhedron

contains
  
  !! Subroutine 'computeVolumeAndCentroid'
  !!
  !! Basic description:
  !!   Computes the volume and volume-weighted centroid of the element.
  !!
  !! Detailed description:
  !!   First estimates the centroid of the element by taking an area-weighted average of the faces' 
  !!   centroids. Using this estimate, the element is decomposed into a number of pyramids whose 
  !!   apices are the estimated centroid of the element. Looping through all the faces, the volume 
  !!   and centroid of each pyramid is computed (see pyramid_class for more details); the centroid 
  !!   of the element is then obtained by taking a volume-weighted average of the pyramids' 
  !!   centroids.
  !!
  !! Arguments:
  !!   faces [in]    -> An array of face structures making the element up.
  !!   vertices [in] -> An array of vertex structures making the element up.
  !!
  !! Error:
  !!   fatalError if the volume is element is negative or infinite.
  !!
  pure subroutine computeComponents(self, faceIdxs, faces, vertices, centroid, volume)
    class(polyhedron), intent(inout)            :: self
    integer(shortInt), dimension(:), intent(in) :: faceIdxs
    type(faceShelf), intent(in)                 :: faces
    type(vertexBox), dimension(:), intent(in)   :: vertices
    real(defReal), dimension(3), intent(out)    :: centroid
    real(defReal), intent(out)                  :: volume
    integer(shortInt)                           :: i, nFaces, absFaceIdx
    real(defReal)                               :: area, sumAreas, sumVolumes
    real(defReal), dimension(3)                 :: sumVolumesCentroid, C
    real(defReal), dimension(:, :), allocatable :: array
    character(*), parameter                     :: Here = 'computeVolumeAndCentroid (element_class.f90)'

    ! Retrieve the number of faces and vertices in the element and initialise variables.
    nFaces = size(faceIdxs)
    C = ZERO
    sumAreas = ZERO
    
    ! Loop through all faces and compute the element's approximate centroid
    ! by performing an area-weighted average of the different faces' centroids.
    do i = 1, nFaces
      ! Retrieve the area of the current face.
      absFaceIdx = abs(faceIdxs(i))
      area = faces % getFaceArea(absFaceIdx)
      ! Update the area-weighted centroid and the sum of faces' areas.
      C = C + faces % getFaceCentroid(absFaceIdx) * area
      sumAreas = sumAreas + area

    end do
    
    ! Using the approximate centroid, compute the actual centroid by performing a volume-weighted
    ! average of the different pyramids' centroids.
    sumVolumes = ZERO
    sumVolumesCentroid = ZERO
    allocate(array(2, 3))
    C = C / sumAreas
    
    ! Loop through all faces (pyramids).
    do i = 1, nFaces
      ! Retrieve the volume of the current pyramid and update the volume-weighted centroid and the sum of volumes.
      absFaceIdx = abs(faceIdxs(i))
      array(1, :) = faces % getFaceNormal(absFaceIdx) * faces % getFaceArea(absFaceIdx)
      array(2, :) = C - faces % getFaceCentroid(absFaceIdx)
      volume = computePyramidVolume(array)

      array(1, :) = 3.0_defReal * faces % getFaceCentroid(absFaceIdx)
      array(2, :) = C
      sumVolumesCentroid = sumVolumesCentroid + computePyramidCentre(array) * volume
      sumVolumes = sumVolumes + volume

    end do
    ! The volume of the element is simply the sum of volumes, while the centroid is the average of
    ! the volume-weighted sum.
    volume = sumVolumes
    centroid = sumVolumesCentroid / sumVolumes

  end subroutine computeComponents
  
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  subroutine kill(self)
    class(polyhedron), intent(inout) :: self
    
    ! Element.
    call kill_super(self)

  end subroutine kill
  
  !! Subroutine 'split'
  !!
  !! Basic description:
  !!   Splits the element into a number of pyramids corresponding to the number of faces in the 
  !!   element.
  !!
  !! Detailed description:
  !!   The element is split into pyramids by subdividing it from its centroid: this then forms the
  !!   apex of each generated pyramid. For each face in the element, triangles are created by
  !!   joining each edge in the face to the common apex. During the creation of triangles a check is
  !!   made to ensure that their normal vectors point in the correct direction (from the pyramid of
  !!   lowest index to that of greatest index, just as OpenFOAM does for elements).
  !!
  !! Arguments:
  !!   faces [in]              -> A faceShelf.
  !!   edges [inout]           -> An edgeShelf.
  !!   vertices [inout]        -> A vertexShelf.
  !!   triangles [inout]       -> A triangleShelf.
  !!   pyramids [inout]        -> A pyramidShelf.
  !!   lastEdgeIdx [inout]     -> Index of the first free item in the edgeShelf.
  !!   lastTriangleIdx [inout] -> Index of the first free item in the triangleShelf.
  !!   lastPyramidIdx [inout]  -> Index of the first free item in the pyramidShelf.
  !!   lastVertexIdx [in]      -> Index of the first free item in the vertexShelf.
  !!
  subroutine split(self, faces, lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx, &
                        newEdges, newFaces, newVertices, tetrahedra, triangles)
    class(polyhedron), intent(inout)              :: self
    type(faceShelf), intent(inout)                :: faces, newFaces
    integer(shortInt), intent(inout)              :: lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx
    type(edgeShelf), intent(inout)                :: newEdges
    type(vertexShelf), intent(inout)              :: newVertices
    type(elementBox), dimension(:), intent(inout) :: tetrahedra
    type(faceBox), dimension(:), intent(inout)    :: triangles
    integer(shortInt)                             :: i, j, k, l, faceIdx, absFaceIdx, faceTriangleIdx, triangleVertexIdx, &
                                                     nFaces, nVertices, edgeIdx, commonTriangleIdx
    integer(shortInt), dimension(:), allocatable  :: faceIdxs, faceTriangleIdxs, vertexIdxs
    integer(shortInt), dimension(2)               :: edgeVertexIdxs
    integer(shortInt), dimension(3)               :: testVertexIdxs
    integer(shortInt), dimension(4)               :: triangleIdxs
    integer(shortInt), dimension(6)               :: edgeIdxs
    real(defReal)                                 :: volume
    real(defReal), dimension(3)                   :: centroid
    real(defReal), dimension(4, 3)                :: array
    real(defReal), dimension(:, :), allocatable   :: vertexCoords
    type(axisAlignedBoundingBox)                  :: boundingBox
    type(vertexBox), dimension(2)                 :: edgeVertices
    type(vertexBox), dimension(3)                 :: triangleVertices
    
    
  
  end subroutine split

end module polyhedron_class