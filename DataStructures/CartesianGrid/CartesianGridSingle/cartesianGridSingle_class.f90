module cartesianGridSingle_class
  
  use cartesianCellSingle_class, only : cartesianCellSingle
  use cartesianGenericProcedures
  use edgeShelf_class,           only : edgeShelf
  use elementShelf_class,        only : elementShelf
  use face_inter,                only : faceSATData
  use faceShelf_class,           only : faceShelf
  use genericProcedures,         only : crossProduct, fatalError
  use numPrecision
  use universalVariables,        only : INF
  use vertexShelf_class,         only : vertexShelf
  
  implicit none
  private
  
  !!
  !!
  type, public :: cartesianGridSingle
    private
    real(defReal)                                              :: wStar = ZERO, alpha = ZERO, l_min = ZERO, &
                                                                  spacingReciprocal = ZERO, spacing = ZERO
    real(defReal), dimension(3)                                :: gridBounds_max = ZERO, gridBounds_min = ZERO, &
                                                                  meshBounds_max = ZERO, meshBounds_min = ZERO
    integer(shortInt), dimension(3)                            :: n_xyz = 0
    type(cartesianCellSingle), dimension(:, :, :), allocatable :: grid
  contains

    ! Build procedures.
    procedure :: init
    procedure :: constructMapping
    procedure :: constructMapping_new
    procedure :: sortAngles
    procedure :: setGridIsOutsideMesh
    procedure :: gridFinitePrecision
    ! Runtime procedures.
    procedure :: getCellIsOutside
    procedure :: getGridBounds_min
    procedure :: getSpacing
    procedure :: getSpacingReciprocal
    procedure :: getGridWStar
    procedure :: getGridPhiCapital
    procedure :: getGridPhi
    procedure :: getGridChi
    procedure :: getGridIsOutsideBounds
    ! Analysis procedures
    procedure :: vertexValenceDistribution
    procedure :: edgeValenceDistribution
  end type cartesianGridSingle

contains

  !!
  !!
  !!
  subroutine init(self, vertices, edges, faces, elements)
    class(cartesianGridSingle), intent(inout)           :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(inout)                     :: edges
    class(faceShelf), intent(inout)                     :: faces
    class(elementShelf), intent(in)                     :: elements
    real(defReal), dimension(6)                         :: extremalCoordinates
    integer(shortInt)                                   :: i
    integer(shortInt), dimension(:), allocatable        :: currEdgeVertexIdxs
    real(defReal), dimension(3)                         :: currEdgeVector, xyz_max, xyz_min!, &
                                                           !centroid !!!!!
    real(defReal)                                       :: currEdgeLength, halfSpacingLessExtraRoom, sinAlpha, sinHalfAlpha

    !-----------------------------------------------------------------------------------------
    ! calculate constants for each face and assign them.
    ! (constant = dot(any point on the plane ⊥ the face, face normal))
    ! (needs to be changed) (set this value in other place and intent(in) not intent(inout))
    !-----------------------------------------------------------------------------------------
    do i = 1, faces % getSize()
      call faces % setFaceConst(i, -dot_product(faces % getFaceNormal(i), faces % getFaceCentroid(i)))

    end do

    !-----------------------------------------------------------------------------------------
    ! set l_min. Concurrently, set edgeLength and edgeUnitVector for all edges
    ! (needs to be changed) (set this value in other place and intent(in) not intent(inout))
    !-----------------------------------------------------------------------------------------
    self % l_min = INF

    do i = 1, edges % getSize()
      currEdgeVertexIdxs = edges % getEdgeVertexIdxs(i)
      currEdgeVector = vertices % getVertexCoordinates(currEdgeVertexIdxs(2)) &
                       - vertices % getVertexCoordinates(currEdgeVertexIdxs(1))
      currEdgeLength = norm2(currEdgeVector)

      call edges % setEdgeUnitVector(i, currEdgeVector / currEdgeLength)
      call edges % setEdgeLength(i, currEdgeLength)
      self % l_min = min(self % l_min, currEdgeLength)

    end do

    ! (needs to be changed) (written for temp operation; can be optimised further; calculate it during l_min calc.?)
    ! (can probably be a separte subroutine on its own, and the return value can be assigned as a grid attribute?)
    !print*, "Average edge length", calculateAvgEdgeLength(edges) 

    !-----------------------------------------------------------------------------------------
    ! set alpha
    !-----------------------------------------------------------------------------------------
    self % alpha = acos(max(findMinFaceAngle(edges, faces), findMinDihedralAngle(edges, faces, elements)))

    !-----------------------------------------------------------------------------------------
    ! set grid dimensions
    !-----------------------------------------------------------------------------------------
    ! Pre-compute sin(alpha) and sin(alpha / 2) and set cartesian cell spacing.
    sinAlpha = sin(self % alpha)
    sinHalfAlpha = sin(HALF * self % alpha)
    self % wStar = self % l_min * min(HALF, sinAlpha) * (ONE - epsilon(ONE))
    self % spacing = TWO * self % wStar * sinAlpha * sinHalfAlpha / (sqrt(THREE) * (ONE + sinAlpha) * (ONE + sinHalfAlpha))
    self % spacingReciprocal = ONE / self % spacing

    ! set n_xyz
    ! set minimum and max xyz-coordinates of the cartesian grid
    extremalCoordinates = vertices % getExtremalCoordinates()
    self % meshBounds_min = extremalCoordinates(1:3)
    self % meshBounds_max = extremalCoordinates(4:6)

    !(needs to be changed)(change the number of spacing for extra room for diff layers)
    do i = 1, 3
      halfSpacingLessExtraRoom = HALF * (self % spacing - &
                                         mod(self % meshBounds_max(i) - self % meshBounds_min(i), self % spacing))
      self % gridBounds_min(i) = self % meshBounds_min(i) - halfSpacingLessExtraRoom
      self % gridBounds_max(i) = self % meshBounds_max(i) + halfSpacingLessExtraRoom
      self % n_xyz(i) = nint((self % gridBounds_max(i) - self % gridBounds_min(i)) * self % spacingReciprocal)
        
    end do

    ! print cartesian grid parameters and mesh quality
    print*, "----------------------------------------------------"
    print*, "/\/\ Cartesian grid parameters and mesh quality /\/\"
    print*, "Minimum angle            : ", self % alpha
    print*, "Minimum edge length      : ", self % l_min
    print*, "No. of vertices          : ", vertices % getSize()
    print*, "No. of edges             : ", edges % getSize()
    print*, "No. of faces             : ", faces % getSize()
    print*, "No. of elements          : ", elements % getSize()
    print*, "Grid spacing             : ", self % spacing
    print*, "Grid size in x           : ", self % n_xyz(1)
    print*, "Grid size in y           : ", self % n_xyz(2)
    print*, "Grid size in z           : ", self % n_xyz(3)
    print*, "Grid lower bounds in xyz : ", self % gridBounds_min
    print*, "Grid upper bounds in xyz : ", self % gridBounds_max
    print*, "----------------------------------------------------"

    ! allocate grid matrix
    allocate(self % grid(self % n_xyz(1), self % n_xyz(2), self % n_xyz(3)))

    !-----------------------------------------------------------------------------------------
    !initialise for patch search
    !-----------------------------------------------------------------------------------------
    call self % constructMapping_new(edges, elements, faces, vertices)
    call self % sortAngles(edges, elements, faces, vertices)
    !call self % gridFinitePrecision(vertices, faces, elements) 
    !call self % setGridIsOutsideMesh()

    print*, "average edge length", calculateAvgEdgeLength(edges)

  end subroutine init

  !!
  !!
  !!
  subroutine constructMapping(self, vertices, edges, faces, elements)
    class(cartesianGridSingle), intent(inout)    :: self
    class(vertexShelf), intent(in)               :: vertices
    class(edgeShelf), intent(inout)              :: edges
    class(faceShelf), intent(in)                 :: faces
    class(elementShelf), intent(in)              :: elements
    integer(shortInt)                            :: i, j, k, l, nElementFaces
    integer(shortInt), dimension(:), allocatable :: currVertexIdxs, currElementFaceIdxs, currFaceEdgeIdxs
    integer(shortInt), dimension(6)              :: AABBIndices
    real(defReal)                                :: circumscribedBallRadius, targetDistance, extraDistance
    real(defReal), dimension(3)                  :: centroid, currEdgeVector, currFaceNormal
    real(defReal), dimension(:, :), allocatable  :: faceNormalSigns

    !----------------------------------------------------------------------------------------------
    ! edge interesection tests
    !----------------------------------------------------------------------------------------------
    !calculate constants for edge interesection tests
    circumscribedBallRadius = HALF * sqrt(THREE) * self % spacing
    targetDistance = self % wStar / (ONE + sin(self % alpha))

    print *, '--- BEGIN ---'

    do i = 1, edges % getSize()
      ! construct box for candidate cells
      currVertexIdxs = edges % getEdgeVertexIdxs(i)
      AABBIndices = constructAABB(vertices, currVertexIdxs, self % gridBounds_min, self % spacing)

      ! calculate edge-only-dependent properties
      currEdgeVector = edges % getEdgeUnitvector(i) * edges % getEdgeLength(i)

      !Loop over all cartesian cells in the box and test if each cell intersects with the edge
      !(needs to be changed) (k and l can be a function of j e.g. k = datum + slope*j so that box is narrowed down)
      do j = AABBIndices(1), AABBIndices(4)
        do k = AABBIndices(2), AABBIndices(5)
          do l = AABBIndices(3), AABBIndices(6)
            ! Store centroid info.
            centroid = self % gridBounds_min + self % spacing * ([real(j, defReal), real(k, defReal), real(l, defReal)] - HALF)
            call self % grid(j, k, l) % cellTestEdgeIntersection(vertices, edges, i, circumscribedBallRadius, targetDistance, &
                                                                 centroid, currEdgeVector, currVertexIdxs, &
                                                                 dot_product(currEdgeVector, currEdgeVector))
              
          end do 

        end do  

      end do

    end do

    print *, '--- FINISHED EDGES ---'

    !----------------------------------------------------------------------------------------------
    ! polyhedron inclusion tests
    !----------------------------------------------------------------------------------------------
    do i = 1, elements % getSize()
      ! construct box for candidate cells
      currVertexIdxs = elements % getElementVertexIdxs(i)
      AABBIndices = constructAABB(vertices, currVertexIdxs, self % gridBounds_min, self % spacing)

      ! calculate element-only-dependent properties
      currElementFaceIdxs = elements % getElementFaceIdxs(i)
      if(allocated(faceNormalSigns)) deallocate(faceNormalSigns)
      nElementFaces = size(currElementFaceIdxs)
      allocate(faceNormalSigns(3, nElementFaces))
      do j = 1, nElementFaces
        faceNormalSigns(:, j) = HALF * sign(ONE, faces % getFaceNormal(currElementFaceIdxs(j))) * self % spacing

      end do


      !Loop over all cartesian cells in the box and test if each cell is entirely included in the polyhedron
      !(needs to be changed) (k and l can be a function of j e.g. k = datum + slope*j so that box is narrowed down)
      do j = AABBIndices(1), AABBIndices(4)
        do k = AABBIndices(2), AABBIndices(5)
          do l = AABBIndices(3), AABBIndices(6)
            ! Store centroid info.
            centroid = self % gridBounds_min + self % spacing * ([real(j, defReal), real(k, defReal), real(l, defReal)] - HALF)
            call self % grid(j, k, l) % cellTestPolyhedronInclusion(faces, currElementFaceIdxs, centroid, faceNormalSigns, i)

          end do

        end do

      end do

    end do

    print *, '--- FINISHED POLYHEDRA INCLUSIONS ---'

    !----------------------------------------------------------------------------------------------
    ! face intersection tests
    !----------------------------------------------------------------------------------------------
    targetDistance = targetDistance * targetDistance
    do i = 1, faces % getSize()
      ! construct box for candidate cells
      currVertexIdxs = faces % getFaceVertexIdxs(i)
      AABBIndices = constructAABB(vertices, currVertexIdxs, self % gridBounds_min, self % spacing)

      ! calculate face-only-dependent properties
      currFaceEdgeIdxs = faces % getFaceEdgeIdxs(i)
      currFaceNormal = faces % getFaceNormal(i)
      extraDistance = HALF * sum(abs(currFaceNormal)) * self % spacing

      !Loop over all cartesian cells in the box and test if each cell intersect with the current face
      !(needs to be changed) (k and l can be a function of j e.g. k = datum + slope*j so that box is narrowed down)
      do j = AABBIndices(1), AABBIndices(4)
        do k = AABBIndices(2), AABBIndices(5)
          do l = AABBIndices(3), AABBIndices(6)
            ! Store centroid info.
            centroid = self % gridBounds_min + self % spacing * ([real(j, defReal), real(k, defReal), real(l, defReal)] - HALF)
            call self % grid(j, k, l) % cellTestFaceIntersection(vertices, edges, faces, currVertexIdxs, extraDistance, &
                                                                 currFaceNormal, centroid, self % spacing, i, currFaceEdgeIdxs, &
                                                                 targetDistance)

          end do

        end do

      end do

    end do

    print *, '--- FINISHED FACES ---'

    ! for cells that intersect more than one face, mappings constructions are performed during face intersection test
    ! for those that intersect exactly one face, mappings constructions are performed here.

    ! loop over all cartesian cells and call relevant subroutine
    do i = 1, self % n_xyz(1)
      do j = 1, self % n_xyz(2)
        do k = 1, self % n_xyz(3)
          call self % grid(i, j, k) % cellConstructMapSingleFace(faces)

        end do

      end do

    end do

    print *, '--- FINISHED SINGLE ---'

  end subroutine constructMapping

  !!
  !!
  !!
  subroutine constructMapping_new(self, edges, elements, faces, vertices)
    class(cartesianGridSingle), intent(inout)    :: self
    type(edgeShelf), intent(inout)               :: edges
    type(elementShelf), intent(in)               :: elements
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    integer(shortInt)                            :: i, j, k, l, debugEdgeIdx
    integer(shortInt), dimension(6)              :: cellIdxs
    integer(shortInt), dimension(:), allocatable :: elementFaceIdxs, elementVertexIdxs, faceVertexIdxs
    real(defReal)                                :: extraDistance, targetDistance
    real(defReal), dimension(2)                  :: tempIdxs
    real(defReal), dimension(3)                  :: centroid, realIdxs
    type(faceSATData)                            :: cache

    ! Pre-compute constants.
    targetDistance = (self % wStar / (ONE + sin(self % alpha))) ** 2

    ! Loop over all faces in the mesh.
    do i = 1, faces % getSize()
      ! Generate SAT data cache for current face.
      cache = faces % computeFaceSATData(i, edges, vertices)
      extraDistance = HALF * sum(abs(cache % faceNormal)) * self % spacing

      ! Generate the indices of the cells contained in the current face's AABB.
      faceVertexIdxs = faces % getFaceVertexIdxs(i)
      cellIdxs = constructAABB(vertices, faceVertexIdxs, self % gridBounds_min, self % spacingReciprocal)

      ! Loop over all cells.
      do l = cellIdxs(3), cellIdxs(6)
        centroid(3) = self % gridBounds_min(3) + self % spacing * (l - HALF)
        do k = cellIdxs(2), cellIdxs(5)
          centroid(2) = self % gridBounds_min(2) + self % spacing * (k - HALF)
          do j = cellIdxs(1), cellIdxs(4)
            centroid(1) = self % gridBounds_min(1) + self % spacing * (j - HALF)
            ! Now test current cell for intersection with the current face.
            call self % grid(j, k, l) % cellTestFaceIntersection_new(i, self % spacing, extraDistance, targetDistance, centroid, &
                                                                     edges, cache, faces, vertices)

          end do

        end do

      end do

    end do

    ! Now loop over all elements in the mesh.
    do i = 1, elements % getSize()
      ! Generate the indices of the cells contained in the current face's AABB.
      elementFaceIdxs = elements % getElementFaceIdxs(i)
      elementVertexIdxs = elements % getElementVertexIdxs(i)
      cellIdxs = constructAABB(vertices, elementVertexIdxs, self % gridBounds_min, self % spacingReciprocal)

      ! Loop over all cells.
      do l = cellIdxs(3), cellIdxs(6)
        tempIdxs(2) = l - HALF
        do k = cellIdxs(2), cellIdxs(5)
          tempIdxs(1) = k - HALF
          do j = cellIdxs(1), cellIdxs(4)
            if(self % grid(j, k, l) % getCellToElementIdx() == 0 .and. all(self % grid(j, k, l) % getFaceIdxs() == 0)) then
              centroid(1) = self % gridBounds_min(1) + self % spacing * (j - HALF)
              centroid(2) = self % gridBounds_min(2) + self % spacing * tempIdxs(1)
              centroid(3) = self % gridBounds_min(3) + self % spacing * tempIdxs(2)
              call self % grid(j, k, l) % cellTestPolyhedronInclusion_new(i, elementFaceIdxs, centroid, elements, faces)

            end if

          end do

        end do

      end do

    end do

  end subroutine constructMapping_new

  !!
  !!
  !!
  subroutine sortAngles(self, edges, elements, faces, vertices)
    class(cartesianGridSingle), intent(inout)    :: self
    type(edgeShelf), intent(inout)               :: edges
    type(elementShelf), intent(in)               :: elements
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    integer(shortInt)                            :: i, idx, j, k, l, pointerIdx, outer2LoopSize, nAngularSectors, nElements
    real(defReal), dimension(2)                  :: temp
    real(defReal), dimension(3)                  :: currEdgeUnitVector, localBasis1, localBasis2, currUnitVector
    integer(shortInt), dimension(:), allocatable :: currFaceEdgeIdxs, currEdgeElementIdxs, currEdgeFaceIdxs, faceIdxsArray, &
                                                    elementIdxsArray, sharingEdgeIdxs, currElementEdgeIdxs, commonEdgeIdxs, &
                                                    elementEdgeVertexIdxs, idxsArray, finalIdxsArray
    integer(shortInt), dimension(2)              :: currEdgeVertexIdxs, currVertexIdxs, face1ElementIdxs, face2ElementIdxs
    real(defReal)                                :: x, y, thetaHat
    real(defReal), dimension(:, :), allocatable  :: anglesArray, finalAnglesArray

    ! Loop through all edges.
    do i = 1, edges % getSize()
      !---------------------------------------------------------------------------------------------------------------
      ! construct 2D local coordinates system (localBasis1,localBasis2) on the plane whose normal is given as the 
      ! current edge's unit vector and contains the second vertex of the edge.
      !---------------------------------------------------------------------------------------------------------------
      currEdgeUnitVector = edges % getEdgeUnitVector(i)
      
      ! construct localBasis1
      if(abs(currEdgeUnitVector(1)) <= abs(currEdgeUnitVector(2)) .and. &
          abs(currEdgeUnitVector(1)) <= abs(currEdgeUnitVector(3))) then
        localBasis1 = [ZERO, currEdgeUnitVector(3), -currEdgeUnitVector(2)]

      elseif(abs(currEdgeUnitVector(2)) <= abs(currEdgeUnitVector(3))) then 
        localBasis1 = [-currEdgeUnitVector(3), ZERO, currEdgeUnitVector(1)]

      else
        localBasis1 = [currEdgeUnitVector(2), -currEdgeUnitVector(1), ZERO]

      end if
      localBasis1 = localBasis1 / norm2(localBasis1)

      ! construct localBasis2 (cross product gives the normalised vector)
      localBasis2 = crossProduct(currEdgeUnitVector, localBasis1)

      ! store localBasis1 and localBasis2 to each associated edge
      call edges % setEdgeLocalBasis1(i, localBasis1)
      call edges % setEdgeLocalBasis2(i, localBasis2)

      !---------------------------------------------------------------------------------------------------------------
      ! construct (unsorted) arrays for angles and associated elementIdxs
      !---------------------------------------------------------------------------------------------------------------
      ! retrieve relevant information
      currEdgeElementIdxs = edges % getEdgeElementIdxs(i)
      currEdgeFaceIdxs = edges % getEdgeFaceIdxs(i)
      currEdgeVertexIdxs = edges % getEdgeVertexIdxs(i)
      nElements = size(currEdgeElementIdxs)

      ! Retrieve all edges connected to the second vertex of the current edge.
      sharingEdgeIdxs = vertices % getVertexEdgeIdxs(currEdgeVertexIdxs(2))

      ! initialise arrays for angle and face index
      if(allocated(anglesArray)) deallocate(anglesArray)
      if(allocated(faceIdxsArray)) deallocate(faceIdxsArray)
      if(allocated(elementIdxsArray)) deallocate(elementIdxsArray)
      allocate(anglesArray(nElements, 2), faceIdxsArray(size(currEdgeFaceIdxs)), elementIdxsArray(nElements))

      ! Loop through all the elements containing the current edge.
      nAngularSectors = 0
      do j = 1, nElements
        ! Retrieve edges in the current element then find common edges with those sharing the second vertex of the current edge.
        currElementEdgeIdxs = elements % getElementEdgeIdxs(currEdgeElementIdxs(j))
        commonEdgeIdxs = findCommon(sharingEdgeIdxs, currElementEdgeIdxs)

        ! Reset idx = 0 then loop through all common edges.
        idx = 0
        do k = 1, size(commonEdgeIdxs)
          ! One of these edges is the original edge itself so skip it.
          if(i == commonEdgeIdxs(k)) cycle
          idx = idx + 1

          ! Get the vertices in the current edge then compute pseudo-angle.
          elementEdgeVertexIdxs = edges % getEdgeVertexIdxs(commonEdgeIdxs(k))
          currUnitVector = edges % getEdgeUnitVector(commonEdgeIdxs(k))
          if(elementEdgeVertexIdxs(1) /= currEdgeVertexIdxs(2)) currUnitVector = -currUnitVector
          x = dot_product(currUnitVector, localBasis1)
          y = dot_product(currUnitVector, localBasis2)
          temp(idx) = sign(ONE - x / (abs(x) + abs(y)), y)

        end do
        anglesArray(j, :) = [minval(temp), maxval(temp)]
        nAngularSectors = nAngularSectors + merge(1, 2, abs(anglesArray(j, 2) - anglesArray(j, 1)) <= TWO)

      end do

      ! Now split angular intervals which are outside the intervals [-2, 0] or [0, 2].
      idxsArray = [(j, j = 1, nAngularSectors)]
      if(allocated(finalAnglesArray)) deallocate(finalAnglesArray)
      if(allocated(finalIdxsArray)) deallocate(finalIdxsArray)
      allocate(finalAnglesArray(nAngularSectors, 2), finalIdxsArray(nAngularSectors))
      idx = 0
      do j = 1, nElements
        if(abs(anglesArray(j, 2) - anglesArray(j, 1)) <= TWO) then
          idx = idx + 1
          finalAnglesArray(idx, :) = anglesArray(j, :)
          finalIdxsArray(idx) = currEdgeElementIdxs(j)

        else
          idx = idx + 1
          finalAnglesArray(idx, :) = [-TWO, anglesArray(j, 1)]
          finalIdxsArray(idx) = currEdgeElementIdxs(j)

          idx = idx + 1
          finalAnglesArray(idx, :) = [anglesArray(j, 2), TWO]
          finalIdxsArray(idx) = currEdgeElementIdxs(j)

        end if

      end do

      ! pass and set elementIdxsArray and anglesArray to each corresponding edge
      call edges % setEdgeAnglesArray(i, finalAnglesArray)
      call edges % setEdgeElementIdxsArray(i, finalIdxsArray)

    end do

  end subroutine sortAngles

  !!
  !!
  !! 
  subroutine setGridIsOutsideMesh(self)
    class(cartesianGridSingle), intent(inout)              :: self
    integer(shortInt)                                      :: i, j, k
    !real(defReal), dimension(3)                            :: centroid !!!

    ! loop over all cartesian cells and call relevant subroutine
    do i = 1, self % n_xyz(1)
      do j = 1, self % n_xyz(2)
        do k = 1, self % n_xyz(3)
          
          !!!!!
          ! centroid(1) = (self % gridBounds_min(1)) + (self % spacing) * (i-0.5)
          ! centroid(2) = (self % gridBounds_min(2)) + (self % spacing) * (j-0.5)
          ! centroid(3) = (self % gridBounds_min(3)) + (self % spacing) * (k-0.5)
          !!!!!

          !!!!
          call self % grid(i,j,k) % setIsOutsideMesh()!centroid)

          ! ! if a cell intersects with neither any edge nor face, and it is not contained in a single polyhedron,
          ! ! then, this cell lies outside the computational domain for the unstructured mesh
          ! if (self % chi(i,j,k) == 0 .AND. self % phi(i,j,k) == 0 .AND. self % phiCapital(i,j,k) == 0) then
          !   ! if lies outside, set the chi value of the cell equal to -1. There are subroutines that test 
          !   ! if (chi != 0), but since this subroutine is called after all of those, they are unafftected.
          !   self % chi(i,j,k) = -1
          ! end if

          !!!!

        end do
      end do 
    end do

  end subroutine setGridIsOutsideMesh

  !!
  !!
  !!
  subroutine gridFinitePrecision(self, vertices, faces, elements)
    class(cartesianGridSingle), intent(inout)    :: self
    class(vertexShelf), intent(in)               :: vertices
    class(faceShelf), intent(in)                 :: faces
    class(elementShelf), intent(in)              :: elements
    integer(shortInt)                            :: i, j, k, l
    integer(shortInt), dimension(:), allocatable :: currVertexIdxs
    integer(shortInt), dimension(6)              :: AABBIndices
    real(defReal), dimension(3)                  :: centroid

    do i = 1, elements % getSize()
      ! construct box for candidate cells
      currVertexIdxs = elements % getElementVertexIdxs(i)
      AABBIndices = constructAABB(vertices, currVertexIdxs, self % gridBounds_min, self % spacing)

      ! calculate element-only-dependent properties (It is extremly rare that the code has to test this finitePrecision.
      ! Hence, we do not pre-calculate these unlike testPolyhedronInclusion.)

      !Loop over all cartesian cells in the box and test if each cell is entirely included in the polyhedron
      !(needs to be changed) (k and l can be a function of j e.g. k = datum + slope*j so that box is narrowed down)
      do j = AABBIndices(1), AABBIndices(4)
        do k = AABBIndices(2), AABBIndices(5)
          do l = AABBIndices(3), AABBIndices(6)
              ! Compute centroid.
              centroid = self % gridBounds_min + self % spacing * ([real(j, defReal), real(k, defReal), real(l, defReal)] - HALF)
              call self % grid(j, k, l) % cellFinitePrecision(faces, elements, centroid, i)

          end do

        end do

      end do

    end do

  end subroutine gridFinitePrecision

  !!
  !!
  !!
  pure function getCellIsOutside(self, idxs) result(isOutside)
    class(cartesianGridSingle), intent(in)      :: self
    integer(shortInt), dimension(3), intent(in) :: idxs
    logical(defBool)                            :: isOutside

    isOutside = self % grid(idxs(1), idxs(2), idxs(3)) % getCellIsOutside()

  end function getCellIsOutside

  !!
  !!
  !!
  pure function getGridBounds_min(self) result(gridBounds_min)
    class(cartesianGridSingle), intent(in)              :: self
    real(defReal), dimension(3)                         :: gridBounds_min

    gridBounds_min = self % gridBounds_min

  end function getGridBounds_min

  !!
  !!
  !!
  elemental function getSpacing(self) result(spacing)
    class(cartesianGridSingle), intent(in) :: self
    real(defReal)                          :: spacing

    spacing = self % spacing

  end function getSpacing

  !!
  !!
  !!
  elemental function getSpacingReciprocal(self) result(spacingReciprocal)
    class(cartesianGridSingle), intent(in) :: self
    real(defReal)                          :: spacingReciprocal

    spacingReciprocal = self % spacingReciprocal

  end function getSpacingReciprocal

  !!
  !!
  !!
  elemental function getGridWStar(self) result(wStar)
    class(cartesianGridSingle), intent(in)              :: self
    real(defReal)                                       :: wStar

    wStar = self % wStar

  end function getGridWStar

  !!
  !!
  !! 
  pure function getGridPhiCapital(self, cellIdxs) result(phiCapital)
    class(cartesianGridSingle), intent(in)      :: self
    integer(shortInt), dimension(3), intent(in) :: cellIdxs
    integer(shortInt)                           :: phiCapital

    phiCapital = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getCellToEdgeIdx()

  end function getGridPhiCapital

  !!
  !!
  !! 
  pure function getGridPhi(self, cellIdxs) result(phi)
    class(cartesianGridSingle), intent(in)      :: self
    integer(shortInt), dimension(3), intent(in) :: cellIdxs
    integer(shortInt)                           :: phi

    phi = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getCellToVertexIdx()

  end function getGridPhi

  !!
  !!
  !! 
  pure function getGridChi(self, cellIdxs) result(chi)
    class(cartesianGridSingle), intent(in)      :: self
    integer(shortInt), dimension(3), intent(in) :: cellIdxs
    integer(shortInt)                           :: chi

    chi = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getCellToElementIdx()

  end function getGridChi

  !!
  !!
  !! (needs to be changed) (possible acceleration?)
  pure function getGridIsOutsideBounds(self, r) result(isOutside)
    class(cartesianGridSingle), intent(in)              :: self
    real(defReal), dimension(3), intent(in)             :: r
    logical                                             :: isOutside

    isOutside = .true.
    if(any(r < self % meshBounds_min) .or. any(self % meshBounds_max < r)) return
    isOutside = .false.

  end function getGridIsOutsideBounds

  !!
  !!
  !!
  subroutine vertexValenceDistribution(self, vertices)
    class(cartesianGridSingle), intent(in)        :: self
    class(vertexShelf), intent(in)                :: vertices
    integer(shortInt)                             :: i, max, sizeArr
    integer(shortInt), dimension(:), allocatable  :: temp, distribution
    
    ! Find maximum valence
    max = 0
    do i = 1, vertices % getSize()
      temp = vertices % getVertexElementIdxs(i)
      if (max < size(temp)) then
        max = size(temp)
      end if
    end do

    ! Allocate distrubution array
    allocate(distribution(max))
    distribution = 0

    ! Find the distribution
    do i = 1, vertices % getSize()
      temp = vertices % getVertexElementIdxs(i)
      sizeArr = size(temp)
      distribution(sizeArr) = distribution(sizeArr) + 1
    end do

    print*, "vertexValence"
    print*, distribution

  end subroutine vertexValenceDistribution

  !!
  !!
  !!
  subroutine edgeValenceDistribution(self, edges)
    class(cartesianGridSingle), intent(in)        :: self
    class(edgeShelf), intent(inout)               :: edges
    integer(shortInt)                             :: i, max, sizeArr
    integer(shortInt), dimension(:), allocatable  :: temp, distribution
    
    ! Find maximum valence
    max = 0
    do i = 1, edges % getSize()
      temp = edges % getEdgeElementIdxs(i)
      if (max < size(temp)) then
        max = size(temp)
      end if
    end do

    ! Allocate distrubution array
    allocate(distribution(max))
    distribution = 0

    ! Find the distribution
    do i = 1, edges % getSize()
      temp = edges % getEdgeElementIdxs(i)
      sizeArr = size(temp)
      distribution(sizeArr) = distribution(sizeArr) + 1
    end do

    print*, "edgeValence"
    print*, distribution

  end subroutine edgeValenceDistribution

end module CartesianGridSingle_class