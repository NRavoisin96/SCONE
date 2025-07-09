module cartesianGrid_class
  
  use universalVariables,           only : ZERO, INF
  use vertexShelf_class,            only : vertexShelf
  use edgeShelf_class,              only : edgeShelf
  use faceShelf_class,              only : faceShelf
  use elementShelf_class,           only : elementShelf
  use genericProcedures,            only : findCommon, crossProduct
  use numPrecision                  
  use cartesianCell_class,          only : cartesianCell
  
  implicit none
  private
  
  !!
  !!
  type, public                                          :: cartesianGrid
    private
    real(defReal)                                       :: spacing = ZERO, alpha = ZERO, l_min = ZERO, &
                                                           wStar = ZERO, spacingReciprocal = ZERO
    real(defReal), dimension(3)                         :: gridBounds_max = ZERO, gridBounds_min = ZERO
    integer(shortInt), dimension(3)                     :: n_xyz = 0
    type(cartesianCell), dimension(:,:,:), allocatable  :: grid

  contains

    ! Build procedures.
    procedure                                    :: init
    procedure                                    :: findMinFaceAngle
    procedure                                    :: findMinDihedralAngle
    procedure                                    :: constructMapping
    procedure                                    :: sortAngles
    procedure                                    :: constructAABB
    procedure                                    :: sortPairs
    procedure                                    :: setGridIsOutside
    ! Runtime procedures.
    procedure                                    :: binarySearchAngle
    procedure                                    :: getGridBounds_min
    procedure                                    :: getSpacingReciprocal
    procedure                                    :: getGridWStar
    procedure                                    :: getGridPhiCapital
    procedure                                    :: getGridPhi
    procedure                                    :: getGridChi
    procedure                                    :: getIsOutside
  end type cartesianGrid

contains

  !!
  !!
  !!
  subroutine init(self, vertices, edges, faces, elements)
    class(cartesianGrid), intent(inout)           :: self
    class(vertexShelf), intent(in)                :: vertices
    class(edgeShelf), intent(inout)               :: edges
    class(faceShelf), intent(inout)               :: faces
    class(elementShelf), intent(in)               :: elements
    real(defReal), dimension(6)                   :: extremalCoordinates
    integer(shortInt)                             :: i!!!, temp, j, k
    integer(shortInt), dimension(:), allocatable  :: currEdgeVertexIdxs
    real(defReal), dimension(3)                   :: currEdgeVector, extraRoom, xyz_max, xyz_min
    real(defReal)                                 :: maxCosValue, tempMaxCosValue, currEdgeLength

    !-----------------------------------------------------------------------------------------
    ! calculate constants for each face and assign them.
    ! (constant = dot(any point on the plane ⊥ the face, face normal))
    ! (needs to be changed) (set this value in other place and intent(in) not intent(inout))
    !-----------------------------------------------------------------------------------------
    do i = 1, faces % getSize()
      call faces % setFaceConst(i, dot_product(faces % getFaceNormal(i), faces % getFaceCentroid(i))*(-1))
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

      call edges % setEdgeUnitVector(i, currEdgeVector/currEdgeLength)
      call edges % setEdgeLength(i, currEdgeLength)

      if (currEdgeLength < self % l_min) self % l_min = currEdgeLength

    end do
    !-----------------------------------------------------------------------------------------
    ! set alpha
    !-----------------------------------------------------------------------------------------
    maxCosValue = self % findMinFaceAngle(edges, faces)

    tempMaxCosValue = self % findMinDihedralAngle(edges, faces, elements)
    if (maxCosValue < tempMaxCosValue) maxCosValue = tempMaxCosValue

    self % alpha = ACOS(maxCosValue)

    !-----------------------------------------------------------------------------------------
    ! set grid dimensions
    !-----------------------------------------------------------------------------------------
    ! set cartesian cell spacing (needs to be changed) (times by 0.9999 for wStar?)
    self % wStar = (self % l_min)*min(0.5d0, SIN(self % alpha))
    self % spacing = 2*(self % wStar)*SIN(self % alpha)*SIN((self % alpha)/2)&
                     /sqrt(3.0d0)/(1+SIN(self % alpha))/(1+SIN((self % alpha)/2))
    self % spacingReciprocal = 1/(self % spacing)

    ! set n_xyz
    ! set minimum and max xyz-coordinates of the cartesian grid
    extremalCoordinates = vertices % getExtremalCoordinates()
    xyz_min = extremalCoordinates(1:3)
    xyz_max = extremalCoordinates(4:6)

    !(needs to be changed)(change the number of spacing for extra room for diff layers)
    do i = 1, 3
        extraRoom(i) = mod(xyz_max(i) - xyz_min(i), self % spacing)
        self % gridBounds_min(i) = xyz_min(i) - (self % spacing - extraRoom(i))/2
        self % gridBounds_max(i) = xyz_max(i) + (self % spacing - extraRoom(i))/2

        self % n_xyz(i) = NINT((self % gridBounds_max(i) - self % gridBounds_min(i))/(self % spacing))
        
    end do

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
    call self % constructMapping(vertices, edges, faces, elements)
    call self % sortAngles(edges, faces, vertices) !!! vertices
    call self % setGridIsOutside()

    ! !!!
    ! do i = 1, self % n_xyz(1)
    !   do j = 1, self % n_xyz(2)
    !     do k = 1, self % n_xyz(3)
    !       temp = self % grid(i,j,k) % getChi()
    !       temp =  self % grid(i,j,k) % getPhi()
    !       temp =  self % grid(i,j,k) % getPhiCapital()
    !     end do
    !   end do 
    ! end do
    ! !!!

  end subroutine init

  !!
  !!
  !!
  function findMinFaceAngle(self, edges, faces) result(maxCosValue)
    class(cartesianGrid), intent(inout)           :: self
    class(edgeShelf), intent(in)                  :: edges
    class(faceShelf), intent(in)                  :: faces
    integer(shortInt)                             :: i, j, k, l ,m, sign1, sign2
    real(defReal)                                 :: currCosValue, maxCosValue
    integer(shortInt), dimension(:), allocatable  :: currFaceEdgeIdxs
    integer(shortInt), dimension(2)               :: currEdge1VertexIdxs, currEdge2VertexIdxs
    !integer(shortInt) :: temp

    !temp = 0
    maxCosValue = -1.0d0

    do i = 1, faces % getSize()
      currFaceEdgeIdxs = faces % getFaceEdgeIdxs(i)

      do j = 1, size(currFaceEdgeIdxs) - 1
        currEdge1VertexIdxs = edges % getEdgeVertexIdxs(currFaceEdgeIdxs(j))

        do k = j + 1, size(currFaceEdgeIdxs)
          currEdge2VertexIdxs = edges % getEdgeVertexIdxs(currFaceEdgeIdxs(k))

          do l = 1, 2

            do m = 1, 2
              if (currEdge1VertexIdxs(l) == currEdge2VertexIdxs(m)) then

                ! correct the direction of unit vector of each edge
                if (l == 1) then
                  sign1 = 1
                else 
                  sign1 = -1
                end if

                if (m == 1) then
                  sign2 = 1
                else 
                  sign2 = -1
                end if

                ! calculate cosine value
                currCosValue = dot_product(edges % getEdgeUnitvector(currFaceEdgeIdxs(j)), &
                                           edges % getEdgeUnitvector(currFaceEdgeIdxs(k)))*sign1*sign2

                !if (maxCosValue < currCosValue) temp = i
                ! update maxCosValue
                if (maxCosValue < currCosValue) maxCosValue = currCosValue

              end if

            end do

          end do

        end do

      end do

    end do

    !print*, temp

  end function findMinFaceAngle

  !!
  !!
  !!
  function findMinDihedralAngle(self, edges, faces, elements) result(maxCosValue)
    class(cartesianGrid), intent(inout)           :: self
    class(edgeShelf), intent(in)                  :: edges
    class(faceShelf), intent(in)                  :: faces
    class(elementShelf), intent(in)               :: elements
    integer(shortInt)                             :: i, j, k
    real(defReal)                                 :: currCosValue, maxCosValue
    integer(shortInt), dimension(:), allocatable  :: currElementFaceIdxs, currElementEdgeIdxs
    integer(shortInt), dimension(2)               :: candidateFaceIdxs, candidateElementIdxs, signArray

    maxCosValue = -1.0d0

    do i = 1, elements % getSize()
      currElementFaceIdxs = abs(elements % getElementFaceIdxs(i))
      currElementEdgeIdxs = elements % getElementEdgeIdxs(i)

      do j = 1, size(currElementEdgeIdxs)
        candidateFaceIdxs = findCommon(currElementFaceIdxs, edges % getEdgeFaceIdxs(currElementEdgeIdxs(j)))

        do k = 1, 2
          candidateElementIdxs = faces % getFaceElementIdxs(candidateFaceIdxs(k))

          if (i < candidateElementIdxs(1) .OR. i < candidateElementIdxs(2)) then
            signArray(k) = 1
          else
            signArray(k) = -1
          end if

        end do

        currCosValue = dot_product(faces % getFaceNormal(candidateFaceIdxs(1)), &
                                   faces % getFaceNormal(candidateFaceIdxs(2)))*signArray(1)*signArray(2)*(-1)

        ! update maxCosValue
        if (maxCosValue < currCosValue) maxCosValue = currCosValue
        
      end do

    end do


  end function findMinDihedralAngle

  !!
  !!
  !!
  subroutine constructMapping(self, vertices, edges, faces, elements)
    class(cartesianGrid), intent(inout)             :: self
    class(vertexShelf), intent(in)                  :: vertices
    class(edgeShelf), intent(inout)                 :: edges
    class(faceShelf), intent(in)                    :: faces
    class(elementShelf), intent(in)                 :: elements
    integer(shortInt)                               :: i, j, k, l
    integer(shortInt), dimension(:), allocatable    :: currVertexIdxs, currElementFaceIdxs, currFaceEdgeIdxs
    integer(shortInt), dimension(6)                 :: AABBIndices
    real(defReal)                                   :: circumscribedBallRadius, targetDistance, a, &
                                                       extraDistance, cellSpacing
    real(defReal), dimension(3)                     :: centroid, currEdgeVector, currFaceNormal
    real(defReal), dimension(:,:), allocatable      :: faceNormalSigns
    
    cellSpacing = self % spacing

    !----------------------------------------------------------------------------------------------
    ! edge interesection tests
    !----------------------------------------------------------------------------------------------
    !calculate constants for edge interesection tests
    circumscribedBallRadius = sqrt(3.0d0)*(self % spacing)/2
    targetDistance = (self % wStar) / (1 + SIN(self % alpha))

    do i = 1, edges % getSize()

        ! construct box for candidate cells
        currVertexIdxs = edges % getEdgeVertexIdxs(i)
        AABBIndices = self % constructAABB(vertices, currVertexIdxs)

        ! calculate edge-only-dependent properties
        currEdgeVector = (edges % getEdgeUnitvector(i))*(edges % getEdgeLength(i))
        a = dot_product(currEdgeVector, currEdgeVector)

        !Loop over all cartesian cells in the box and test if each cell intersects with the edge
        !(needs to be changed) (k and l can be a function of j e.g. k = datum + slope*j so that box is narrowed down)
        do j = AABBIndices(1), AABBIndices(4)
            do k = AABBIndices(2), AABBIndices(5)
                do l = AABBIndices(3), AABBIndices(6)

                    ! (needs to be changed) (store centroid info)
                    centroid(1) = (self % gridBounds_min(1)) + (self % spacing) * (j-0.5)
                    centroid(2) = (self % gridBounds_min(2)) + (self % spacing) * (k-0.5)
                    centroid(3) = (self % gridBounds_min(3)) + (self % spacing) * (l-0.5)

                    call self % grid(j,k,l) % testEdgeIntersection(vertices, edges, i, circumscribedBallRadius, &
                                                    targetDistance, centroid, currEdgeVector, currVertexIdxs, a)
                  
                end do 
            end do    
        end do

    end do

    !----------------------------------------------------------------------------------------------
    ! polyhedron inclusion tests
    !----------------------------------------------------------------------------------------------
    allocate(faceNormalSigns(3, 2))
    do i = 1, elements % getSize()

        ! construct box for candidate cells
        currVertexIdxs = elements % getElementVertexIdxs(i)
        AABBIndices = self % constructAABB(vertices, currVertexIdxs)

        ! calculate element-only-dependent properties
        currElementFaceIdxs = elements % getElementFaceIdxs(i)

        deallocate(faceNormalSigns)
        allocate(faceNormalSigns(3, size(currElementFaceIdxs)))
        do j = 1, size(currElementFaceIdxs)
          currFaceNormal = faces % getFaceNormal(currElementFaceIdxs(j))
          do k = 1, 3
            if (currFaceNormal(k) > 0) then
              faceNormalSigns(k, j) = 1
            else
              faceNormalSigns(k, j) = -1
            end if
          end do
        end do
        faceNormalSigns = faceNormalSigns * (self % spacing)/2


        !Loop over all cartesian cells in the box and test if each cell is entirely included in the polyhedron
        !(needs to be changed) (k and l can be a function of j e.g. k = datum + slope*j so that box is narrowed down)
        do j = AABBIndices(1), AABBIndices(4)
            do k = AABBIndices(2), AABBIndices(5)
                do l = AABBIndices(3), AABBIndices(6)

                    ! (needs to be changed) (store centroid info)
                    centroid(1) = (self % gridBounds_min(1)) + (self % spacing) * (j-0.5)
                    centroid(2) = (self % gridBounds_min(2)) + (self % spacing) * (k-0.5)
                    centroid(3) = (self % gridBounds_min(3)) + (self % spacing) * (l-0.5)

                    call self % grid(j,k,l) % testPolyhedronInclusion(faces, currElementFaceIdxs, centroid, &
                                                                      faceNormalSigns, i)
                  
                end do 
            end do    
        end do

    end do

    !----------------------------------------------------------------------------------------------
    ! face intersection tests
    !----------------------------------------------------------------------------------------------
    targetDistance = targetDistance**2
    do i = 1, faces % getSize()

      ! construct box for candidate cells
      currVertexIdxs = faces % getFaceVertexIdxs(i)
      AABBIndices = self % constructAABB(vertices, currVertexIdxs)

      ! calculate face-only-dependent properties
      currFaceEdgeIdxs = faces % getFaceEdgeIdxs(i)
      currFaceNormal = faces % getFaceNormal(i)
      extraDistance = (abs(currFaceNormal(1)) + abs(currFaceNormal(2)) + abs(currFaceNormal(3))) &
                      * (self%spacing) / 2               

      !Loop over all cartesian cells in the box and test if each cell intersect with the current face
      !(needs to be changed) (k and l can be a function of j e.g. k = datum + slope*j so that box is narrowed down)
      do j = AABBIndices(1), AABBIndices(4)
          do k = AABBIndices(2), AABBIndices(5)
              do l = AABBIndices(3), AABBIndices(6)

                  ! (needs to be changed) (store centroid info)
                  centroid(1) = (self % gridBounds_min(1)) + (self % spacing) * (j-0.5)
                  centroid(2) = (self % gridBounds_min(2)) + (self % spacing) * (k-0.5)
                  centroid(3) = (self % gridBounds_min(3)) + (self % spacing) * (l-0.5)

                  call self % grid(j,k,l) % testFaceIntersection(vertices, edges, faces, &
                                            currVertexIdxs, extraDistance, currFaceNormal, &
                                            centroid, cellSpacing, i, currFaceEdgeIdxs, targetDistance)
                
              end do 
          end do    
      end do

    end do

    ! for cells that intersec more than one face, mappings constructions are performed during face intersection test
    ! for those that intersect exactly one face, mappings constructions are performed here.

    !!!
    print*, "beginning single Face case"
    !!!

    ! loop over all cartesian cells and call relevant subroutine
    do i = 1, self % n_xyz(1)
      do j = 1, self % n_xyz(2)
        do k = 1, self % n_xyz(3)
          call self % grid(i,j,k) % constructMapSingleFace(faces)
        end do
      end do 
    end do

    !!!
    print*, "ending single Face case"
    !!!

  end subroutine constructMapping

  !!
  !!
  !!
  subroutine sortAngles(self, edges, faces, vertices) !!! vertices
    class(cartesianGrid), intent(inout)           :: self
    class(vertexShelf), intent(in)                  :: vertices !!!
    class(edgeShelf), intent(inout)               :: edges
    class(faceShelf), intent(in)                  :: faces
    integer(shortInt)                             :: i, j, k, l, m, n, currPhiCapital, v_e, pointerIdx, &
                                                     outer2LoopSize, currEdgeVertex1Idx
    real(defReal), dimension(3)                   :: currEdgeUnitVector, localBasis1, localBasis2, currUnitVector
    integer(shortInt), dimension(:), allocatable  :: currEdgeFaceIdxs, currFaceEdgeIdxs, faceIdxsArray, &
                                                     elementIdxsArray
    integer(shortInt), dimension(2)               :: currEdgeVertexIdxs, currVertexIdxs, face1ElementIdxs, &
                                                     face2ElementIdxs
    real(defReal)                                 :: x, y, thetaHat
    real(defReal), dimension(:), allocatable      :: anglesArray

    !!!
    !logical                                       :: flagBoundaryFace
    !!!

    !!!
    ! integer(shortInt)                             :: temp
    ! temp = 0
    ! print*, "beginning sorting angles"
    !!!

    ! loop through all cartesian cells so that only [edge index]s, where there exists at least one [cell index] s.t. 
    ! phiCapital([cell index]) = [edge index], are used.
    ! (needs to be changed) (since all edges (most likely) are assigned for phiCapital mapping anyways, just loop through all edges?)
    do i = 1, self % n_xyz(1)
      do j = 1, self % n_xyz(2)
        do k = 1, self % n_xyz(3)

          currPhiCapital = self % grid(i,j,k) % getPhiCapital()

          ! continue only if the current cell contains valid edge index mapping of phiCapital
          if (currPhiCapital /= 0) then
            ! continue only if the current edge index (= phiCapital) has not been used for sorting angles yet
            if (.NOT. edges % isAllocatedEdgeAnglesArray(currPhiCapital)) then

              !!!
              ! temp = temp + 1
              !!!

              !---------------------------------------------------------------------------------------------------------------
              ! construct 2D local cooridnate system (localBasis1,localBasis2) on the plane whose normal is given as the current edge's unit vector
              ! and contains the second vertex of the edge.
              !---------------------------------------------------------------------------------------------------------------
              currEdgeUnitVector = edges % getEdgeUnitVector(currPhiCapital)
              
              ! construct localBasis1
              if (abs(currEdgeUnitVector(1)) <= abs(currEdgeUnitVector(2)) .AND. &
                  abs(currEdgeUnitVector(1)) <= abs(currEdgeUnitVector(3))) then
                    localBasis1 = [0.0d0, currEdgeUnitVector(3), -currEdgeUnitVector(2)]
              elseif (abs(currEdgeUnitVector(2)) <= abs(currEdgeUnitVector(3))) then 
                    localBasis1 = [-currEdgeUnitVector(3), 0.0d0, currEdgeUnitVector(1)]
              else
                    localBasis1 = [currEdgeUnitVector(2), -currEdgeUnitVector(1), 0.0d0]
              end if
              localBasis1 = localBasis1 / norm2(localBasis1)

              ! construct localBasis2 (cross product gives the normalised vector)
              localBasis2 = crossProduct(currEdgeUnitVector, localBasis1)

              ! store localBasis1 and localBasis2 to each associated edge
              call edges % setEdgeLocalBasis1(currPhiCapital, localBasis1)
              call edges % setEdgeLocalBasis2(currPhiCapital, localBasis2)

              !---------------------------------------------------------------------------------------------------------------
              ! construct (unsorted) arrays for angles and associated elementIdxs
              !---------------------------------------------------------------------------------------------------------------
              ! retrieve relevant information
              currEdgeFaceIdxs = edges % getEdgeFaceIdxs(currPhiCapital)
              currEdgeVertexIdxs = edges % getEdgeVertexIdxs(currPhiCapital)
              v_e = currEdgeVertexIdxs(2)
              currEdgeVertex1Idx = currEdgeVertexIdxs(1)

              ! initialise arrays for angle and face index
              if (allocated(anglesArray)) deallocate(anglesArray)
              if (allocated(faceIdxsArray)) deallocate(faceIdxsArray)
              if (allocated(elementIdxsArray)) deallocate(elementIdxsArray)
              allocate(anglesArray(size(currEdgeFaceIdxs)))
              allocate(faceIdxsArray(size(currEdgeFaceIdxs)))
              allocate(elementIdxsArray(size(currEdgeFaceIdxs)))
              !elementIdxsArray(:) = 0

              ! loop through all faces attached to the current edge (= currPhiCapital)
              outer1: do l = 1, size(currEdgeFaceIdxs)
                currFaceEdgeIdxs = faces % getFaceEdgeIdxs(currEdgeFaceIdxs(l))
                
                ! loop through all edges of the current face (of the current edge = currPhiCapital)
                inner1: do m = 1, size(currFaceEdgeIdxs)
                  currVertexIdxs = edges % getEdgeVertexIdxs(currFaceEdgeIdxs(m))

                  ! test if this current edge of the current face is the one that is connected to edge = currPhiCapital
                  ! and that this current edge is not currPhiCapital itself
                  if (ANY(currVertexIdxs == v_e)) then
                    if (.NOT. ANY(currVertexIdxs == currEdgeVertex1Idx)) then

                      ! retrieve and correct the orientation of the unit vector of the edge
                      if (currVertexIdxs(1) == v_e) then
                        currUnitVector = edges % getEdgeUnitVector(currFaceEdgeIdxs(m))
                      else
                        currUnitVector = edges % getEdgeUnitVector(currFaceEdgeIdxs(m))*(-1)
                      end if

                      ! calculate the 2D local coordinates
                      x = dot_product(currUnitVector, localBasis1)
                      y = dot_product(currUnitVector, localBasis2)

                      ! calculate pseudo angle
                      thetaHat = SIGN(1 - (x / (abs(x) + abs(y))), y)

                      ! add the calculated angle and face index to each corresponding arrays
                      anglesArray(l) = thetaHat
                      faceIdxsArray(l) = currEdgeFaceIdxs(l)

                      ! once the calculation is performed for current face, move on to the next face of the current edge (= currPhiCapital)
                      ! (l = l + 1)
                      exit inner1

                    end if
                  end if

                end do inner1

              end do outer1

              ! perform index sorting on anglesArray and faceIdxsArray
              call self % sortPairs(anglesArray, faceIdxsArray)

              ! construct a sorted array for element index
              ! pointer of index for sorted array assignment
              pointerIdx = 1

              !!!
              ! initialise flag whether any face is a boundary face
              !flagBoundaryFace = .FALSE.
              !!!

              ! loop through all neighbouring faces
              outer2LoopSize = size(faceIdxsArray)
              outer2: do l = 1, outer2LoopSize
                face1ElementIdxs = faces % getFaceElementIdxs(faceIdxsArray(l))
                face2ElementIdxs = faces % getFaceElementIdxs(faceIdxsArray(mod(l, outer2LoopSize) + 1))

                !!!
                ! checking if any face is a boundary face
                !if (size(faces % getFaceElementIdxs(faceIdxsArray(l))) == 1) flagBoundaryFace = .TRUE.
                !!!

                ! if boundary face, set the second element index = 0, which can be used when constructing elementIdxsArray. 
                if (size(faces % getFaceElementIdxs(faceIdxsArray(l))) == 1) face1ElementIdxs(2) = 0
                if (size(faces % getFaceElementIdxs(faceIdxsArray(mod(l, outer2LoopSize) + 1))) == 1) face2ElementIdxs(2) = 0

                ! find the common element indicies from 2 X arrays of size 2
                ! (needs checking) (check and throw fatal error if there are two common element idxs-for boundary edge?)
                ! (There can be only one element attached to the face for boundary faces?)
                ! (If it is a boundary edge with only two faces attached, no need to calculate this)
                ! (During face intersection, if a cell intersects one or two boundary faces (of the same element))
                ! (then we can directly assign "chi" rather than phi and phiCapital)
                ! (and this removes the issue disscussed above because only [edge idx] s.t. there exists [cell idx] s.t.)
                ! (phiCapital([cell idx]) = [edgeIdx] are used for this subroutine)
                ! (this applies for boundary edges with any number of faces attached to it)
                ! (needs to be changed) (due to the possible acceleration above)
                middle2: do m = 1, 2
                  inner2: do n = 1, 2
                    if (face1ElementIdxs(m) == face2ElementIdxs(n)) then
                      elementIdxsArray(pointerIdx) = face1ElementIdxs(m)
                      pointerIdx = pointerIdx + 1
                      exit middle2
                    end if
                  end do inner2
                end do middle2

              end do outer2

              !!!
              ! if any of the faces is a boundary face, set the last element index = 0
              ! (needs to be changed) (it might not be the last element index in some cases?)
              ! (initially set all indices = 0 so that if face1ElementIdxs(m) != face2ElementIdxs(n))
              ! (for all m and n, space outside the unstructured mesh has element index = 0)
              ! (for this move "pointerIdx = pointerIdx + 1" out of loops) (have tried but not work?)
              ! if (flagBoundaryFace) elementIdxsArray(size(elementIdxsArray)) = 0
              ! solved problem by: if (size(faces % getFaceElementIdxs(faceIdxsArray(l))) == 1) face1ElementIdxs(2) = 0
              !!!

              ! pass and set elementIdxsArray and anglesArray to each corresponding edge
              call edges % setEdgeAnglesArray(currPhiCapital, anglesArray)
              call edges % setEdgeElementIdxsArray(currPhiCapital, elementIdxsArray)

              ! set if the current edge is a boundary edge (containing boundary face with a single element index)
              ! (needs to be changed) (can be accelerated further?)
              ! if (ANY(elementIdxsArray == 0)) then
              !   call edges % setEdgeIsBoundary(currPhiCapital)
              ! end if
              !!! not needed anymore

              ! outside of this subroutine
              !! construct: self % grid(i,j,k) % getCellPhiCapital()
              !! construct: edges % setEdge2dBasis(currPhiCapital, u, v)
              !setEdgeAnglesArray / setEdgeElementIdxsArray


              !!!
              ! print*, "-------------------------------------------------------------------"
              ! do l = 1, size(currEdgeFaceIdxs)
              !   currFaceEdgeIdxs = faces % getFaceVertexIdxs(currEdgeFaceIdxs(l))
              !   print*, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
              !   print*, "Face index", currEdgeFaceIdxs(l)
              !   do m = 1, size(currFaceEdgeIdxs)
              !     print*, vertices % getVertexCoordinates(currFaceEdgeIdxs(m))
              !   end do
              !   print*, faces % getFaceElementIdxs(currEdgeFaceIdxs(l))
              ! end do
              ! print*, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
              ! print*, localBasis1
              ! print*, localBasis2
              ! print*, anglesArray
              ! print*, faceIdxsArray
              ! print*, elementIdxsArray
              ! !print*, edges % getEdgeIsBoundary(currPhiCapital)
              ! print*, "-------------------------------------------------------------------"

              ! print*, edges % getEdgeAnglesArray(currPhiCapital)
              ! print*, edges % getEdgeElementIdxsArray(currPhiCapital)
              ! print*, currEdgeUnitVector
              ! print*, localBasis1
              ! print*, localBasis2
              ! print*, currEdgeFaceIdxs
              ! print*, currEdgeVertexIdxs
              !!!







            end if
          end if




        end do
      end do 
    end do

    !!!
    print*, "ending sorting angles"
    ! print*, temp
    ! print*, SIGN(1 - (COS(0.0) / (abs(COS(0.0)) + abs(SIN(0.0)))), SIN(0.0))
    ! print*, SIGN(1 - (COS(3.14/4) / (abs(COS(3.14/4)) + abs(SIN(3.14/4)))), SIN(3.14/4))
    ! print*, SIGN(1 - (COS(3.14/2) / (abs(COS(3.14/2)) + abs(SIN(3.14/2)))), SIN(3.14/4))
    ! print*, SIGN(1 - (COS(3*3.14/4) / (abs(COS(3*3.14/4)) + abs(SIN(3*3.14/4)))), SIN(3*3.14/4))
    ! print*, SIGN(1 - (COS(3.14) / (abs(COS(3.14)) + abs(SIN(3.14)))), SIN(3.14))
    ! print*, SIGN(1 - (COS(5*3.14/4) / (abs(COS(5*3.14/4)) + abs(SIN(5*3.14/4)))), SIN(5*3.14/4))
    ! print*, SIGN(1 - (COS(6*3.14/4) / (abs(COS(6*3.14/4)) + abs(SIN(6*3.14/4)))), SIN(6*3.14/4))
    ! print*, SIGN(1 - (COS(7*3.14/4) / (abs(COS(7*3.14/4)) + abs(SIN(7*3.14/4)))), SIN(7*3.14/4))
    !!!



  end subroutine sortAngles

  !!
  !!
  !! AABBIndices = [xmin, ymin, zmin, xmax, ymax, zmax]
  function constructAABB(self, vertices, currVertexIdxs) result(AABBIndices)
    class(cartesianGrid), intent(inout)              :: self
    class(vertexShelf), intent(in)                   :: vertices
    integer(shortInt), dimension(:), intent(in)      :: currVertexIdxs
    integer(shortInt), dimension(6)                  :: AABBIndices
    real(defreal), dimension(3)                      :: xyz_min, xyz_max, currVertexCoords
    integer(shortInt)                                :: i, j

    ! initialise xyz_min and xyz_max using the first vertex
    currVertexCoords = vertices % getVertexCoordinates(currVertexIdxs(1))
    xyz_max = currVertexCoords
    xyz_min = currVertexCoords

    ! find xyz_min and xyz_max 
    do i = 2, size(currVertexIdxs)
        currVertexCoords = vertices % getVertexCoordinates(currVertexIdxs(i))

        do j = 1, 3
            if (xyz_min(j) > currVertexCoords(j)) then
                 xyz_min(j) = currVertexCoords(j)
            elseif (xyz_max(j) < currVertexCoords(j)) then
                xyz_max(j) = currVertexCoords(j)
            end if
        end do

    end do

    ! find AABBIndices
    do i = 1, 3
        AABBIndices(i) = ceiling((xyz_min(i) - self % gridBounds_min(i))/(self % spacing)) 
        AABBIndices(3+i) = ceiling((xyz_max(i) - self % gridBounds_min(i))/(self % spacing)) 
    end do

  end function constructAABB

  !! insertion sort O(N^2). There exists cheaper sorting algorithm (quickSort O(N logN)) 
  !! but for N < 8 (which is mostly the case in FEM), insertion sort is better because quicksort
  !! needs extra procedures such as selecting pivots, ....
  !!
  !! sorts arraysReal so that its values are increasing with index. Array Int are 
  !! sorted using the exactly the same swaps made during arrayReal sorting process.
  ! (needs to be changed) (move to genericProcedures)
  subroutine sortPairs(self, arrayReal, arrayInt)
    class(cartesianGrid), intent(inout)              :: self
    real(defReal), dimension(:), intent(inout)       :: arrayReal
    integer(shortInt), dimension(:), intent(inout)   :: arrayInt
    integer(shortInt)                                :: i, j
    real(defReal)                                    :: key_real
    integer(shortInt)                                :: key_Int

    do i = 2, size(arrayReal)
       key_real = arrayReal(i);  key_Int = arrayInt(i)
       j = i - 1
       do while (j >= 1 .and. arrayReal(j) > key_real)
          arrayReal(j+1) = arrayReal(j)
          arrayInt(j+1) = arrayInt(j)
          j      = j - 1
       end do
       arrayReal(j+1) = key_real
       arrayInt(j+1) = key_Int
    end do

  end subroutine sortPairs

  !!
  !!
  !! 
  subroutine setGridIsOutside(self)
    class(cartesianGrid), intent(inout)              :: self
    integer(shortInt)                                :: i, j, k

    ! loop over all cartesian cells and call relevant subroutine
    do i = 1, self % n_xyz(1)
      do j = 1, self % n_xyz(2)
        do k = 1, self % n_xyz(3)
          call self % grid(i,j,k) % setIsOutside()
        end do
      end do 
    end do

  end subroutine setGridIsOutside

  !!
  !!
  !! (needs to be changed) (the binarySearch subroutine in generic procedure does not allow)
  !! ("value" to be outside the bounds of "array". So rewritten. Can be moved to genericProcedure Later)
  !! (needs to be changed) (for boundary edges, value > array(size(array)) or value < array(1) cases might not)
  !! (work using this subroutine)
  pure function binarySearchAngle(self, array, value) result(idx)
    class(cartesianGrid), intent(in)              :: self
    real(defReal), dimension(:), intent(in)       :: array
    real(defReal), intent(in)                     :: value
    integer(shortInt)                             :: idx, bottom, top, i

    ! in case of value being outside the ranges of "array", manually assign idx = size(array).
    if (value > array(size(array))) then
      idx = size(array)
      return
    elseif (value < array(1)) then
      idx = size(array)
      return
    end if

    ! Find Top and Bottom Index Array
    bottom = 1
    top = size(array)

    do i = 1,100
      !Calculate mid point
      idx = (top + bottom)*0.5

      ! Termination condition
      if (bottom == idx) return

      ! Binary Step
      if (array(idx) <= value) then
        bottom = idx
      else
        top = idx
      end if
    end do

  end function binarySearchAngle

  !!
  !!
  !!
  pure function getGridBounds_min(self) result(gridBounds_min)
    class(cartesianGrid), intent(in)              :: self
    real(defReal), dimension(3)                   :: gridBounds_min

    gridBounds_min = self % gridBounds_min

  end function getGridBounds_min

  !!
  !!
  !!
  elemental function getSpacingReciprocal(self) result(spacingReciprocal)
    class(cartesianGrid), intent(in)              :: self
    real(defReal)                                 :: spacingReciprocal

    spacingReciprocal = self % spacingReciprocal

  end function getSpacingReciprocal

  !!
  !!
  !!
  elemental function getGridWStar(self) result(wStar)
    class(cartesianGrid), intent(in)              :: self
    real(defReal)                                 :: wStar

    wStar = self % wStar

  end function getGridWStar

  !!
  !!
  !! 
  pure function getGridPhiCapital(self, cellIdxs) result(phiCapital)
    class(cartesianGrid), intent(in)              :: self
    integer(shortInt), dimension(3), intent(in)   :: cellIdxs
    integer(shortInt)                             :: phiCapital

    phiCapital = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getPhiCapital()

  end function getGridPhiCapital

  !!
  !!
  !! 
  pure function getGridPhi(self, cellIdxs) result(phi)
    class(cartesianGrid), intent(in)              :: self
    integer(shortInt), dimension(3), intent(in)   :: cellIdxs
    integer(shortInt)                             :: phi

    phi = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getPhi()

  end function getGridPhi

  !!
  !!
  !! 
  pure function getGridChi(self, cellIdxs) result(chi)
    class(cartesianGrid), intent(in)              :: self
    integer(shortInt), dimension(3), intent(in)   :: cellIdxs
    integer(shortInt)                             :: chi

    chi = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getChi()

  end function getGridChi

  !!
  !!
  !! (needs to be changed) (possible acceleration?)
  pure function getIsOutside(self, r) result(isOutside)
    class(cartesianGrid), intent(in)              :: self
    real(defReal), dimension(3), intent(in)       :: r
    logical                                       :: isOutside
    integer(shortInt)                             :: i

    isOutside = .FALSE.

    do i = 1, 3
      if (r(i) > self % gridBounds_max(i) .OR. r(i) < self % gridBounds_min(i)) then
        isOutside = .TRUE.
        return
      end if
    end do

  end function getIsOutside

end module CartesianGrid_class