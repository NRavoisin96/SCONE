module cartesianGridCoarsest_class

  use universalVariables,              only : ZERO, INF, NUDGE
  use vertexShelf_class,               only : vertexShelf
  use edgeShelf_class,                 only : edgeShelf
  use faceShelf_class,                 only : faceShelf
  use elementShelf_class,              only : elementShelf
  use numPrecision      
  use cartesianCellCoarsest_class,     only : cartesianCellCoarsest
  use cartesianGenericProcedures
  use genericProcedures,               only : fatalError

  implicit none
  private

  !!
  !!
  !! stores information for all layers because non-coarsest layers have multiple of them per layer due to the adaptive-refinement 
  type, public                                                  :: cartesianGridCoarsest
    private
    real(defReal)                                               :: alpha = ZERO, l_min = ZERO, wStar = ZERO
    real(defReal), dimension(:), allocatable                    :: spacing, spacingInv
    real(defReal), dimension(3)                                 :: gridBounds_max = ZERO, gridBounds_min = ZERO, &
                                                                   meshBounds_max = ZERO, meshBounds_min = ZERO
    integer(shortInt), dimension(:,:), allocatable              :: n_xyz, shift, mask, nSub_xyz !!!!!(last one)
    type(cartesianCellCoarsest), dimension(:,:,:), allocatable  :: grid
    integer(shortInt)                                           :: n_layers
    ! (needs to be changed) (one-dimensinoalise matrixes if possible)

  contains

    ! Build procedures
    procedure                               :: init
    procedure                               :: constructMapping
    procedure                               :: setGridIsOutsideMesh
    procedure                               :: refineGrid
    ! Runtime procedures
    procedure                               :: getGridBounds_min
    procedure                               :: getSpacingReciprocal
    procedure                               :: getGridWStar
    procedure                               :: getGridIsOutsideBounds
    procedure                               :: getGridChi
    procedure                               :: getGridPhi
    procedure                               :: getGridPhiCapital
    ! Analysis procedures
    procedure                               :: getNumberOfCells

  end type cartesianGridCoarsest

contains

  !!
  !!
  !!
  subroutine init(self, vertices, edges, faces, elements)
    class(cartesianGridCoarsest), intent(inout)         :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(inout)                     :: edges
    class(faceShelf), intent(inout)                     :: faces
    class(elementShelf), intent(in)                     :: elements
    real(defReal), dimension(6)                         :: extremalCoordinates
    integer(shortInt)                                   :: i, j, diffExponent, exponentGapLayers, k
    integer(shortInt), dimension(:), allocatable        :: currEdgeVertexIdxs, numberOfCellsAnalysis
    real(defReal), dimension(3)                         :: currEdgeVector, extraRoom, n_xyzTarget, &
                                                           spacingComparison
    real(defReal)                                       :: maxCosValue, tempMaxCosValue, currEdgeLength, &
                                                           residual
    integer(shortInt), dimension(:,:), allocatable      :: minExponent

    !!!
    !integer(shortInt) :: baseIntegerCoord, globalIdx
    !!!

    ! (needs to be changed) (change it so that it can be read from the inputfile?)
    ! (Currently, n_layers can be only either 2 or 3; to be updated)
    self % n_layers = 4
    allocate(minExponent(self % n_layers,3))
    allocate(self % spacing(self % n_layers))
    allocate(self % spacingInv(self % n_layers))
    allocate(self % n_xyz(self % n_layers,3))
    allocate(self % shift(self % n_layers,3))
    allocate(self % mask(self % n_layers,3))
    allocate(self % nSub_xyz(self % n_layers,3))   !!!!!

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

    ! (needs to be changed) (written for temp operation; can be optimised further; calculate it during l_min calc.?)
    ! (can probably be a separte subroutine on its own, and the return value can be assigned as a grid attribute?)
    !print*, "Average edge length", calculateAvgEdgeLength(edges) 

    !-----------------------------------------------------------------------------------------
    ! set alpha
    !-----------------------------------------------------------------------------------------
    maxCosValue = findMinFaceAngle(edges, faces)

    tempMaxCosValue = findMinDihedralAngle(edges, faces, elements)
    if (maxCosValue < tempMaxCosValue) maxCosValue = tempMaxCosValue

    self % alpha = ACOS(maxCosValue)

    !-----------------------------------------------------------------------------------------
    ! calculate grid dimensions
    !-----------------------------------------------------------------------------------------
    !---------------------------------------------
    ! Repeatative process for the coarsest layer
    !---------------------------------------------
    ! (needs to be changed) (can combine the procedure for coarsest and finest layers)
    ! set the target cartesian grid spacing for the coarsest layer
    ! important: we set and keep the spacing in all directions the same
    self % spacing(1) = (self % l_min)*0.5*0.2
    self % spacingInv(1) = 1/(self % spacing(1))

    ! calculate true cartesian grid spacing and number of cells for the coarsest layer
    extremalCoordinates = vertices % getExtremalCoordinates()
    self % meshBounds_min = extremalCoordinates(1:3)
    self % meshBounds_max = extremalCoordinates(4:6)
    outer1: do i = 1, 3 
      ! extend the cartesian grid dimension slightly further from the bounds of the mesh
      self % gridBounds_min(i) = self % meshBounds_min(i) - NUDGE
      self % gridBounds_max(i) = self % meshBounds_max(i) + NUDGE

      ! calculate the target number of cells for the coarsest layer
      n_xyzTarget(i) = (self % gridBounds_max(i) - self % gridBounds_min(i)) / self % spacing(1)
      
      ! to use cheaper indexing calc. procedure, the number of cells must be a power of 2
      ! Hence, find the lowest possible j such that n_xyzTarget(i) <= 2**j   
      inner1: do j = 1, 30
        if (n_xyzTarget(i) <= 2**j) then
          minExponent(1,i) = j
          exit inner1
        end if
      end do inner1 

      ! set the true spacing and the number of cells for the coarsest layer
      self % n_xyz(1,i) = 2**minExponent(1,i)

      ! extend gridBounds so that spacing and n_xyz are kept the same.
      ! residual is added to the max bound only to increase the change of this extra space lying
      ! outside of the mesh bounds so that they do not have to be refined in sub-layers.
      !residual = self % n_xyz(1,i)*self % spacing(1) - (self % gridBounds_max(i) - self % gridBounds_min(i))
      residual = (self % n_xyz(1,i) - n_xyzTarget(i))*self % spacing(1)
      self % gridBounds_max(i) = self % gridBounds_max(i) + residual 

      ! print*, "spacing                         : ", self % spacing(1)
      ! print*, "No.                             : ", self % n_xyz(1,i)
      ! print*, "No. target                      : ", n_xyzTarget(i)
      ! print*, "Grid size after                 : ", (self%gridBounds_max(i)-self%gridBounds_min(i))
      ! print*, "Grid size before                : ", (self%gridBounds_max(i)-self%gridBounds_min(i)-residual)
      ! print*, "No. times spacing               : ", self%n_xyz(1,i)*self%spacing(1)
      ! print*, "residual                        : ", residual
      ! print*, "Grid lower bounds in xyz        : ", self % gridBounds_min
      ! print*, "Grid upper bounds in xyz        : ", self % gridBounds_max      
      ! Through fatal error in case grid bounds not set appropriately
      if (self%n_xyz(1,i)*self%spacing(1) /= (self%gridBounds_max(i)-self%gridBounds_min(i))) then
            call fatalError("Calculation of grid dimensions for the coarsest", "grid bounds not set appropriately")
      end if

    end do outer1

    !---------------------------------------------
    ! Repeatative process for the finest layer
    !---------------------------------------------
    ! calculate target spacing for the finest layer
    ! (needs to be changed) (times by 0.9999 for wStar?)
    self % wStar = (self % l_min)*min(0.5d0, SIN(self % alpha))
    self % spacing(self % n_layers) = 2*(self % wStar)*SIN(self % alpha)*SIN((self % alpha)/2)&
                                      /sqrt(3.0d0)/(1+SIN(self % alpha))/(1+SIN((self % alpha)/2))
    self % spacingInv(self % n_layers) = 1/(self % spacing(self % n_layers))

    ! calculate the target number of cells for the finest layer in x-direction
    ! Note that diffExponent will be the same for all directions (xyz) because spaings are the same in all directions
    ! for all layers, and number of cells are always a power of 2. Hence, the calculation is done for x-diretion only, 
    ! and the same result is applied to all other directions.
    n_xyzTarget(1) = (self % gridBounds_max(1) - self % gridBounds_min(1)) / self % spacing(self % n_layers)

    ! to use cheaper indexing calc. procedure, the number of cells must be a power of 2
    ! Hence, find the lowest possible j such that n_xyz(i) <= 2**j   
    loop1: do j = 1, 30
      if (n_xyzTarget(1) <= 2**j) then

        ! ensure the gap between the exponent of the coarsest and finest layers is large enough to accomodate the specified n_layers
        if (j - minExponent(1,1) < self % n_layers - 1) then
          diffExponent = self % n_layers - 1
        else
          diffExponent = j - minExponent(1,1)
        end if

        exit loop1

      end if
    end do loop1

    ! apply the result from x-direction to all xyz
    do i = 1, 3

      ! set minExponent
      minExponent(self % n_layers,i) = minExponent(1,i) + diffExponent

      ! set the number of cells for the finest layer
      self % n_xyz(self % n_layers,i) = 2**minExponent(self % n_layers,i)
      
    end do

    ! set the true spacing for the finest layer
    self % spacing(self % n_layers) = (self % gridBounds_max(1) - self % gridBounds_min(1)) / self % n_xyz(self % n_layers,1)
    self % spacingInv(self % n_layers) = 1/(self % spacing(self % n_layers))

    ! though fatal error for inappropirate cell number
    do i = 1, 3
      spacingComparison(i) = (self % gridBounds_max(i) - self % gridBounds_min(i)) / self % n_xyz(self % n_layers,i)
    end do
    if (spacingComparison(1) /= spacingComparison(2)) then
      call fatalError("Calculation of grid dimensions for the finest", "cell number not set appropriately")
    end if
    if (spacingComparison(1) /= spacingComparison(3)) then
      call fatalError("Calculation of grid dimensions for the finest", "cell number not set appropriately")
    end if

    !---------------------------------------------
    ! Repeat the same process for intermediate layers
    !---------------------------------------------
    if (self % n_layers > 2) then 
      ! find (approximately) equally spaced exponents for intermediate layers
      ! (needs to be changed) (parametric study needed to fine the best spacing strategy)
      exponentGapLayers = INT(diffExponent/(self % n_layers - 1))

      ! set minExponent, spacing and spacingInv for each of intermediate layers
      do i = 2, self % n_layers-1

        ! set minExponent
        do j = 1, 3
          minExponent(i,j) = minExponent(1,j) + exponentGapLayers*(i-1)
          self % n_xyz(i,j) = 2**minExponent(i,j)
        end do
        
        ! set grid spacing
        self % spacing(i) = (self % gridBounds_max(1) - self % gridBounds_min(1)) / (self % n_xyz(i,1))
        self % spacingInv(i) = 1/(self % spacing(i))  

        ! through fatal error in case of inappropriate spacings
        if (self % spacing(i) /= (self%gridBounds_max(2)-self%gridBounds_min(2))/(self % n_xyz(i,2))) then
          call fatalError("Calculation of grid dimensions for intermediate", "cell number not set appropriately")
        end if
        if (self % spacing(i) /= (self%gridBounds_max(3)-self%gridBounds_min(3))/(self % n_xyz(i,3))) then
          call fatalError("Calculation of grid dimensions for intermediate", "cell number not set appropriately")
        end if

      end do

    end if

    !-----------------------------------------------------------------------------------------
    ! Extra pre-calculation to help efficient indexing calculation procedure for each layer
    !-----------------------------------------------------------------------------------------

    ! calculate constants to be used for the cheaper indexing calc. procedure
    ! First, initialise for the coarsest layer
    do i = 1, 3
      self % shift(1,i) = minExponent(self % n_layers,i) - minExponent(1,i)
      self % mask(1,i) = 0 ! can be any integer (not used)
    end do

    ! then for all other layers
    do i = 2, self % n_layers
      do j = 1, 3
        self % shift(i,j) = minExponent(self % n_layers,j) - minExponent(i,j)
        self % mask(i,j) = ((self % n_xyz(i,j))/(self % n_xyz(i-1,j))) - 1
      end do
    end do

    !!!!!
    !-----------------------------------------------------------------------------------------
    ! Calculate the number of sub-division at each layer descending down the grid.
    !-----------------------------------------------------------------------------------------
    ! Each row refers to each layer, and column to dimension. Hence, the second row represents
    ! the number of sub-division when going from the coarsest to the first intermediate layer. 
    ! First row is the dummy one to avoid having to subtract current layer by 1 during in-cycle.
    self % nSub_xyz(:,:) = self % mask(:,:) + 1
    !!!!!

    !-----------------------------------------------------------------------------------------
    ! print cartesian grid parameters and mesh quality
    !-----------------------------------------------------------------------------------------
    print*, "----------------------------------------------------"
    print*, "/\/\ Cartesian grid parameters and mesh quality /\/\"
    print*, "Minimum angle                       : ", self % alpha
    print*, "Minimum edge length                 : ", self % l_min
    print*, "No. of vertices                     : ", vertices % getSize()
    print*, "No. of edges                        : ", edges % getSize()
    print*, "No. of faces                        : ", faces % getSize()
    print*, "No. of elements                     : ", elements % getSize()
    print*, "No. of layers                       : ", self % n_layers
    print*, "Grid spacing for coarsest           : ", self % spacing(1)
    if (self % n_layers > 2) then
      do i = 2, self % n_layers - 1 
        print*, "Grid spacing for intermediate       : ", self % spacing(i)
      end do
    end if
    print*, "Grid spacing for finest             : ", self % spacing(self % n_layers)
    print*, "Grid size in xyz for coarsest       : ", self % n_xyz(1,:)
    if (self % n_layers > 2) then
      do i = 2, self % n_layers - 1 
        print*, "Grid size in xyz for intermediate   : ", self % n_xyz(i,:)
      end do
    end if
    print*, "Grid size in xyz for finest         : ", self % n_xyz(self % n_layers,:)
    print*, "Grid lower bounds in xyz            : ", self % gridBounds_min
    print*, "Grid upper bounds in xyz            : ", self % gridBounds_max
    print*, "----------------------------------------------------"

    ! for temporary debugging
    ! print*, "/////////////////////////////////////////"
    ! print*, "spacing of each layer"
    ! print*, self % spacing
    ! print*, "/////////////////////////////////////////"
    ! print*, "grid size of each layer in each direction"
    ! print*, self % n_xyz
    ! print*, "/////////////////////////////////////////"

    ! for temporary debugging
    ! print*, self % n_xyz(2,:)
    ! baseIntegerCoord = floor((0.71-self % gridBounds_min(1)) * self%spacingInv(self%n_layers))
    ! print*, baseIntegerCoord
    ! do i = 1, 3
    !   print*, "-----------------------------------------"
    !   print*, i
    !   print*, self%shift(i,:)
    !   print*, self%mask(i,:)
    !   globalIdx = ishft(baseIntegerCoord,-self%shift(i,1))
    !   print*, globalIdx
    !   print*, iand(globalIdx, self%mask(i,1)) + 1
    ! end do

    !-----------------------------------------------------------------------------------------
    !initialise for patch search
    !-----------------------------------------------------------------------------------------
    allocate(self % grid(self % n_xyz(1,1), self % n_xyz(1,2), self % n_xyz(1,3)))

    call self % constructMapping(vertices, edges, faces, elements)
    call self % setGridIsOutsideMesh()
    call self % refineGrid(vertices, edges, faces, elements)

    ! Print the number of cells for each layer (for extra analysis).
    ! print*, "Starting the procedure to calculate number of cells for each type"
    ! numberOfCellsAnalysis = self % getNumberOfCells()
    ! print*, numberOfCellsAnalysis
    ! call fatalError("Init, cartesianGridCoarsest_class.f90", "Terminating after printing &
    !                 the number of cells for each type. If not intended, comment these lines.")

    ! ! temporary for debugging
    ! do i = 1, self % n_xyz(1,1)
    !   do j = 1, self % n_xyz(1,2)
    !     do k = 1, self % n_xyz(1,3)
    !       print*, "--------------------------------"
    !       print*,  self % grid(i,j,k) % getChi()
    !       print*, self % grid(i,j,k) % getElementIdxs()
    !     end do
    !   end do
    ! end do

    !!!!!
    ! print*,"***************************************"
    ! print*, "Coarsest layer global index", getGlobalIdx([154,37,175], self%shift(1,:))
    ! print*, "Intermediate layer local index", getlocalIdx([154,37,175], self%shift(2,:), self%mask(2,:))
    ! print*, "Finest layer local index",getlocalIdxFinest([154,37,175], self%shift(3,:), self%mask(3,:))
    !print*, "gridChi", self%getGridChi([154,37,175])
    !!!!!

  end subroutine init

  !!
  !!
  !!
  subroutine constructMapping(self, vertices, edges, faces, elements)
    class(cartesianGridCoarsest), intent(inout)           :: self
    class(vertexShelf), intent(in)                        :: vertices
    class(edgeShelf), intent(inout)                       :: edges
    class(faceShelf), intent(in)                          :: faces
    class(elementShelf), intent(in)                       :: elements
    integer(shortInt)                                     :: i, j, k, l
    integer(shortInt), dimension(:), allocatable          :: currVertexIdxs, currElementFaceIdxs
    integer(shortInt), dimension(6)                       :: AABBIndices
    real(defReal), dimension(3)                           :: centroid, currFaceNormal
    real(defReal), dimension(:,:), allocatable            :: faceNormalSigns

    !----------------------------------------------------------------------------------------------
    ! polyhedron inclusion tests
    !----------------------------------------------------------------------------------------------
    allocate(faceNormalSigns(3, 2))
    do i = 1, elements % getSize()

        ! construct box for candidate cells
        currVertexIdxs = elements % getElementVertexIdxs(i)
        AABBIndices = constructAABB(vertices, currVertexIdxs, self % gridBounds_min, self % spacing(1))

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
        faceNormalSigns = faceNormalSigns * (self % spacing(1))/2


        !Loop over all cartesian cells in the box and test if each cell is entirely included in the polyhedron
        !(needs to be changed) (k and l can be a function of j e.g. k = datum + slope*j so that box is narrowed down)
        do j = AABBIndices(1), AABBIndices(4)
            do k = AABBIndices(2), AABBIndices(5)
                do l = AABBIndices(3), AABBIndices(6)

                    ! (needs to be changed) (store centroid info)
                    centroid(1) = (self % gridBounds_min(1)) + (self % spacing(1)) * (j-0.5)
                    centroid(2) = (self % gridBounds_min(2)) + (self % spacing(1)) * (k-0.5)
                    centroid(3) = (self % gridBounds_min(3)) + (self % spacing(1)) * (l-0.5)

                    call self % grid(j,k,l) % cellTestPolyhedronInclusion(faces, currElementFaceIdxs, centroid, &
                                                                          faceNormalSigns, i)

                end do 
            end do    
        end do

    end do

  end subroutine constructMapping

  !!
  !!
  !!
  subroutine setGridIsOutsideMesh(self)
    class(cartesianGridCoarsest), intent(inout)         :: self
    integer(shortInt)                                   :: i, j, k

    do i = 1, self % n_xyz(1,1)
      do j = 1, self % n_xyz(1,2)
        do k = 1, self % n_xyz(1,3)
          call self % grid(i,j,k) % setIsOutsideMesh()
        end do
      end do
    end do

  end subroutine setGridIsOutsideMesh

  !!
  !!
  !!
  subroutine refineGrid(self, vertices, edges, faces, elements)
    class(cartesianGridCoarsest), intent(inout)         :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(inout)                     :: edges
    class(faceShelf), intent(inout)                     :: faces
    class(elementShelf), intent(in)                     :: elements
    integer(shortInt)                                   :: i, j, k
    real(defReal), dimension(3)                         :: newGridBoundsMin

    do i = 1, self % n_xyz(1,1)
      do j = 1, self % n_xyz(1,2)
        do k = 1, self % n_xyz(1,3)

          ! (needs to be changed) (store newGridBoundsMin info)
          newGridBoundsMin(1) = (self % gridBounds_min(1)) + (self % spacing(1)) * (i-1)
          newGridBoundsMin(2) = (self % gridBounds_min(2)) + (self % spacing(1)) * (j-1)
          newGridBoundsMin(3) = (self % gridBounds_min(3)) + (self % spacing(1)) * (k-1)

          call self % grid(i,j,k) % refineCell(vertices, edges, faces, elements, self % spacing, &
                                               self % spacingInv, self % n_xyz, self % n_layers, &
                                               newGridBoundsMin, self % alpha, self % wStar)
        end do
      end do
    end do

  end subroutine refineGrid

  !!
  !!
  !!
  function getGridBounds_min(self) result(gridBounds_min)
    class(cartesianGridCoarsest), intent(in)            :: self
    real(defReal), dimension(3)                         :: gridBounds_min

    gridBounds_min = self % gridBounds_min

  end function getGridBounds_min

  !!
  !!
  !!
  function getSpacingReciprocal(self) result(spacingReciprocal)
    class(cartesianGridCoarsest), intent(in)            :: self
    real(defReal)                                       :: spacingReciprocal

    spacingReciprocal = self % spacingInv(self % n_layers)

  end function getSpacingReciprocal

  !!
  !!
  !!
  function getGridWStar(self) result(wStar)
    class(cartesianGridCoarsest), intent(in)            :: self
    real(defReal)                                       :: wStar

    wStar = self % wStar

  end function getGridWStar

  !!
  !!
  !! (needs to be changed) (possible acceleration?)
  !! (needs to be changed) (meshBounds here are AABB of the entire mesh. Hence,
  !! there can be a problem if mesh bounds are not perfectly cubic but somewhat irregular.
  !! more rigorous analysis required.)
  function getGridIsOutsideBounds(self, r) result(isOutside)
    class(cartesianGridCoarsest), intent(in)            :: self
    real(defReal), dimension(3), intent(in)             :: r
    logical                                             :: isOutside
    integer(shortInt)                                   :: i

    isOutside = .FALSE.

    do i = 1, 3
      if (r(i) > self % meshBounds_max(i) .OR. r(i) < self % meshBounds_min(i)) then
        isOutside = .TRUE.
        return
      end if
    end do

  end function getGridIsOutsideBounds

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! bit-trick (not saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! !!
  ! !!
  ! !! 
  ! function getGridChi(self, baseIntegerCoord) result(chi)
  !   class(cartesianGridCoarsest), intent(in)            :: self
  !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
  !   integer(shortInt)                                   :: chi
  !   integer(shortInt), dimension(3)                     :: cellIdxs

  !   cellIdxs = getGlobalIdx(baseIntegerCoord,self%shift(1,:))

  !   chi = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getChi(baseIntegerCoord, &
  !                                                                     self%shift,self%mask)

  ! end function getGridChi

  ! !!
  ! !!
  ! !! 
  ! function getGridPhi(self, baseIntegerCoord) result(phi)
  !   class(cartesianGridCoarsest), intent(in)            :: self
  !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
  !   integer(shortInt)                                   :: phi
  !   integer(shortInt), dimension(3)                     :: cellIdxs

  !   cellIdxs = getGlobalIdx(baseIntegerCoord,self%shift(1,:))

  !   phi = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getPhi(baseIntegerCoord, &
  !                                                                     self%shift,self%mask)

  ! end function getGridPhi

  ! !!
  ! !!
  ! !! 
  ! function getGridPhiCapital(self, baseIntegerCoord) result(phiCapital)
  !   class(cartesianGridCoarsest), intent(in)            :: self
  !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
  !   integer(shortInt)                                   :: phiCapital
  !   integer(shortInt), dimension(3)                     :: cellIdxs

  !   cellIdxs = getGlobalIdx(baseIntegerCoord,self%shift(1,:))

  !   phiCapital = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getPhiCapital(&
  !                                             baseIntegerCoord, self%shift,self%mask)

  ! end function getGridPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! bit-trick (saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !!
  !!
  !! 
  function getGridChi(self, baseIntegerCoord, cellIdxsMat) result(chi)
    class(cartesianGridCoarsest), intent(in)                       :: self
    integer(shortInt), dimension(3), intent(in)                    :: baseIntegerCoord
    integer(shortInt), dimension(:,:), allocatable, intent(inout)  :: cellIdxsMat
    integer(shortInt)                                              :: chi

    if (.NOT. allocated(cellIdxsMat)) allocate(cellIdxsMat(self%n_layers,3))
    cellIdxsMat(1,:) = getGlobalIdx(baseIntegerCoord,self%shift(1,:))

    chi = self % grid(cellIdxsMat(1,1), cellIdxsMat(1,2), cellIdxsMat(1,3)) % &
                    getChi(baseIntegerCoord, self%shift,self%mask, cellIdxsMat)

  end function getGridChi

  !!
  !!
  !! 
  function getGridPhi(self, cellIdxsMat) result(phi)
    class(cartesianGridCoarsest), intent(in)            :: self
    integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
    integer(shortInt)                                   :: phi

    phi = self % grid(cellIdxsMat(1,1), cellIdxsMat(1,2), cellIdxsMat(1,3)) % &
                                                            getPhi(cellIdxsMat)

  end function getGridPhi

  !!
  !!
  !! 
  function getGridPhiCapital(self, cellIdxsMat) result(phiCapital)
    class(cartesianGridCoarsest), intent(in)            :: self
    integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
    integer(shortInt)                                   :: phiCapital

    phiCapital = self % grid(cellIdxsMat(1,1), cellIdxsMat(1,2), cellIdxsMat(1,3)) &
                                                        % getPhiCapital(cellIdxsMat)

  end function getGridPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Non bit-trick (saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! !!
  ! !!
  ! !! 
  ! function getGridChi(self, r, cellIdxsMat) result(chi)
  !   class(cartesianGridCoarsest), intent(in)                       :: self
  !   real(defReal), dimension(3), intent(in)                        :: r                                 
  !   integer(shortInt), dimension(:,:), allocatable, intent(inout)  :: cellIdxsMat
  !   integer(shortInt)                                              :: chi

  !   if (.NOT. allocated(cellIdxsMat)) allocate(cellIdxsMat(self%n_layers,3))
  !   cellIdxsMat(1,:) = ceiling((r(:) - self%gridBounds_min(:))*(self%spacingInv(1)))

  !   chi = self % grid(cellIdxsMat(1,1), cellIdxsMat(1,2), cellIdxsMat(1,3)) % &
  !         getChi(r, self%gridBounds_min, self%spacingInv, cellIdxsMat, self%nSub_xyz)

  ! end function getGridChi

  ! !!
  ! !!
  ! !! 
  ! function getGridPhi(self, cellIdxsMat) result(phi)
  !   class(cartesianGridCoarsest), intent(in)            :: self
  !   integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
  !   integer(shortInt)                                   :: phi

  !   phi = self % grid(cellIdxsMat(1,1), cellIdxsMat(1,2), cellIdxsMat(1,3)) % &
  !                                                           getPhi(cellIdxsMat)

  ! end function getGridPhi

  ! !!
  ! !!
  ! !! 
  ! function getGridPhiCapital(self, cellIdxsMat) result(phiCapital)
  !   class(cartesianGridCoarsest), intent(in)            :: self
  !   integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
  !   integer(shortInt)                                   :: phiCapital

  !   phiCapital = self % grid(cellIdxsMat(1,1), cellIdxsMat(1,2), cellIdxsMat(1,3)) &
  !                                                       % getPhiCapital(cellIdxsMat)

  ! end function getGridPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  !!
  !!
  !!
  function getNumberOfCells(self) result(numberOfCells)
    class(cartesianGridCoarsest), intent(in)         :: self
    integer(shortInt), dimension(:), allocatable     :: numberOfCells
    integer(shortInt)                                :: i, j, k, l
    integer(shortInt), dimension(:), allocatable     :: output
    integer(shortInt), dimension(:,:), allocatable   :: localNxyz

    ! Allocate array and matrix
    allocate(numberOfCells(self % n_layers + 1))
    numberOfCells = 0
    allocate(localNxyz(self % n_layers - 1, 3))

    ! Initialise localNxyz. It starts from the second layer in the first row.
    ! Three different columns represent three different directions. 
    do i = 1, self % n_layers - 1 
      do j = 1, 3
        localNxyz(i,j) = self % n_xyz(i+1, j)/self % n_xyz(i, j)
      end do
    end do

    ! Loop through all Cartesian cells in the coarsest layer. If needed, it descend down the layers.
    do i = 1, self % n_xyz(1,1)
      do j = 1, self % n_xyz(1,2)
        do k = 1, self % n_xyz(1,3)

          ! Retrive an array containing the number of cells for each type.
          output = self % grid(i,j,k) % cellGetNumberOfCells(localNxyz, self % n_layers)

          ! Update the number of cells for each type.
          ! The first integer represents the number of cells in coarsest layer
          ! The second integer represents the number of cells in the second coarsest layer
          ! The second last integer represent the number of cells in the finest layer
          ! The last integer represents the number of cells in the finest layer without phiCaptial mapping.
          do l = 1, self % n_layers + 1
            numberOfCells(l) = numberOfCells(l) + output(l)
          end do

        end do
      end do
    end do
    
  end function getNumberOfCells
    
end module cartesianGridCoarsest_class
