module cartesianCellCoarsest_class
  
  use numPrecision
  use universalVariables,               only : ZERO
  use vertexShelf_class,                only : vertexShelf
  use edgeShelf_class,                  only : edgeShelf
  use faceShelf_class,                  only : faceShelf
  use elementShelf_class,               only : elementShelf
  use cartesianInitProcedures
  use cartesianGridSubLayer_inter,      only : cartesianGridSubLayer
  use cartesianGridFinest_class,        only : cartesianGridFinest
  use cartesianGridIntermediate1_class, only : cartesianGridIntermediate1
  use cartesianGenericProcedures,       only : quickSort, sortByHighestFrequency
  use genericProcedures,                only : append, fatalError
  use ragged3dMatrix_class,             only : ragged3d
  use dynamic2dMatSet_class,            only : dynamic2dMatSet

  implicit none
  private
  
  !!
  !!
  type, public                                          :: cartesianCellCoarsest
    private
    class(cartesianGridSubLayer), pointer               :: subGrid => null()
    integer(shortInt)                                   :: chi = 0 !-1 !"""
    integer(shortInt), dimension(:), allocatable        :: intersectedFaceIdxs !""(remove first)
    ! (needs to be changed) (cell centre should be a property to avoid repetative calc.)
    ! (due to limited memory, this is calculated each time needed)
    ! (try to avoid adding properties tho due to memory)

  contains

    ! Build procedures.
    procedure                                    :: cellTestPolyhedronInclusion
    procedure                                    :: cellFinitePrecision
    procedure                                    :: cellTestFaceIntersection
    procedure                                    :: setIsOutsideMesh
    procedure                                    :: refineCell
    procedure                                    :: refineCell2
    procedure                                    :: setChiFace !"""
    ! Runtime procedures.
    procedure                                    :: getChi
    procedure                                    :: getPhi
    procedure                                    :: getPhiCapital
    !procedure                                    :: getElementIdxs
    !procedure                                    :: getCandElemIdxs
    ! Analysis procedures.
    procedure                                    :: cellgetNumberOfCells

  end type cartesianCellCoarsest

contains

  !!
  !!
  !!
  subroutine cellTestPolyhedronInclusion(self, faces, currElementFaceIdxs, centroid, &
                                              faceNormalSigns, elementIdx)
    class(cartesianCellCoarsest), intent(inout)         :: self
    class(faceShelf), intent(in)                        :: faces
    integer(shortInt), dimension(:), intent(in)         :: currElementFaceIdxs
    real(defReal), dimension(3), intent(in)             :: centroid
    real(defReal), dimension(:,:), intent(in)           :: faceNormalSigns
    integer(shortInt), intent(in)                       :: elementIdx

    ! if current cell is found to be entirely contained within a polyhedron, there cannot not be other polyhedra
    if (self % chi > 0) return

    ! Otherwise, continue testing and appending the array of candidate element indices 
    call testPolyhedronInclusion2Coarsest(faces, currElementFaceIdxs, centroid, faceNormalSigns, elementIdx, self % chi)
    !""call append(self % candidateElementIdxs, elementIdx)

  end subroutine cellTestPolyhedronInclusion

  !!
  !!
  !!
  subroutine cellFinitePrecision(self, faces, elements, currElementFaceIdxs, centroid, &
                                 elementIdx)
    class(cartesianCellCoarsest), intent(inout)         :: self
    class(faceShelf), intent(in)                        :: faces
    class(elementShelf), intent(in)                     :: elements
    integer(shortInt), dimension(:), intent(in)         :: currElementFaceIdxs
    real(defReal), dimension(3), intent(in)             :: centroid
    integer(shortInt), intent(in)                       :: elementIdx
    integer(shortInt), dimension(:), allocatable        :: elementFaceIdxs, elementIdxs
    integer(shortInt)                                   :: i, faceIdx

    if (self % chi > 0) then
      return
    elseif (self % chi == 0) then
      if (.NOT. allocated(self % intersectedFaceIdxs)) then

      do i = 1, size(currElementFaceIdxs)

        if (dot_product(faces % getFaceNormal(currElementFaceIdxs(i)), centroid) &
            + faces % getFaceConst(currElementFaceIdxs(i)) > 0) then
            ! if TRUE, then centroid of this cell lies outside of the polyhedron
            return
        end if 

      end do

      ! If survived to this point, this cell is contained within a single mesh element, but
      ! chi mapping info was not assigned due to floating-point-finite-precision.
      self % chi = elementIdx




      end if
    end if
  end subroutine cellFinitePrecision

  !!
  !!
  !!
  subroutine cellTestFaceIntersection(self, vertices, edges, faces, currVertexIdxs, extraDistance, currFaceNormal, &
                                      centroid, cellSpacing, faceIdx, currFaceEdgeIdxs)
    class(cartesianCellCoarsest), intent(inout)         :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(in)                        :: edges
    class(faceShelf), intent(in)                        :: faces
    integer(shortInt), dimension(:), intent(in)         :: currVertexIdxs, currFaceEdgeIdxs
    real(defReal), dimension(3), intent(in)             :: currFaceNormal, centroid
    real(defReal), intent(in)                           :: extraDistance, cellSpacing
    integer(shortInt), intent(in)                       :: faceIdx

    if (self % chi > 0) then
      return
    elseif (self % chi < 0) then
      self % chi = 0
    end if
    
    call testFaceIntersectionCoarsest(vertices, edges, faces, currVertexIdxs, extraDistance, currFaceNormal, &
                                  centroid, cellSpacing, faceIdx, currFaceEdgeIdxs, self % intersectedFaceIdxs)

  end subroutine cellTestFaceIntersection

  ! !!
  ! !!
  ! !!
  ! subroutine setIsOutsideMesh(self)
  !   class(cartesianCellCoarsest), intent(inout)         :: self

  !   ! if candidateElementIdxs is not allocated, this cell does not intersect AABB of any polyhedron.
  !   ! Hence, this cell lies outside of mesh
  !   if (.NOT. allocated(self % candidateElementIdxs)) self % chi = -1

  ! end subroutine setIsOutsideMesh

  !!
  !!
  !!
  subroutine setIsOutsideMesh(self, no) !"""
    class(cartesianCellCoarsest), intent(inout)         :: self
    integer(shortInt), intent(in)                       :: no

    ! if candidateElementIdxs is not allocated, this cell does not intersect AABB of any polyhedron.
    ! Hence, this cell lies outside of mesh
    !@@
    !if (.NOT. allocated(self % candidateElementIdxs)) then
      if (.NOT. allocated(self % intersectedFaceIdxs)) then
        !self % chi = no
      end if
    !end if
    !@@

  end subroutine setIsOutsideMesh

  !!
  !!
  !!
  subroutine refineCell(self, vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                             newGridBoundsMin, alpha, wStar, circumscribedBallRadius, targetDistance, &
                             targetDistanceSqr)
    class(cartesianCellCoarsest), intent(inout)         :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(inout)                     :: edges
    class(faceShelf), intent(inout)                     :: faces
    class(elementShelf), intent(in)                     :: elements
    real(defReal), dimension(:), intent(in)             :: spacing, spacingInv
    integer(shortInt), dimension(:,:), intent(in)       :: n_xyz
    integer(shortInt), intent(in)                       :: n_layers
    real(defReal), dimension(3), intent(in)             :: newGridBoundsMin
    real(defReal), intent(in)                           :: alpha, wStar, circumscribedBallRadius, targetDistance, &
                                                           targetDistanceSqr
    integer(shortInt), dimension(:), allocatable        :: candidateElementIdxs, currElementFaceIdxs
    real(defReal), dimension(:), allocatable            :: extraDistanceArr
    integer(shortInt)                                   :: i, j, numberOfElements, numberOfFaces
    type(ragged3d)                                      :: normalSignsMat
    real(defReal), dimension(:,:), allocatable          :: currNormalSignsMat
    real(defReal)                                       :: halfSpacing

    !if the current cell is not entirely contained within a polyhedron nor outside of the mesh domain, refine the cell
    if (self % chi == 0) then

      ! If (very rarely) the cell is inside the mesh domain but intersects with no faces, then specify that this cell lies outside
      if (.NOT. allocated(self % intersectedFaceIdxs)) then
        self % chi = faces % getSize() + 1 !-1 !"""
        return
      end if

      halfSpacing = 0.5*spacing(2)

      ! First, setup element- and face-specific parameters so that they do not have to be recalculated for multiple cells at each layer.
      ! Because, we perform DFS for initialisation, this is benefitial (Otherwise, each cell in the same level need to
      ! recalculate these.) There is no memory issue because each travelsal down a path require (reduced list of) these
      ! parameters at each layer only once. It is important to calculate these during "refine" one by one because otherwise,
      ! all cells in the coarsest layer store these values.

      !@@ (remove)
      ! Begin with face-specific paremeters
      if (allocated(extraDistanceArr)) deallocate(extraDistanceArr)
      allocate(extraDistanceArr(size(self % intersectedFaceIdxs)))
      do i = 1, size(self % intersectedFaceIdxs)
        extraDistanceArr(i) = faces % getFaceExtraDistance(self % intersectedFaceIdxs(i))
      end do
      extraDistanceArr = extraDistanceArr * spacing(2)

      ! Then, element-specific parameters
      if (size(self % intersectedFaceIdxs) > 1) then
        if (allocated(candidateElementIdxs)) deallocate(candidateElementIdxs)
      !@@

        ! Preparation before setting up normalSignsMat
        candidateElementIdxs = faces % getFaceElementIdxs(self % intersectedFaceIdxs)
        
        numberOfElements = size(candidateElementIdxs)
        call quickSort(candidateElementIdxs, 1, numberOfElements)
        call normalSignsMat % kill()
        call normalSignsMat % init()

        ! Construct currNormalSignsMat for each face of the element and append to normalSignsMat
        do i = 1, numberOfElements
          currElementFaceIdxs = elements % getElementFaceIdxs(candidateElementIdxs(i))
          numberOfFaces = size(currElementFaceIdxs)

          if (allocated(currNormalSignsMat)) deallocate(currNormalSignsMat)
          allocate(currNormalSignsMat(3,numberOfFaces))
          do j = 1, numberOfFaces
            currNormalSignsMat(:,j) = faces % getFaceNormalSigns(currElementFaceIdxs(j))
          end do

          currNormalSignsMat(:,:) = currNormalSignsMat(:,:)*halfSpacing
          call normalSignsMat % append(currNormalSignsMat)

        end do

      elseif (size(self % intersectedFaceIdxs) == 1) then
        ! If there is a singe intersected face, we perform special version of polyhedron inclusion tests.
        ! So no need to construct normalSignsMat (initialisation is still required)
        call normalSignsMat % kill()
        call normalSignsMat % init()
    
      end if

      ! Determine and set the next layer
      if (n_layers > 2) then 
        allocate(cartesianGridIntermediate1:: self % subGrid)
      else
        allocate(cartesianGridFinest:: self % subGrid)
      end if

      ! print*, "intersectedFaceIdxs", self % intersectedFaceIdxs
      ! print*, "extraDistanceArr", extraDistanceArr
      !call fatalError("done", "done")


      !""
      call self % subgrid % init(vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                                 2, candidateElementIdxs, newGridBoundsMin, alpha, wStar)
      ! call self % subgrid % initt(vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
      !                             2, self % intersectedFaceIdxs, newGridBoundsMin, alpha, wStar, &
      !                             extraDistanceArr, candidateElementIdxs, normalSignsMat, &
      !                             circumscribedBallRadius, targetDistance, targetDistanceSqr)
      !""

    end if 

    ! deallocate candidateElementIdxs to save memory. Otherwise, memory could explode
    !""
    !if (allocated(self % candidateElementIdxs)) deallocate(self % candidateElementIdxs)
    if (allocated(self % intersectedFaceIdxs)) deallocate(self % intersectedFaceIdxs)
    !""

  end subroutine refineCell

  !!
  !!
  !!
  subroutine refineCell2(self, vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                             newGridBoundsMin, alpha, wStar, circumscribedBallRadius, targetDistance, &
                             targetDistanceSqr)
    class(cartesianCellCoarsest), intent(inout)         :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(inout)                     :: edges
    class(faceShelf), intent(inout)                     :: faces
    class(elementShelf), intent(in)                     :: elements
    real(defReal), dimension(:), intent(in)             :: spacing, spacingInv
    integer(shortInt), dimension(:,:), intent(in)       :: n_xyz
    integer(shortInt), intent(in)                       :: n_layers
    real(defReal), dimension(3), intent(in)             :: newGridBoundsMin
    real(defReal), intent(in)                           :: alpha, wStar, circumscribedBallRadius, targetDistance, &
                                                           targetDistanceSqr
    integer(shortInt), dimension(:), allocatable        :: candidateElementIdxs, currElementFaceIdxs, removedFaceIdxsInArr
    real(defReal), dimension(:), allocatable            :: extraDistanceArr
    integer(shortInt)                                   :: i, j, numberOfElements, numberOfFaces
    real(defReal), dimension(:,:), allocatable          :: currNormalSignsMat, faceNormalSigns
    real(defReal)                                       :: halfSpacing
    type(dynamic2dMatSet)                               :: normalSignsMat
    real(defReal), dimension(3)                         :: centroid
    logical(defBool)                                    :: isOut

    !if the current cell is not entirely contained within a polyhedron nor outside of the mesh domain, refine the cell
    if (self % chi == 0) then

      ! If (very rarely) the cell is inside the mesh domain but intersects with no faces, then specify that this cell lies outside
      if (.NOT. allocated(self % intersectedFaceIdxs)) then
        !@@
        self % chi = -(faces % getSize() + 1) !-1 !"""
        return


        !@@
      end if

      halfSpacing = 0.5*spacing(1)

      ! First, setup element- and face-specific parameters so that they do not have to be recalculated for multiple cells at each layer.
      ! Because, we perform DFS for initialisation, this is benefitial (Otherwise, each cell in the same level need to
      ! recalculate these.) There is no memory issue because each travelsal down a path require (reduced list of) these
      ! parameters at each layer only once. It is important to calculate these during "refine" one by one because otherwise,
      ! all cells in the coarsest layer store these values.

      ! Then, element-specific parameters
      ! if (size(self % intersectedFaceIdxs) > 1) then
        if (allocated(candidateElementIdxs)) deallocate(candidateElementIdxs)
      

        ! Preparation before setting up normalSignsMat

        ! Sort the candidate element indices such that when we retrive duplicatable list of element indices
        ! from faceElementIdxs of intersectedFaces, the most frequent element comes first because it is mostly 
        ! likely that refined cells will be contained within this element
        do i = 1, size(self % intersectedFaceIdxs)
          call append(candidateElementIdxs, faces % getFaceElementIdxs(self % intersectedFaceIdxs(i)))
        end do
        candidateElementIdxs = sortByHighestFrequency(candidateElementIdxs)
        
        numberOfElements = size(candidateElementIdxs)
        call normalSignsMat % kill()
        call normalSignsMat % init(numberOfElements)

        ! Construct currNormalSignsMat for each face of the element and append to normalSignsMat
        do i = 1, numberOfElements
          currElementFaceIdxs = elements % getElementFaceIdxs(candidateElementIdxs(i))
          numberOfFaces = size(currElementFaceIdxs)

          if (allocated(currNormalSignsMat)) deallocate(currNormalSignsMat)
          allocate(currNormalSignsMat(3,numberOfFaces))
          do j = 1, numberOfFaces
            currNormalSignsMat(:,j) = faces % getFaceNormalSigns(currElementFaceIdxs(j))
          end do
          currNormalSignsMat(:,:) = currNormalSignsMat(:,:)*halfSpacing
          call normalSignsMat % append(currNormalSignsMat, currElementFaceIdxs)

        end do

      !@@ if ther is a single intersectedFaceIdxs, then chi is non-zero, the code does not reach up to here.
      ! elseif (size(self % intersectedFaceIdxs) == 1) then
      !   ! If there is a singe intersected face, we perform special version of polyhedron inclusion tests.
      !   ! So no need to construct normalSignsMat
      !   call normalSignsMat % kill()
    
      ! end if

      ! Determine and set the next layer
      if (n_layers > 2) then 
        allocate(cartesianGridIntermediate1:: self % subGrid)
      else
        allocate(cartesianGridFinest:: self % subGrid)
      end if

      ! Update normalSignsMat
      do i = 1, 3
        centroid(i) = newGridBoundsMin(i) + spacing(1)*0.5
      end do

      do i = 1, numberOfElements
        if (allocated(faceNormalSigns)) deallocate(faceNormalSigns)
        if (allocated(currElementFaceIdxs)) deallocate(currElementFaceIdxs)
        if (allocated(removedFaceIdxsInArr)) deallocate(removedFaceIdxsInArr)
        call normalSignsMat % get_copy(i, faceNormalSigns, currElementFaceIdxs)
        call testPolyhedronInclusion2New(faces, currElementFaceIdxs, centroid, faceNormalSigns, candidateElementIdxs(i), &
                                      self % chi, removedFaceIdxsInArr, isOut, 1)

        if (isOut) then
          call normalSignsMat % delete(i)
        else
          if (allocated(removedFaceIdxsInArr)) then
            call normalSignsMat % delete_columns(i, removedFaceIdxsInArr)
            !print*, "Yes"
            !print*, size(removedFaceIdxsInArr)
          else
            !print*, "Not"
          end if
        end if
      
      end do
      call normalSignsMat % scale(spacingInv(1)*spacing(2))

      ! if (normalSignsMat % is_singleton()) then
      !   print*, "here"
      ! end if 
      !@@@

      ! call self % subgrid % init(vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
      !                            2, candidateElementIdxs, newGridBoundsMin, alpha, wStar)

      ! call self % subgrid % initt(vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
      !                             2, self % intersectedFaceIdxs, newGridBoundsMin, alpha, wStar, &
      !                             extraDistanceArr, candidateElementIdxs, normalSignsMat, &
      !                             circumscribedBallRadius, targetDistance, targetDistanceSqr)

      call self % subgrid % init2(vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                                 2, candidateElementIdxs, newGridBoundsMin, alpha, wStar, normalSignsMat)

      ! call fatalError("as", "as")
      !@@@

    end if 

    ! deallocate candidateElementIdxs to save memory. Otherwise, memory could explode
    !""
    !if (allocated(self % candidateElementIdxs)) deallocate(self % candidateElementIdxs)
    if (allocated(self % intersectedFaceIdxs)) deallocate(self % intersectedFaceIdxs)
    call normalSignsMat % kill()
    !""

  end subroutine refineCell2

  !!
  !!
  !!
  subroutine setChiFace(self, faces) !"""
    class(cartesianCellCoarsest), intent(inout)         :: self
    class(faceShelf), intent(in)                        :: faces

    if (allocated(self % intersectedFaceIdxs)) then
      if (size(self % intersectedFaceIdxs) == 1) then

        if (self % chi > 0) return
        if (self % chi == -(faces % getSize() + 1)) return
        self % chi = -(self % intersectedFaceIdxs(1))
        ! print*, self % chi

      end if
    end if

  end subroutine setChiFace

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! bit-trick (not saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! !!
  ! !!
  ! !!
  ! function getChi(self, baseIntegerCoord, shift, mask) result(chi)
  !   class(cartesianCellCoarsest), intent(in)            :: self
  !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
  !   integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
  !   integer(shortInt)                                   :: chi

  !   if (self % chi /= 0) then
  !     chi = self % chi
  !   else
  !     chi = self % subGrid % getGridChi(baseIntegerCoord, shift, mask, 2)
  !   end if

  ! end function getChi

  ! !!
  ! !!
  ! !!
  ! function getPhi(self, baseIntegerCoord, shift, mask) result(phi)
  !   class(cartesianCellCoarsest), intent(in)            :: self
  !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
  !   integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
  !   integer(shortInt)                                   :: phi

  !   phi = self % subGrid % getGridphi(baseIntegerCoord, shift, mask, 2)

  ! end function getPhi

  ! !!
  ! !!
  ! !!
  ! function getPhiCapital(self, baseIntegerCoord, shift, mask) result(phiCapital)
  !   class(cartesianCellCoarsest), intent(in)            :: self
  !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
  !   integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
  !   integer(shortInt)                                   :: phiCapital

  !   phiCapital = self % subGrid % getGridphiCapital(baseIntegerCoord, shift, mask, 2)

  ! end function getPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! bit-trick (saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !!
  !!
  !!
  function getChi(self, baseIntegerCoord, shift, mask, cellIdxsMat) result(chi)
    class(cartesianCellCoarsest), intent(in)            :: self
    integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
    integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
    integer(shortInt), dimension(:,:), intent(inout)    :: cellIdxsMat
    integer(shortInt)                                   :: chi

    if (self % chi /= 0) then
      chi = self % chi
    else
      chi = self % subGrid % getGridChi(baseIntegerCoord, shift, mask, 2, cellIdxsMat)
    end if

  end function getChi

  !!
  !!
  !!
  function getPhi(self, cellIdxsMat) result(phi)
    class(cartesianCellCoarsest), intent(in)            :: self
    integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
    integer(shortInt)                                   :: phi

    phi = self % subGrid % getGridphi(cellIdxsMat, 2)

  end function getPhi

  !!
  !!
  !!
  function getPhiCapital(self, cellIdxsMat) result(phiCapital)
    class(cartesianCellCoarsest), intent(in)            :: self
    integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
    integer(shortInt)                                   :: phiCapital

    phiCapital = self % subGrid % getGridphiCapital(cellIdxsMat, 2)

  end function getPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Non bit-trick (saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! !!
  ! !!
  ! !!
  ! function getChi(self, r, gridBounds_min, spacingInv, cellIdxsMat, nSub_xyz) result(chi)
  !   class(cartesianCellCoarsest), intent(in)            :: self
  !   real(defReal), dimension(3), intent(in)             :: r, gridBounds_min                                    
  !   real(defReal), dimension(:), intent(in)             :: spacingInv
  !   integer(shortInt), dimension(:,:), intent(inout)    :: cellIdxsMat
  !   integer(shortInt), dimension(:,:), intent(in)       :: nSub_xyz
  !   integer(shortInt)                                   :: chi

  !   if (self % chi /= 0) then
  !     chi = self % chi
  !   else
  !     chi = self % subGrid % getGridChi(r, gridBounds_min, spacingInv, 2, cellIdxsMat, nSub_xyz)
  !   end if

  ! end function getChi

  ! !!
  ! !!
  ! !!
  ! function getPhi(self, cellIdxsMat) result(phi)
  !   class(cartesianCellCoarsest), intent(in)            :: self
  !   integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
  !   integer(shortInt)                                   :: phi

  !   phi = self % subGrid % getGridphi(cellIdxsMat, 2)

  ! end function getPhi

  ! !!
  ! !!
  ! !!
  ! function getPhiCapital(self, cellIdxsMat) result(phiCapital)
  !   class(cartesianCellCoarsest), intent(in)            :: self
  !   integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
  !   integer(shortInt)                                   :: phiCapital

  !   phiCapital = self % subGrid % getGridphiCapital(cellIdxsMat, 2)

  ! end function getPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!   !!
!   !!
!   !!
!   function getElementIdxs(self) result(Idxs)
!     class(cartesianCellCoarsest), intent(in)                       :: self
!     integer(shortInt), dimension(:), allocatable                   :: Idxs

!     if (allocated(self % candidateElementIdxs)) then
!       Idxs = self % candidateElementIdxs
!     else
!       if (self % chi /= -1) then 
!         call fatalError("here", "here")
!       end if
!         allocate(Idxs(3))
!       Idxs(:) = 999
!     end if 
! end function getElementIdxs

! !!
! !!
! !!
! function getCandElemIdxs(self) result(arr)
!   class(cartesianCellCoarsest), intent(in)            :: self
!   integer(shortInt), dimension(:), allocatable        :: arr

!   arr = self%CandidateElementIdxs

! end function getCandElemIdxs

!!
!!
!!
function cellGetNumberOfCells(self, localNxyz, n_layers) result(output)
  class(cartesianCellCoarsest), intent(in)            :: self
  integer(shortInt), dimension(:,:), intent(in)       :: localNxyz
  integer(shortInt), intent(in)                       :: n_layers
  integer(shortInt), dimension(:), allocatable        :: output

  ! Allocate and initialise output array
  allocate(output(n_layers+1))
  output = 0

  if (self % chi /= 0) then
    output(1) = 1
  else
    output = self % subgrid % getNumberOfCells(localNxyz, n_layers, 2)
  end if
  
end function cellGetNumberOfCells



end module cartesianCellCoarsest_class