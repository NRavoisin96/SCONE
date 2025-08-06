module cartesianGridIntermediate2_class

  use vertexShelf_class,                       only : vertexShelf
  use edgeShelf_class,                         only : edgeShelf
  use faceShelf_class,                         only : faceShelf
  use elementShelf_class,                      only : elementShelf
  use cartesianGridSubLayer_inter,             only : cartesianGridSubLayer
  use cartesianCellIntermediate2_class,        only : cartesianCellIntermediate2
  use numPrecision   
  use cartesianGenericProcedures

  implicit none
  private


  !!
  !!
  !!
  type, public, extends(cartesianGridSubLayer)                      :: cartesianGridIntermediate2
    private
    type(cartesianCellIntermediate2), dimension(:,:,:), allocatable  :: grid

  contains

    ! Build procedures
    procedure                                    :: init
    procedure                                    :: constructMapping
    procedure                                    :: refineGrid
    ! Runtime procedures
    procedure                                    :: getGridChi
    procedure                                    :: getGridPhi
    procedure                                    :: getGridPhiCapital
    ! Analysis procedures
    procedure                                    :: getNumberOfCells

  end type cartesianGridIntermediate2

contains

  !!
  !!
  !!
  subroutine init(self, vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                       currLayer, candidateElementIdxs, gridBoundsMin, alpha, wStar)
    class(cartesianGridIntermediate2), intent(inout)     :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(inout)                     :: edges
    class(faceShelf), intent(inout)                     :: faces
    class(elementShelf), intent(in)                     :: elements
    real(defReal), dimension(:), intent(in)             :: spacing, spacingInv
    integer(shortInt), dimension(:,:), intent(in)       :: n_xyz
    integer(shortInt), intent(in)                       :: n_layers, currLayer
    integer(shortInt), dimension(:), intent(in)         :: candidateElementIdxs
    real(defReal), dimension(3), intent(in)             :: gridBoundsMin
    real(defReal), intent(in)                           :: alpha, wStar
    integer(shortInt), dimension(3)                     :: localNxyz
    integer(shortInt)                                   :: i

    ! calculate local number of cells 
    do i = 1, 3
      localNxyz(i) = n_xyz(currLayer,i)/n_xyz(currLayer-1,i)
    end do

    ! construct mapping and refine further
    allocate(self % grid(localNxyz(1), localNxyz(2), localNxyz(3)))
    call self % constructMapping(vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                                 currLayer, candidateElementIdxs, gridBoundsMin, localNxyz)
    call self % refineGrid(vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                           currLayer, candidateElementIdxs, gridBoundsMin, localNxyz, alpha, wStar)

  end subroutine

  !!
  !!
  !!
  subroutine constructMapping(self, vertices, edges, faces, elements, spacing, spacingInv, n_xyz, &
                                   n_layers, currLayer, candidateElementIdxs, gridBoundsMin, localNxyz)
    class(cartesianGridIntermediate2), intent(inout)     :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(inout)                     :: edges
    class(faceShelf), intent(inout)                     :: faces
    class(elementShelf), intent(in)                     :: elements
    real(defReal), dimension(:), intent(in)             :: spacing, spacingInv
    real(defReal), dimension(3), intent(in)             :: gridBoundsMin
    integer(shortInt), dimension(:,:), intent(in)       :: n_xyz
    integer(shortInt), intent(in)                       :: n_layers, currLayer
    integer(shortInt), dimension(3), intent(in)         :: localNxyz
    integer(shortInt)                                   :: i, j, k, l
    integer(shortInt), dimension(:), allocatable        :: currElementFaceIdxs
    real(defReal), dimension(3)                         :: centroid, currFaceNormal
    real(defReal), dimension(:,:), allocatable          :: faceNormalSigns
    integer(shortInt), dimension(:), intent(in)         :: candidateElementIdxs

    !----------------------------------------------------------------------------------------------
    ! polyhedron inclusion tests
    !----------------------------------------------------------------------------------------------
    allocate(faceNormalSigns(3, 2))
    do i = 1, size(candidateElementIdxs)

        ! calculate element-only-dependent properties
        currElementFaceIdxs = elements % getElementFaceIdxs(candidateElementIdxs(i))

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
        faceNormalSigns = faceNormalSigns * (spacing(currLayer))/2


        !Loop over all cartesian cells in the box and test if each cell is entirely included in the polyhedron
        !(needs to be changed) (k and l can be a function of j e.g. k = datum + slope*j so that box is narrowed down)
        do j = 1, localNxyz(1)
            do k = 1, localNxyz(2)
                do l = 1, localNxyz(3)

                    ! (needs to be changed) (store centroid info)
                    centroid(1) = (gridBoundsMin(1)) + (spacing(currLayer)) * (j-0.5)
                    centroid(2) = (gridBoundsMin(2)) + (spacing(currLayer)) * (k-0.5)
                    centroid(3) = (gridBoundsMin(3)) + (spacing(currLayer)) * (l-0.5)

                    call self % grid(j,k,l) % cellTestPolyhedronInclusion(faces, currElementFaceIdxs, centroid, &
                                                                          faceNormalSigns, candidateElementIdxs(i))

                end do 
            end do    
        end do

    end do

  end subroutine constructMapping 

  !!
  !!
  !!
  subroutine refineGrid(self, vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                             currLayer, candidateElementIdxs, gridBoundsMin, localNxyz, alpha, wStar)
    class(cartesianGridIntermediate2), intent(inout)     :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(inout)                     :: edges
    class(faceShelf), intent(inout)                     :: faces
    class(elementShelf), intent(in)                     :: elements
    real(defReal), dimension(:), intent(in)             :: spacing, spacingInv
    integer(shortInt), dimension(:,:), intent(in)       :: n_xyz
    integer(shortInt), intent(in)                       :: n_layers, currLayer
    integer(shortInt), dimension(:), intent(in)         :: candidateElementIdxs
    real(defReal), dimension(3), intent(in)             :: gridBoundsMin
    integer(shortInt), dimension(3), intent(in)         :: localNxyz
    real(defReal), intent(in)                           :: alpha, wStar
    integer(shortInt)                                   :: i, j, k
    real(defReal), dimension(3)                         :: newGridBoundsMin

    do i = 1, localNxyz(1)
      do j = 1, localNxyz(2)
        do k = 1, localNxyz(3)

          ! (needs to be changed) (store newGridBoundsMin info)
          newGridBoundsMin(1) = (gridBoundsMin(1)) + (spacing(currLayer)) * (i-1)
          newGridBoundsMin(2) = (gridBoundsMin(2)) + (spacing(currLayer)) * (j-1)
          newGridBoundsMin(3) = (gridBoundsMin(3)) + (spacing(currLayer)) * (k-1)

          call self % grid(i,j,k) % refineCell(vertices, edges, faces, elements, spacing, spacingInv, n_xyz, &
                                               n_layers, currLayer, candidateElementIdxs, newGridBoundsMin, &
                                               alpha, wStar)
                                               
        end do
      end do
    end do

  end subroutine refineGrid

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! bit-trick (not saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! !!
  ! !!
  ! !! 
  ! function getGridChi(self, baseIntegerCoord, shift, mask, currLayer) result(chi)
  !   class(cartesianGridIntermediate2), intent(in)        :: self
  !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
  !   integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: chi
  !   integer(shortInt), dimension(3)                     :: cellIdxs

  !   cellIdxs = getLocalIdx(baseIntegerCoord, shift(currLayer,:), mask(currLayer,:))

  !   chi = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getChi(baseIntegerCoord, &
  !                                                                   shift,mask,currLayer)

  ! end function getGridChi

  ! !!
  ! !!
  ! !! 
  ! function getGridPhi(self, baseIntegerCoord, shift, mask, currLayer) result(phi)
  !   class(cartesianGridIntermediate2), intent(in)        :: self
  !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
  !   integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: phi
  !   integer(shortInt), dimension(3)                     :: cellIdxs

  !   cellIdxs = getLocalIdx(baseIntegerCoord, shift(currLayer,:), mask(currLayer,:))

  !   phi = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getPhi(baseIntegerCoord, &
  !                                                                   shift,mask,currLayer)

  ! end function getGridPhi

  ! !!
  ! !!
  ! !! 
  ! function getGridPhiCapital(self, baseIntegerCoord, shift, mask, currLayer) result(phiCapital)
  !   class(cartesianGridIntermediate2), intent(in)        :: self
  !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
  !   integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: phiCapital
  !   integer(shortInt), dimension(3)                     :: cellIdxs

  !   cellIdxs = getLocalIdx(baseIntegerCoord, shift(currLayer,:), mask(currLayer,:))

  !   phiCapital = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getPhiCapital(&
  !                                             baseIntegerCoord, shift,mask,currLayer)

  ! end function getGridPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! bit-trick (saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !!
  !!
  !! 
  function getGridChi(self, baseIntegerCoord, shift, mask, currLayer, cellIdxsMat) result(chi)
    class(cartesianGridIntermediate2), intent(in)        :: self
    integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
    integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
    integer(shortInt), dimension(:,:), intent(inout)    :: cellIdxsMat
    integer(shortInt), intent(in)                       :: currLayer
    integer(shortInt)                                   :: chi

    cellIdxsMat(currLayer,:) = getLocalIdx(baseIntegerCoord, shift(currLayer,:), mask(currLayer,:))

    chi = self % grid(cellIdxsMat(currLayer,1), cellIdxsMat(currLayer,2), cellIdxsMat(currLayer,3)) &
                                         % getChi(baseIntegerCoord,shift,mask,currLayer, cellIdxsMat)

  end function getGridChi

  !!
  !!
  !! 
  function getGridPhi(self, cellIdxsMat, currLayer) result(phi)
    class(cartesianGridIntermediate2), intent(in)        :: self
    integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
    integer(shortInt), intent(in)                       :: currLayer
    integer(shortInt)                                   :: phi

    phi = self % grid(cellIdxsMat(currLayer,1), cellIdxsMat(currLayer,2), cellIdxsMat(currLayer,3)) &
                                                                      % getPhi(cellIdxsMat,currLayer)

  end function getGridPhi

  !!
  !!
  !! 
  function getGridPhiCapital(self, cellIdxsMat, currLayer) result(phiCapital)
    class(cartesianGridIntermediate2), intent(in)        :: self
    integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
    integer(shortInt), intent(in)                       :: currLayer
    integer(shortInt)                                   :: phiCapital

    phiCapital = self % grid(cellIdxsMat(currLayer,1), cellIdxsMat(currLayer,2), cellIdxsMat(currLayer,3)) &
                                                                     % getPhiCapital(cellIdxsMat ,currLayer)

  end function getGridPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Non bit-trick (saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! !!
  ! !!
  ! !! 
  ! function getGridChi(self, r, gridBounds_min, spacingInv, currLayer, cellIdxsMat, nSub_xyz) result(chi)
  !   class(cartesianGridIntermediate2), intent(in)        :: self
  !   real(defReal), dimension(3), intent(in)             :: r, gridBounds_min                                    
  !   real(defReal), dimension(:), intent(in)             :: spacingInv
  !   integer(shortInt), dimension(:,:), intent(inout)    :: cellIdxsMat
  !   integer(shortInt), dimension(:,:), intent(in)       :: nSub_xyz
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: chi

  !   cellIdxsMat(currLayer,:) = floor((r(:) - gridBounds_min(:))*(spacingInv(currLayer)))
  !   cellIdxsMat(currLayer,:) = mod(cellIdxsMat(currLayer,:),nSub_xyz(currLayer,:)) + 1

  !   chi = self % grid(cellIdxsMat(currLayer,1), cellIdxsMat(currLayer,2), cellIdxsMat(currLayer,3)) &
  !                           % getChi(r, gridBounds_min, spacingInv, currLayer, cellIdxsMat, nSub_xyz)

  ! end function getGridChi

  ! !!
  ! !!
  ! !! 
  ! function getGridPhi(self, cellIdxsMat, currLayer) result(phi)
  !   class(cartesianGridIntermediate2), intent(in)        :: self
  !   integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: phi

  !   phi = self % grid(cellIdxsMat(currLayer,1), cellIdxsMat(currLayer,2), cellIdxsMat(currLayer,3)) &
  !                                                                     % getPhi(cellIdxsMat,currLayer)

  ! end function getGridPhi

  ! !!
  ! !!
  ! !! 
  ! function getGridPhiCapital(self, cellIdxsMat, currLayer) result(phiCapital)
  !   class(cartesianGridIntermediate2), intent(in)        :: self
  !   integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: phiCapital

  !   phiCapital = self % grid(cellIdxsMat(currLayer,1), cellIdxsMat(currLayer,2), cellIdxsMat(currLayer,3)) &
  !                                                                    % getPhiCapital(cellIdxsMat ,currLayer)

  ! end function getGridPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  !!
  !!
  !!
  function getNumberOfCells(self, localNxyz, n_layers, currLayer) result(report)
    class(cartesianGridIntermediate2), intent(in)        :: self
    integer(shortInt), dimension(:,:), intent(in)       :: localNxyz
    integer(shortInt), intent(in)                       :: n_layers, currLayer
    integer(shortInt), dimension(:), allocatable        :: output, report
    integer(shortInt)                                   :: i, j, k, l

    ! Allocate and initialise output array
    allocate(report(n_layers+1))
    report = 0

    do i = 1, localNxyz(currLayer-1,1)
      do j = 1, localNxyz(currLayer-1,2)
        do k = 1, localNxyz(currLayer-1,3)

          output = self % grid(i,j,k) % cellGetNumberOfCells(localNxyz, n_layers, currLayer)

          do l = 1, n_layers + 1
            report(l) = report(l) + output(l)
          end do

        end do
      end do
    end do

  end function getNumberOfCells

end module cartesianGridIntermediate2_class
