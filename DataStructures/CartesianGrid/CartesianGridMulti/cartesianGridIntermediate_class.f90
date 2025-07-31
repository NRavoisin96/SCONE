module cartesianGridIntermediate_class

  use vertexShelf_class,                       only : vertexShelf
  use edgeShelf_class,                         only : edgeShelf
  use faceShelf_class,                         only : faceShelf
  use elementShelf_class,                      only : elementShelf
  use cartesianGridSubLayer_inter,             only : cartesianGridSubLayer
  use cartesianCellIntermediate_class,         only : cartesianCellIntermediate
  use numPrecision   
  use cartesianGenericProcedures

  implicit none
  private


  !!
  !!
  !!
  type, public, extends(cartesianGridSubLayer)                      :: cartesianGridIntermediate
    private
    type(cartesianCellIntermediate), dimension(:,:,:), allocatable  :: grid

  contains

    ! Build procedures
    procedure                                    :: init
    procedure                                    :: constructMapping
    procedure                                    :: refineGrid
    ! Runtime procedures
    procedure                                    :: getGridChi
    procedure                                    :: getGridPhi
    procedure                                    :: getGridPhiCapital

  end type cartesianGridIntermediate

contains

  !!
  !!
  !!
  subroutine init(self, vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                       currLayer, candidateElementIdxs, gridBoundsMin, alpha, wStar)
    class(cartesianGridIntermediate), intent(inout)     :: self
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
    class(cartesianGridIntermediate), intent(inout)     :: self
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
    class(cartesianGridIntermediate), intent(inout)     :: self
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

  !!
  !!
  !! 
  function getGridChi(self, baseIntegerCoord, shift, mask, currLayer) result(chi)
    class(cartesianGridIntermediate), intent(in)        :: self
    integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
    integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
    integer(shortInt), intent(in)                       :: currLayer
    integer(shortInt)                                   :: chi
    integer(shortInt), dimension(3)                     :: cellIdxs

    cellIdxs = getLocalIdx(baseIntegerCoord, shift(currLayer,:), mask(currLayer,:))

    chi = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getChi(baseIntegerCoord, &
                                                                    shift,mask,currLayer)

  end function getGridChi

  !!
  !!
  !! 
  function getGridPhi(self, baseIntegerCoord, shift, mask, currLayer) result(phi)
    class(cartesianGridIntermediate), intent(in)        :: self
    integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
    integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
    integer(shortInt), intent(in)                       :: currLayer
    integer(shortInt)                                   :: phi
    integer(shortInt), dimension(3)                     :: cellIdxs

    cellIdxs = getLocalIdx(baseIntegerCoord, shift(currLayer,:), mask(currLayer,:))

    phi = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getPhi(baseIntegerCoord, &
                                                                    shift,mask,currLayer)

  end function getGridPhi

  !!
  !!
  !! 
  function getGridPhiCapital(self, baseIntegerCoord, shift, mask, currLayer) result(phiCapital)
    class(cartesianGridIntermediate), intent(in)        :: self
    integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
    integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
    integer(shortInt), intent(in)                       :: currLayer
    integer(shortInt)                                   :: phiCapital
    integer(shortInt), dimension(3)                     :: cellIdxs

    cellIdxs = getLocalIdx(baseIntegerCoord, shift(currLayer,:), mask(currLayer,:))

    phiCapital = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getPhiCapital(&
                                              baseIntegerCoord, shift,mask,currLayer)

  end function getGridPhiCapital

end module cartesianGridIntermediate_class
