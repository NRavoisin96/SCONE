module cartesianCellIntermediate3_class
  
  use numPrecision
  use universalVariables,              only : ZERO
  use vertexShelf_class,               only : vertexShelf
  use edgeShelf_class,                 only : edgeShelf
  use faceShelf_class,                 only : faceShelf
  use cartesianInitProcedures
  use cartesianGridSubLayer_inter,     only : cartesianGridSubLayer
  use cartesianGridIntermediate4_class,only : cartesianGridIntermediate4
  use cartesianGridFinest_class,       only : cartesianGridFinest
  use genericProcedures,               only : append, fatalError

  implicit none
  private
  
  !!
  !!
  type, public                                          :: cartesianCellIntermediate3
    private
    class(cartesianGridSubLayer), pointer               :: subGrid => null()
    integer(shortInt)                                   :: chi = 0
    ! (needs to be changed) (cell centre should be a property to avoid repetative calc.)
    ! (due to limited memory, this is calculated each time needed)
    ! (try to avoid adding properties tho due to memory)

  contains

    ! (needs to be changed) (If we mark cells who lie outside of mesh in upper layers,)
    ! (we can say a huge memory and a bit of initialisation time and in-cycle time.)
    ! (see OneNote for more detail: Ideas>To do)

    ! Build procedures.
    procedure                                    :: cellTestPolyhedronInclusion
    procedure                                    :: refineCell
    ! Runtime procedures.
    procedure                                    :: getChi
    procedure                                    :: getPhi
    procedure                                    :: getPhiCapital
    ! Analysis procedures
    procedure                                    :: cellGetNumberOfCells
  end type cartesianCellIntermediate3

contains

  !!
  !!
  !!
  subroutine cellTestPolyhedronInclusion(self, faces, currElementFaceIdxs, centroid, &
                                              faceNormalSigns, elementIdx)
    class(cartesianCellIntermediate3), intent(inout)     :: self
    class(faceShelf), intent(in)                        :: faces
    integer(shortInt), dimension(:), intent(in)         :: currElementFaceIdxs
    real(defReal), dimension(3), intent(in)             :: centroid
    real(defReal), dimension(:,:), intent(in)           :: faceNormalSigns
    integer(shortInt), intent(in)                       :: elementIdx

    ! if current cell is found to be entirely contained within a polyhedron, there cannot not be other polyhedra
    if (self % chi /= 0) return

    ! Otherwise, continue testing and appending the array of candidate element indices 
    call testPolyhedronInclusion(faces, currElementFaceIdxs, centroid, faceNormalSigns, elementIdx, self % chi)

  end subroutine cellTestPolyhedronInclusion

  !!
  !!
  !!
  subroutine refineCell(self, vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                             currLayer, candidateElementIdxs, newGridBoundsMin, alpha, wStar)
    class(cartesianCellIntermediate3), intent(inout)     :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(inout)                     :: edges
    class(faceShelf), intent(inout)                     :: faces
    class(elementShelf), intent(in)                     :: elements
    real(defReal), dimension(:), intent(in)             :: spacing, spacingInv
    integer(shortInt), dimension(:,:), intent(in)       :: n_xyz
    integer(shortInt), intent(in)                       :: n_layers, currLayer
    integer(shortInt), dimension(:), intent(in)         :: candidateElementIdxs
    real(defReal), dimension(3), intent(in)             :: newGridBoundsMin
    real(defReal), intent(in)                           :: alpha, wStar

    !if the current cell is not entirely contained within a polyhedron, then refine the grid
    if (self % chi == 0) then

      ! (needs to be changed) (pointers cannot point to the same class)
      if (n_layers > currLayer+1) then 
        allocate(cartesianGridIntermediate4:: self % subGrid)
      else
        allocate(cartesianGridFinest:: self % subGrid)
      end if

      call self % subgrid % init(vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                                 currLayer + 1, candidateElementIdxs, newGridBoundsMin, alpha, wStar)

    end if 

  end subroutine refineCell

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! bit-trick (not saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! !!
  ! !!
  ! !!
  ! function getChi(self, baseIntegerCoord, shift, mask, currLayer) result(chi)
  !   class(cartesianCellIntermediate3), intent(in)        :: self
  !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
  !   integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: chi

  !   if (self % chi /= 0) then
  !     chi = self % chi
  !   else
  !     chi = self % subGrid % getGridChi(baseIntegerCoord, shift, mask, currLayer+1)
  !   end if

  ! end function getChi

  ! !!
  ! !!
  ! !!
  ! function getPhi(self, baseIntegerCoord, shift, mask, currLayer) result(phi)
  !   class(cartesianCellIntermediate3), intent(in)        :: self
  !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
  !   integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: phi

  !   phi = self % subGrid % getGridPhi(baseIntegerCoord, shift, mask, currLayer+1)

  ! end function getPhi

  ! !!
  ! !!
  ! !!
  ! function getPhiCapital(self, baseIntegerCoord, shift, mask, currLayer) result(phiCapital)
  !   class(cartesianCellIntermediate3), intent(in)        :: self
  !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
  !   integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: phiCapital

  !   phiCapital = self % subGrid % getGridPhiCapital(baseIntegerCoord, shift, mask, currLayer+1)

  ! end function getPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! bit-trick (saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !!
  !!
  !!
  function getChi(self, baseIntegerCoord, shift, mask, currLayer, cellIdxsMat) result(chi)
    class(cartesianCellIntermediate3), intent(in)        :: self
    integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
    integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
    integer(shortInt), dimension(:,:), intent(inout)    :: cellIdxsMat
    integer(shortInt), intent(in)                       :: currLayer
    integer(shortInt)                                   :: chi

    if (self % chi /= 0) then
      chi = self % chi
    else
      chi = self % subGrid % getGridChi(baseIntegerCoord, shift, mask, currLayer+1, cellIdxsMat)
    end if

  end function getChi

  !!
  !!
  !!
  function getPhi(self, cellIdxsMat, currLayer) result(phi)
    class(cartesianCellIntermediate3), intent(in)        :: self
    integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
    integer(shortInt), intent(in)                       :: currLayer
    integer(shortInt)                                   :: phi

    phi = self % subGrid % getGridPhi(cellIdxsMat, currLayer+1)

  end function getPhi

  !!
  !!
  !!
  function getPhiCapital(self, cellIdxsMat, currLayer) result(phiCapital)
    class(cartesianCellIntermediate3), intent(in)        :: self
    integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
    integer(shortInt), intent(in)                       :: currLayer
    integer(shortInt)                                   :: phiCapital

    phiCapital = self % subGrid % getGridPhiCapital(cellIdxsMat, currLayer+1)

  end function getPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Non bit-trick (saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! !!
  ! !!
  ! !!
  ! function getChi(self, r, gridBounds_min, spacingInv, currLayer, cellIdxsMat, nSub_xyz) result(chi)
  !   class(cartesianCellIntermediate3), intent(in)        :: self
  !   real(defReal), dimension(3), intent(in)             :: r, gridBounds_min                                    
  !   real(defReal), dimension(:), intent(in)             :: spacingInv
  !   integer(shortInt), dimension(:,:), intent(inout)    :: cellIdxsMat
  !   integer(shortInt), dimension(:,:), intent(in)       :: nSub_xyz
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: chi

  !   if (self % chi /= 0) then
  !     chi = self % chi
  !   else
  !     chi = self % subGrid % getGridChi(r, gridBounds_min, spacingInv, currLayer+1, cellIdxsMat, nSub_xyz)
  !   end if

  ! end function getChi

  ! !!
  ! !!
  ! !!
  ! function getPhi(self, cellIdxsMat, currLayer) result(phi)
  !   class(cartesianCellIntermediate3), intent(in)        :: self
  !   integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: phi

  !   phi = self % subGrid % getGridPhi(cellIdxsMat, currLayer+1)

  ! end function getPhi

  ! !!
  ! !!
  ! !!
  ! function getPhiCapital(self, cellIdxsMat, currLayer) result(phiCapital)
  !   class(cartesianCellIntermediate3), intent(in)        :: self
  !   integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
  !   integer(shortInt), intent(in)                       :: currLayer
  !   integer(shortInt)                                   :: phiCapital

  !   phiCapital = self % subGrid % getGridPhiCapital(cellIdxsMat, currLayer+1)

  ! end function getPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  !!
  !!
  !!
  function cellGetNumberOfCells(self, localNxyz, n_layers, currLayer) result(output)
    class(cartesianCellIntermediate3), intent(in)        :: self
    integer(shortInt), dimension(:,:), intent(in)       :: localNxyz
    integer(shortInt), intent(in)                       :: n_layers, currLayer
    integer(shortInt), dimension(:), allocatable        :: output

    ! Allocate and initialise output array
    allocate(output(n_layers+1))
    output = 0

    if (self % chi /= 0) then
      output(currLayer) = 1
    else
      output = self % subgrid % getNumberOfCells(localNxyz, n_layers, currLayer + 1)
    end if
    
  end function cellGetNumberOfCells


end module cartesianCellIntermediate3_class