module cartesianCellCoarsest_class
  
  use numPrecision
  use universalVariables,              only : ZERO
  use vertexShelf_class,               only : vertexShelf
  use edgeShelf_class,                 only : edgeShelf
  use faceShelf_class,                 only : faceShelf
  use cartesianInitProcedures
  use cartesianGridSubLayer_inter,     only : cartesianGridSubLayer
  use cartesianGridFinest_class,       only : cartesianGridFinest
  use cartesianGridIntermediate_class, only : cartesianGridIntermediate
  use genericProcedures,               only : append, fatalError

  implicit none
  private
  
  !!
  !!
  type, public                                          :: cartesianCellCoarsest
    private
    class(cartesianGridSubLayer), pointer               :: subGrid => null()
    integer(shortInt)                                   :: chi = 0
    integer(shortInt), dimension(:), allocatable        :: candidateElementIdxs
    ! (needs to be changed) (cell centre should be a property to avoid repetative calc.)
    ! (due to limited memory, this is calculated each time needed)
    ! (try to avoid adding properties tho due to memory)

  contains

    ! Build procedures.
    procedure                                    :: cellTestPolyhedronInclusion
    procedure                                    :: setIsOutsideMesh
    procedure                                    :: refineCell
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
    if (self % chi /= 0) return

    ! Otherwise, continue testing and appending the array of candidate element indices 
    call testPolyhedronInclusion(faces, currElementFaceIdxs, centroid, faceNormalSigns, elementIdx, self % chi)
    call append(self % candidateElementIdxs, elementIdx)

  end subroutine cellTestPolyhedronInclusion


  !!
  !!
  !!
  subroutine setIsOutsideMesh(self)
    class(cartesianCellCoarsest), intent(inout)         :: self

    ! if candidateElementIdxs is not allocated, this cell does not intersect AABB of any polyhedron.
    ! Hence, this cell lies outside of mesh
    if (.NOT. allocated(self % candidateElementIdxs)) self % chi = -1

  end subroutine setIsOutsideMesh

  !!
  !!
  !!
  subroutine refineCell(self, vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                             newGridBoundsMin, alpha, wStar)
    class(cartesianCellCoarsest), intent(inout)         :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(inout)                     :: edges
    class(faceShelf), intent(inout)                     :: faces
    class(elementShelf), intent(in)                     :: elements
    real(defReal), dimension(:), intent(in)             :: spacing, spacingInv
    integer(shortInt), dimension(:,:), intent(in)       :: n_xyz
    integer(shortInt), intent(in)                       :: n_layers
    real(defReal), dimension(3), intent(in)             :: newGridBoundsMin
    real(defReal), intent(in)                           :: alpha, wStar

    !if the current cell is not entirely contained within a polyhedron, then refine the grid
    if (self % chi == 0) then

      if (n_layers > 2) then 
        allocate(cartesianGridIntermediate:: self % subGrid)
      else
        allocate(cartesianGridFinest:: self % subGrid)
      end if

      call self % subgrid % init(vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                                 2, self % candidateElementIdxs, newGridBoundsMin, alpha, wStar)

    end if 

    ! deallocate candidateElementIdxs to save memory. Otherwise, memory could explode
    if (allocated(self % candidateElementIdxs)) deallocate(self % candidateElementIdxs)

  end subroutine refineCell

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