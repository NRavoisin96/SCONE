module cartesianGridSubLayer_inter

  use vertexShelf_class,               only : vertexShelf
  use edgeShelf_class,                 only : edgeShelf
  use faceShelf_class,                 only : faceShelf
  use elementShelf_class,              only : elementShelf
  use numPrecision,                    only : shortInt, defReal
  use ragged3dMatrix_class,            only : ragged3d
  use dynamic2dMatSet_class,           only : dynamic2dMatSet

  implicit none
  private

! abstract base cartesianGrid data type for intermediate and finest layers
! using abstract data type allows switching between 2 and non-2 layered grids
  type, public, abstract                          :: cartesianGridSubLayer
  contains
    procedure(init), deferred                     :: init
    procedure(init2), deferred                    :: init2
    procedure(initt), deferred                    :: initt
    procedure(getGridChi), deferred               :: getGridChi
    procedure(getGridPhi), deferred               :: getGridPhi
    procedure(getGridPhiCapital), deferred        :: getGridPhiCapital
    procedure(getNumberOfCells), deferred         :: getNumberOfCells
  end type cartesianGridSubLayer

  abstract interface

    !!
    !!
    !!
    subroutine init(self, vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                         currLayer, candidateElementIdxs, gridBoundsMin, alpha, wStar)
      import                                                 cartesianGridSubLayer, vertexShelf, edgeShelf, faceShelf, &
                                                             elementShelf, defReal, shortInt
      class(cartesianGridSubLayer), intent(inout)         :: self
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

     end subroutine

    !!
    !!
    !!
    subroutine init2(self, vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                         currLayer, candidateElementIdxs, gridBoundsMin, alpha, wStar, normalSignsMat, &
                         circumscribedBallRadius, targetDistance, targetDistanceSqr)
      import                                                 cartesianGridSubLayer, vertexShelf, edgeShelf, faceShelf, &
                                                             elementShelf, defReal, shortInt, dynamic2dMatSet
      class(cartesianGridSubLayer), intent(inout)         :: self
      class(vertexShelf), intent(in)                      :: vertices
      class(edgeShelf), intent(inout)                     :: edges
      class(faceShelf), intent(inout)                     :: faces
      class(elementShelf), intent(in)                     :: elements
      real(defReal), dimension(:), intent(in)             :: spacing, spacingInv
      integer(shortInt), dimension(:,:), intent(in)       :: n_xyz
      integer(shortInt), intent(in)                       :: n_layers, currLayer
      integer(shortInt), dimension(:), intent(in)         :: candidateElementIdxs
      real(defReal), dimension(3), intent(in)             :: gridBoundsMin
      real(defReal), intent(in)                           :: alpha, wStar, circumscribedBallRadius, targetDistance, &
                                                             targetDistanceSqr
      type(dynamic2dMatSet), intent(in)                   :: normalSignsMat

     end subroutine 

  !!
  !!
  !!
  subroutine initt(self, vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                   currLayer, intersectedFaceIdxs, gridBoundsMin, alpha, wStar, &
                   extraDistanceArr, candidateElementIdxs, normalSignsMat, &
                   circumscribedBallRadius, targetDistance, targetDistanceSqr)
      import                                               cartesianGridSubLayer, vertexShelf, edgeShelf, faceShelf, &
                                                           elementShelf, defReal, shortInt, ragged3d
    class(cartesianGridSubLayer), intent(inout)         :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(inout)                     :: edges
    class(faceShelf), intent(inout)                     :: faces
    class(elementShelf), intent(in)                     :: elements
    real(defReal), dimension(:), intent(in)             :: spacing, spacingInv
    integer(shortInt), dimension(:,:), intent(in)       :: n_xyz
    integer(shortInt), intent(in)                       :: n_layers, currLayer
    integer(shortInt), dimension(:), intent(in)         :: intersectedFaceIdxs, candidateElementIdxs
    real(defReal), dimension(3), intent(in)             :: gridBoundsMin
    real(defReal), intent(in)                           :: alpha, wStar, circumscribedBallRadius, targetDistance, &
                                                           targetDistanceSqr
    real(defReal), dimension(:), intent(in)             :: extraDistanceArr
    type(ragged3d), intent(in)                          :: normalSignsMat

  end subroutine initt

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! bit-trick (not saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    ! !!
    ! !!
    ! !!
    ! function getGridChi(self, baseIntegerCoord, shift, mask, currLayer) result(chi)
    !   import                                                 cartesianGridSubLayer, shortInt
    !   class(cartesianGridSubLayer), intent(in)            :: self
    !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
    !   integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
    !   integer(shortInt), intent(in)                       :: currLayer
    !   integer(shortInt)                                   :: chi

    ! end function getGridChi

    ! !!
    ! !!
    ! !!
    ! function getGridPhi(self, baseIntegerCoord, shift, mask, currLayer) result(phi)
    !   import                                                 cartesianGridSubLayer, shortInt
    !   class(cartesianGridSubLayer), intent(in)            :: self
    !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
    !   integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
    !   integer(shortInt), intent(in)                       :: currLayer
    !   integer(shortInt)                                   :: phi

    ! end function getGridPhi

    ! !!
    ! !!
    ! !!
    ! function getGridPhiCapital(self, baseIntegerCoord, shift, mask, currLayer) result(phiCapital)
    !   import                                                 cartesianGridSubLayer, shortInt
    !   class(cartesianGridSubLayer), intent(in)            :: self
    !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
    !   integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
    !   integer(shortInt), intent(in)                       :: currLayer
    !   integer(shortInt)                                   :: phiCapital

    ! end function getGridPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! bit-trick (saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    ! !!
    ! !!
    ! !!
    ! function getGridChi(self, baseIntegerCoord, shift, mask, currLayer, cellIdxsMat) result(chi)
    !   import                                                 cartesianGridSubLayer, shortInt
    !   class(cartesianGridSubLayer), intent(in)            :: self
    !   integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
    !   integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
    !   integer(shortInt), dimension(:,:), intent(inout)    :: cellIdxsMat
    !   integer(shortInt), intent(in)                       :: currLayer
    !   integer(shortInt)                                   :: chi

    ! end function getGridChi

    ! !!
    ! !!
    ! !!
    ! function getGridPhi(self, cellIdxsMat, currLayer) result(phi)
    !   import                                                 cartesianGridSubLayer, shortInt
    !   class(cartesianGridSubLayer), intent(in)            :: self
    !   integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
    !   integer(shortInt), intent(in)                       :: currLayer
    !   integer(shortInt)                                   :: phi

    ! end function getGridPhi

    ! !!
    ! !!
    ! !!
    ! function getGridPhiCapital(self, cellIdxsMat, currLayer) result(phiCapital)
    !   import                                                 cartesianGridSubLayer, shortInt
    !   class(cartesianGridSubLayer), intent(in)            :: self
    !   integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
    !   integer(shortInt), intent(in)                       :: currLayer
    !   integer(shortInt)                                   :: phiCapital

    ! end function getGridPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Non bit-trick (saving)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    !!
    !!
    !!
    function getGridChi(self, r, gridBounds_min, spacingInv, currLayer, cellIdxsMat, nSub_xyz) result(chi)
      import                                                 cartesianGridSubLayer, shortInt, defReal
      class(cartesianGridSubLayer), intent(in)            :: self
      real(defReal), dimension(3), intent(in)             :: r, gridBounds_min                                    
      real(defReal), dimension(:), intent(in)             :: spacingInv
      integer(shortInt), dimension(:,:), intent(inout)    :: cellIdxsMat
      integer(shortInt), dimension(:,:), intent(in)       :: nSub_xyz
      integer(shortInt), intent(in)                       :: currLayer
      integer(shortInt)                                   :: chi

    end function getGridChi

    !!
    !!
    !!
    function getGridPhi(self, cellIdxsMat, currLayer) result(phi)
      import                                                 cartesianGridSubLayer, shortInt
      class(cartesianGridSubLayer), intent(in)            :: self
      integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
      integer(shortInt), intent(in)                       :: currLayer
      integer(shortInt)                                   :: phi

    end function getGridPhi

    !!
    !!
    !!
    function getGridPhiCapital(self, cellIdxsMat, currLayer) result(phiCapital)
      import                                                 cartesianGridSubLayer, shortInt
      class(cartesianGridSubLayer), intent(in)            :: self
      integer(shortInt), dimension(:,:), intent(in)       :: cellIdxsMat
      integer(shortInt), intent(in)                       :: currLayer
      integer(shortInt)                                   :: phiCapital

    end function getGridPhiCapital

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    !!
    !!
    !!
    function getNumberOfCells(self, localNxyz, n_layers, currLayer) result(report)
      import                                                 cartesianGridSubLayer, shortInt
      class(cartesianGridSubLayer), intent(in)            :: self
      integer(shortInt), dimension(:,:), intent(in)       :: localNxyz
      integer(shortInt), intent(in)                       :: n_layers, currLayer
      integer(shortInt), dimension(:), allocatable        :: output, report
      integer(shortInt)                                   :: i, j, k, l

    end function getNumberOfCells

  end interface

end module cartesianGridSubLayer_inter
