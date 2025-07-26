module cartesianGridSubLayer_inter

  use vertexShelf_class,               only : vertexShelf
  use edgeShelf_class,                 only : edgeShelf
  use faceShelf_class,                 only : faceShelf
  use elementShelf_class,              only : elementShelf
  use numPrecision

  implicit none
  private

! abstract base cartesianGrid data type for intermediate and finest layers
! using abstract data type allows switching between 2 and non-2 layered grids
  type, public, abstract                          :: cartesianGridSubLayer
  contains
    procedure(init), deferred                     :: init
    procedure(getGridChi), deferred               :: getGridChi
    procedure(getGridPhi), deferred               :: getGridPhi
    procedure(getGridPhiCapital), deferred        :: getGridPhiCapital
  end type cartesianGridSubLayer

  abstract interface

    !!
    !!
    !!
    pure subroutine init(self, vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
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
    pure function getGridChi(self, baseIntegerCoord, shift, mask, currLayer) result(chi)
      import                                                 cartesianGridSubLayer, shortInt
      class(cartesianGridSubLayer), intent(in)            :: self
      integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
      integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
      integer(shortInt), intent(in)                       :: currLayer
      integer(shortInt)                                   :: chi

    end function getGridChi

    !!
    !!
    !!
    pure function getGridPhi(self, baseIntegerCoord, shift, mask, currLayer) result(phi)
      import                                                 cartesianGridSubLayer, shortInt
      class(cartesianGridSubLayer), intent(in)            :: self
      integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
      integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
      integer(shortInt), intent(in)                       :: currLayer
      integer(shortInt)                                   :: phi

    end function getGridPhi

    !!
    !!
    !!
    pure function getGridPhiCapital(self, baseIntegerCoord, shift, mask, currLayer) result(phiCapital)
      import                                                 cartesianGridSubLayer, shortInt
      class(cartesianGridSubLayer), intent(in)            :: self
      integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
      integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
      integer(shortInt), intent(in)                       :: currLayer
      integer(shortInt)                                   :: phiCapital

    end function getGridPhiCapital

  end interface

end module cartesianGridSubLayer_inter
