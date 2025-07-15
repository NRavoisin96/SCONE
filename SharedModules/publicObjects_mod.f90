module publicObjects

  use numPrecision
  use RNG_class,               only : RNG
  use universalVariables,      only : ZERO

  implicit none
  public

  !!
  !!
  !!
  type :: basicEdgeInfo
    integer(shortInt)               :: idx = 0
    integer(shortInt), dimension(2) :: vertexIdxs = 0
  end type basicEdgeInfo

  !!
  !!
  !!
  type :: basicElementInfo
    integer(shortInt)                            :: idx = 0, parentIdx = 0
    integer(shortInt), dimension(:), allocatable :: edgeIdxs, faceIdxs, vertexIdxs
  end type basicElementInfo

  !!
  !!
  !!
  type :: basicFaceInfo
    integer(shortInt)                            :: idx = 0, parentIdx = 0
    logical(defBool)                             :: isBoundary = .false.
    integer(shortInt), dimension(:), allocatable :: edgeIdxs, vertexIdxs
  end type basicFaceInfo

  !!
  !!
  !!
  type :: basicVertexInfo
    integer(shortInt)           :: idx = 0
    real(defReal), dimension(3) :: coordinates = ZERO
  end type basicVertexInfo

  !!
  !!
  !!
  type :: meshLocalIdInfo
    integer(shortInt)                            :: localId = 0
    integer(shortInt), dimension(:), allocatable :: elementIdxs
  end type meshLocalIdInfo

  !!
  !!
  !!
  type :: particleData
    integer(shortInt)   :: matIdx = 0
    real(defReal)       :: E = ZERO
    class(RNG), pointer :: rand => null()
  end type particleData

end module publicObjects