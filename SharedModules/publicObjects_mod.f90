module publicObjects

  use edge_class,              only : edgeBox
  use face_class,              only : orientatedFaceBox
  use numPrecision
  use RNG_class,               only : RNG
  use topologicalObject_inter, only : topologicalObject
  use universalVariables,      only : ZERO
  use vertex_class,            only : vertexBox

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
  type :: buildEdgeInfo
    integer(shortInt)             :: idx = 0
    type(vertexBox), dimension(2) :: vertices
  end type buildEdgeInfo

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