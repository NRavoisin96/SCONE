module publicObjects

  use numPrecision
  use RNG_class,               only : RNG
  use universalVariables,      only : INF, ZERO

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
  type :: coordData
    real(defReal)                  :: d = INF, dMax = ZERO
    real(defReal), dimension(3)    :: r = ZERO, u = ZERO
    real(defReal), dimension(3, 3) :: rotationMatrix = ZERO
    integer(shortInt)              :: cellIdx = 0, elementIdx = 0, faceIdx = 0, localId = 1, meshIdx = 0, &
                                      surfaceIdx = 0, universeIdx = 0, universeRootId = 0, updateLevel = 0
    logical(defBool)               :: isInside = .false., isRotated = .false.
  end type coordData

  !!
  !!
  !!
  type :: intersectionTestPayload
    real(defReal), dimension(3) :: r = ZERO, u = ZERO
    real(defReal)               :: dMax = ZERO
  end type intersectionTestPayload

  !!
  !!
  !!
  type :: intersectionTestResult
    logical(defBool) :: intersects = .false.
    real(defReal)    :: d = INF
  end type intersectionTestResult

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

contains
  !!
  !!
  !!
  pure function newCoordData(r, u, cellIdx, localId, surfaceIdx, universeIdx, universeRootId, dMax) result(data)
    real(defReal), dimension(3), intent(in) :: r, u
    integer(shortInt), intent(in), optional :: cellIdx, localId, surfaceIdx, universeIdx, universeRootId
    real(defReal), intent(in), optional     :: dMax
    type(coordData)                         :: data

    data % r = r
    data % u = u / norm2(u)
    if (present(cellIdx)) data % cellIdx = cellIdx
    if (present(localId)) data % localId = localId
    if (present(surfaceIdx)) data % surfaceIdx = surfaceIdx
    if (present(universeIdx)) data % universeIdx = universeIdx
    if (present(universeRootId)) data % universeRootId = universeRootId
    if (present(dMax)) data % dMax = dMax

  end function newCoordData

  !!
  !!
  !!
  pure function newIntersectionTestPayload(r, u, dMax) result(payload)
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal), intent(in)               :: dMax
    type(intersectionTestPayload)           :: payload

    payload % r = r
    payload % u = u
    payload % dMax = dMax

  end function newIntersectionTestPayload

  !!
  !!
  !!
  pure subroutine resetIntersectionTestResult(result)
    class(intersectionTestResult), intent(inout) :: result

    result % intersects = .false.
    result % d = INF

  end subroutine resetIntersectionTestResult

end module publicObjects