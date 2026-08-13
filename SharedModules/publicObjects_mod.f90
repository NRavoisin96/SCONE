module publicObjects

  use numPrecision
  use ratint
  use RNG_class,         only : RNG
  use universalVariables

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
    type(ratint_t), dimension(3) :: ratintCoordinates  !!!!!!!!!!!!!!!!!!!!!!
  end type basicVertexInfo

  !!
  !!
  !!
  type :: coordData
    integer(shortInt)                     :: cellIdx = 0, elementIdx = 0, faceIdx = 0, front = 0, localId = 1, meshIdx = 0, &
                                             surfaceIdx = 0, universeIdx = 0, universeRootId = 0, updateLevel = 0
    integer(shortInt), dimension(VALENCE) :: currentFaceIdxs = 0
    logical(defBool)                      :: isInside = .false., isRotated = .false.
    real(defReal)                         :: d = INF, dMax = ZERO
    real(defReal), dimension(3)           :: r = ZERO, u = ZERO
    real(defReal), dimension(3, 3)        :: rotationMatrix = ZERO
  end type coordData

  !!
  !!
  !!
  type :: intersectionTestPayload
    integer(shortInt)                     :: front = 0
    integer(shortInt), dimension(VALENCE) :: currentFaceIdxs = 0
    real(defReal)                         :: dMax = ZERO
    real(defReal), dimension(3)           :: r = ZERO, u = ZERO
  end type intersectionTestPayload

  !!
  !!
  !!
  type :: rationalIntersectionTestPayload
    integer(shortInt)                     :: front = 0
    integer(shortInt), dimension(VALENCE) :: currentFaceIdxs = 0
    type(ratint_t)                        :: dMax
    type(ratint_t), dimension(3)          :: r, u
  end type rationalIntersectionTestPayload

  !!
  !!
  !!
  type :: intersectionTestResult
    logical(defBool) :: intersects = .false., needsRescue = .false.
    real(defReal)    :: d = INF
    type(ratint_t)   :: d_rational
  end type intersectionTestResult

  !!
  !!
  !!
  type :: meshBoundaryConditionInfo
    integer(shortInt), dimension(N_BC_TYPES)     :: boundaryConditions = [INTERNAL_TRANSPORT_BC, INTERNAL_TEMPERATURE_BC]
    integer(shortInt), dimension(:), allocatable :: faceIdxs
    real(defReal), dimension(N_BC_TYPES)         :: boundaryValues = ZERO
  end type meshBoundaryConditionInfo

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
    type(RNG), pointer :: rand => null()
  end type particleData

  !!
  !!
  !!
  type :: transportObjectStateCoordUpdateData
    integer(shortInt)           :: geometryIdx = 0, lowestCellIdx = 0, lowestElementIdx = 0, lowestMeshIdx = 0, &
                                   materialIdx = 0, uniqueId = 0
    real(defReal), dimension(3) :: rGlobal = ZERO, uGlobal = ZERO
  end type transportObjectStateCoordUpdateData

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
  pure function newRationalIntersectionTestPayload(r, u, dMax) result(payload)
    type(ratint_t), intent(in)               :: dMax
    type(ratint_t), dimension(3), intent(in) :: r, u
    type(rationalIntersectionTestPayload)    :: payload

    payload % r = r
    payload % u = u
    payload % dMax = dMax

  end function newRationalIntersectionTestPayload

  !!
  !!
  !!
  pure subroutine resetIntersectionTestResult(result)
    class(intersectionTestResult), intent(inout) :: result

    result % intersects = .false.
    result % d = INF

  end subroutine resetIntersectionTestResult

end module publicObjects