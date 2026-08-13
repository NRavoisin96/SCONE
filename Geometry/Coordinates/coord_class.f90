module coord_class

  use genericProcedures,  only : areEqual, numToChar
  use numPrecision
  use publicObjects,      only : coordData
  use universalVariables, only : HARDCODED_MAX_NEST, NUDGE, VALENCE

  implicit none
  private

  

  !!
  !! Co-ordinates in a single geometry level
  !!
  !! Co-ordinates are considered valid if:
  !!   * dir is normalised to 1.0 (norm2(dir) ~= 1.0)
  !!   * uniIdx, uniRootId & localId are set to +ve values
  !!
  !! Public Members:
  !!   r          -> Position
  !!   rEnd       -> Pre-computed end position (reserved for mesh tracking)
  !!   dir        -> Direction
  !!   isRotated  -> Is rotated wrt previous (higher by 1) level
  !!   rotMat     -> Rotation matrix wrt previous level
  !!   uniIdx     -> Index of the occupied universe
  !!   uniRootId  -> Location of the occupied universe in geometry graph
  !!   localId    -> Local cell in the occupied universe
  !!   cellIdx    -> Index of the occupied cell in cellShelf. 0 is cell is local to the universe
  !!   elementIdx -> Index of the occupied element in meshShelf (for mesh universe).
  !!
  !! Interface:
  !!   isValid    -> Returns .true. if coordinates are valid
  !!   display    -> Prints coordinates to the console
  !!   kill       -> Returns to uninitialised state
  !!
  type, public :: coord
    private
    integer(shortInt)                     :: cellIdx = 0, elementIdx = 0, front = 0, localId = 0, meshIdx = 0, &
                                             universeIdx = 0, universeRootId = 0
    integer(shortInt), dimension(VALENCE) :: currentFaceIdxs = 0
    logical(defBool)                      :: isRotated = .false.
    real(defReal), dimension(3)           :: r = ZERO, dir = ZERO
    real(defReal), dimension(3, 3)        :: rotMat = ZERO
  contains
    procedure :: display
    procedure :: getCellIdx
    procedure :: getData
    procedure :: getDirection
    procedure :: getElementIdx
    procedure :: getIsRotated
    procedure :: getLocalId
    procedure :: getMeshIdx
    procedure :: getPosition
    procedure :: getRotationMatrix
    procedure :: getUniverseIdx
    procedure :: getUniverseRootId
    procedure :: isValid
    procedure :: kill
    procedure :: offsetPosition
    procedure :: rotateComponents
    procedure :: setCellIdx
    procedure :: setDirection
    procedure :: setElementIdx
    procedure :: setIsRotated
    procedure :: setLocalId
    procedure :: setMeshIdx
    procedure :: setPosition
    procedure :: setPositionAndDirection
    procedure :: setRotationMatrix
    procedure :: setUniverseIdx
    procedure :: setUniverseRootId
    procedure :: updateFromData
  end type coord

contains
  !!
  !! Print to screen contents of the coord
  !!
  subroutine display(self)
    class(coord), intent(in) :: self

    print *, 'Position: ', self % r
    print *, 'Direction: ', self % dir
    print *, 'Universe index: ', numToChar(self % universeIdx)
    print *, 'Local id: ', numToChar(self % localId)
    print *, 'Universe root id: ', numToChar(self % universeRootId)

  end subroutine display

  !!
  !!
  !!
  elemental function getCellIdx(self) result(cellIdx)
    class(coord), intent(in) :: self
    integer(shortInt)        :: cellIdx

    cellIdx = self % cellIdx

  end function getCellIdx

  !!
  !!
  !!
  elemental function getData(self) result(data)
    class(coord), intent(in) :: self
    type(coordData)          :: data

    data % cellIdx = self % cellIdx
    data % elementIdx = self % elementIdx
    data % front = self % front
    data % localId = self % localId
    data % meshIdx = self % meshIdx
    data % universeIdx = self % universeIdx
    data % universeRootId = self % universeRootId
    data % currentFaceIdxs = self % currentFaceIdxs
    data % r = self % r
    data % u = self % dir

  end function getData

  !! Function 'getDirection'
  !!
  !! Basic description:
  !!   Returns the direction of the coordinates.
  !!
  !! Result:
  !!   u -> 3D direction of the coordinates.
  !!
  pure function getDirection(self) result(u)
    class(coord), intent(in)    :: self
    real(defReal), dimension(3) :: u

    u = self % dir

  end function getDirection

  !!
  !!
  !!
  elemental function getElementIdx(self) result(elementIdx)
    class(coord), intent(in) :: self
    integer(shortInt)        :: elementIdx

    elementIdx = self % elementIdx

  end function getElementIdx

  !!
  !!
  !!
  elemental function getIsRotated(self) result(isRotated)
    class(coord), intent(in) :: self
    logical(defBool)         :: isRotated

    isRotated = self % isRotated

  end function getIsRotated

  !!
  !!
  !!
  elemental function getLocalId(self) result(localId)
    class(coord), intent(in) :: self
    integer(shortInt)        :: localId

    localId = self % localId

  end function getLocalId

  !!
  !!
  !!
  elemental function getMeshIdx(self) result(meshIdx)
    class(coord), intent(in) :: self
    integer(shortInt)        :: meshIdx

    meshIdx = self % meshIdx

  end function getMeshIdx

  !! Function 'getPosition'
  !!
  !! Basic description:
  !!   Returns the position of the coordinates.
  !!
  !! Result:
  !!   r -> 3D position of the coordinates.
  !!
  pure function getPosition(self) result(r)
    class(coord), intent(in)    :: self
    real(defReal), dimension(3) :: r

    r = self % r

  end function getPosition

  !!
  !!
  !!
  pure function getRotationMatrix(self) result(rotationMatrix)
    class(coord), intent(in)       :: self
    real(defReal), dimension(3, 3) :: rotationMatrix

    rotationMatrix = self % rotMat

  end function getRotationMatrix

  !!
  !!
  !!
  elemental function getUniverseIdx(self) result(universeIdx)
    class(coord), intent(in) :: self
    integer(shortInt)        :: universeIdx

    universeIdx = self % universeIdx

  end function getUniverseIdx

  !!
  !!
  !!
  elemental function getUniverseRootId(self) result(universeRootId)
    class(coord), intent(in) :: self
    integer(shortInt)        :: universeRootId

    universeRootId = self % universeRootId

  end function getUniverseRootId

  !!
  !! Returns .true. if coordinates are valid
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   True if coord is valid. See type doc-comment for definition of valid.
  !!
  elemental function isValid(self) result(correct)
    class(coord), intent(in) :: self
    logical(defBool)         :: correct

    ! Direction vector is normalised within floating point tolerance
    correct = areEqual(norm2(self % dir), ONE)
    correct = correct .and. 0 < self % localId
    correct = correct .and. 0 < self % universeIdx
    correct = correct .and. 0 < self % universeRootId

  end function isValid

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(coord), intent(inout) :: self

    self % cellIdx = 0
    self % dir = ZERO
    self % elementIdx = 0
    self % isRotated = .false.
    self % localId = 0
    self % meshIdx = 0
    self % universeIdx = 0
    self % universeRootId = 0
    self % r = ZERO
    self % rotMat = ZERO

  end subroutine kill

  !!
  !!
  !!
  pure subroutine offsetPosition(self, offset)
    class(coord), intent(inout)             :: self
    real(defReal), dimension(3), intent(in) :: offset

    self % r = self % r - offset

  end subroutine offsetPosition

  !!
  !!
  !!
  elemental subroutine rotateComponents(self)
    class(coord), intent(inout) :: self

    self % r = matmul(self % rotMat, self % r)
    self % dir = matmul(self % rotMat, self % dir)

  end subroutine rotateComponents

  !!
  !!
  !!
  elemental subroutine setCellIdx(self, cellIdx)
    class(coord), intent(inout)   :: self
    integer(shortInt), intent(in) :: cellIdx

    self % cellIdx = cellIdx

  end subroutine setCellIdx

  !!
  !!
  !!
  pure subroutine setDirection(self, u)
    class(coord), intent(inout)             :: self
    real(defReal), dimension(3), intent(in) :: u

    self % dir = u

  end subroutine setDirection

  !!
  !!
  !!
  elemental subroutine setElementIdx(self, elementIdx)
    class(coord), intent(inout)   :: self
    integer(shortInt), intent(in) :: elementIdx

    self % elementIdx = elementIdx

  end subroutine setElementIdx

  !!
  !!
  !!
  elemental subroutine setIsRotated(self, isRotated)
    class(coord), intent(inout)  :: self
    logical(defBool), intent(in) :: isRotated

    self % isRotated = isRotated

  end subroutine setIsRotated

  !!
  !!
  !!
  elemental subroutine setLocalId(self, localId)
    class(coord), intent(inout)   :: self
    integer(shortInt), intent(in) :: localId

    self % localId = localId

  end subroutine setLocalId

  !!
  !!
  !!
  elemental subroutine setMeshIdx(self, meshIdx)
    class(coord), intent(inout)   :: self
    integer(shortInt), intent(in) :: meshIdx

    self % meshIdx = meshIdx

  end subroutine setMeshIdx

  !!
  !!
  !!
  pure subroutine setPosition(self, r)
    class(coord), intent(inout)             :: self
    real(defReal), dimension(3), intent(in) :: r

    self % r = r

  end subroutine setPosition

  !!
  !!
  !!
  pure subroutine setPositionAndDirection(self, r, u)
    class(coord), intent(inout)             :: self
    real(defReal), dimension(3), intent(in) :: r, u

    self % r = r
    self % dir = u

  end subroutine setPositionAndDirection

  !!
  !!
  !!
  pure subroutine setRotationMatrix(self, rotationMatrix)
    class(coord), intent(inout)                :: self
    real(defReal), dimension(3, 3), intent(in) :: rotationMatrix

    self % rotMat = rotationMatrix

  end subroutine setRotationMatrix

  !!
  !!
  !!
  elemental subroutine setUniverseIdx(self, universeIdx)
    class(coord), intent(inout)   :: self
    integer(shortInt), intent(in) :: universeIdx

    self % universeIdx = universeIdx

  end subroutine setUniverseIdx

  !!
  !!
  !!
  elemental subroutine setUniverseRootId(self, universeRootId)
    class(coord), intent(inout)   :: self
    integer(shortInt), intent(in) :: universeRootId

    self % universeRootId = universeRootId

  end subroutine setUniverseRootId

  !!
  !!
  !!
  elemental subroutine updateFromData(self, data)
    class(coord), intent(inout) :: self
    type(coordData), intent(in) :: data

    self % cellIdx = data % cellIdx
    self % elementIdx = data % elementIdx
    self % front = data % front
    self % localId = data % localId
    self % meshIdx = data % meshIdx
    self % universeIdx = data % universeIdx
    self % universeRootId = data % universeRootId
    self % currentFaceIdxs = data % currentFaceIdxs
    self % isRotated = data % isRotated
    self % r = data % r
    self % dir = data % u
    self % rotMat = data % rotationMatrix

  end subroutine updateFromData

end module coord_class