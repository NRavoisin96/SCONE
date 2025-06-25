module coord_class

  use numPrecision
  use universalVariables, only : HARDCODED_MAX_NEST, NUDGE
  use genericProcedures,  only : areEqual, numToChar

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
  !!   isLeaving  -> Is leaving 
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
  type, public                    :: coord
    private
    real(defReal), dimension(3)   :: r = ZERO, rEnd = ZERO, dir = ZERO
    logical(defBool)              :: isRotated = .false., isLeaving = .false., nudgeEndPosition = .false.
    real(defReal), dimension(3,3) :: rotMat = ZERO
    integer(shortInt)             :: uniIdx = 0, uniRootId = 0, localId = 0, cellIdx = 0, elementIdx = 0, &
                                     parentElementIdx = 0
  contains
    procedure :: display
    procedure :: getCellIdx
    procedure :: getDirection
    procedure :: getElementIdx
    procedure :: getEndPosition
    procedure :: getIsRotated
    procedure :: getLocalId
    procedure :: getParentElementIdx
    procedure :: getPosition
    procedure :: getPositionToNudge
    procedure :: getRotationMatrix
    procedure :: getUniIdx
    procedure :: getUniRootId
    procedure :: isValid
    procedure :: kill
    procedure :: nudgePosition
    procedure :: offsetPosition
    procedure :: rotateComponents
    procedure :: setCellIdx
    procedure :: setDirection
    procedure :: setElementIdx
    procedure :: setEndPosition
    procedure :: setIsRotated
    procedure :: setLocalId
    procedure :: setNudgeEndPosition
    procedure :: setParentElementIdx
    procedure :: setPosition
    procedure :: setRotationMatrix
    procedure :: setUniIdx
    procedure :: setUniRootId
  end type coord

contains
  !!
  !! Print to screen contents of the coord
  !!
  subroutine display(self)
    class(coord), intent(in) :: self

    print *, "R: ", self % r
    print *, "U: ", self % dir
    print *, "UniIdx: ", numToChar(self % uniIdx), " LocalId: ", numToChar(self % localId), &
             "UniRootId", numToChar(self % uniRootId)

  end subroutine display

  !!
  !!
  !!
  elemental function getCellIdx(self) result(cellIdx)
    class(coord), intent(in) :: self
    integer(shortInt)        :: cellIdx

    cellIdx = self % cellIdx

  end function getCellIdx

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
  pure function getEndPosition(self) result(rEnd)
    class(coord), intent(in)    :: self
    real(defReal), dimension(3) :: rEnd

    rEnd = self % rEnd

  end function getEndPosition

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
  elemental function getParentElementIdx(self) result(parentElementIdx)
    class(coord), intent(in) :: self
    integer(shortInt)        :: parentElementIdx

    parentElementIdx = self % parentElementIdx

  end function getParentElementIdx

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
  pure function getPositionToNudge(self) result(r)
    class(coord), intent(in)    :: self
    real(defReal), dimension(3) :: r

    if (self % nudgeEndPosition) then
      r = self % rEnd

    else
      r = self % r

    end if

  end function getPositionToNudge

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
  elemental function getUniIdx(self) result(uniIdx)
    class(coord), intent(in) :: self
    integer(shortInt)        :: uniIdx

    uniIdx = self % uniIdx

  end function getUniIdx

  !!
  !!
  !!
  elemental function getUniRootId(self) result(uniRootId)
    class(coord), intent(in) :: self
    integer(shortInt)        :: uniRootId

    uniRootId = self % uniRootId

  end function getUniRootId

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

    correct = correct .and. self % uniIdx  > 0
    correct = correct .and. self % localId > 0
    correct = correct .and. self % uniRootId > 0

  end function isValid

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(coord), intent(inout) :: self

    self % r = ZERO
    self % rEnd = ZERO
    self % dir = ZERO
    self % isLeaving = .false.
    self % isRotated = .false.
    self % nudgeEndPosition = .false.
    self % rotMat = ZERO
    self % uniIdx = 0
    self % uniRootId = 0
    self % localId = 0
    self % cellIdx = 0
    self % elementIdx = 0
    self % parentElementIdx = 0

  end subroutine kill

  !!
  !!
  !!
  pure subroutine nudgePosition(self, u)
    class(coord), intent(inout)                       :: self
    real(defReal), dimension(3), intent(in), optional :: u
    real(defReal), dimension(3)                       :: dir

    if (present(u)) then 
      dir = u

    else
      dir = self % dir

    end if

    if (self % nudgeEndPosition) then
      self % rEnd = self % rEnd + dir * NUDGE

    else
      self % r = self % r + dir * NUDGE

    end if

  end subroutine nudgePosition

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
  pure subroutine setEndPosition(self, rEnd)
    class(coord), intent(inout)             :: self
    real(defReal), dimension(3), intent(in) :: rEnd

    self % rEnd = rEnd

  end subroutine setEndPosition

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
  elemental subroutine setNudgeEndPosition(self, nudgeEndPosition)
    class(coord), intent(inout)  :: self
    logical(defBool), intent(in) :: nudgeEndPosition

    self % nudgeEndPosition = nudgeEndPosition

  end subroutine setNudgeEndPosition

  !!
  !!
  !!
  elemental subroutine setParentElementIdx(self, parentElementIdx)
    class(coord), intent(inout)   :: self
    integer(shortInt), intent(in) :: parentElementIdx

    self % parentElementIdx = parentElementIdx

  end subroutine setParentElementIdx

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
  pure subroutine setRotationMatrix(self, rotationMatrix)
    class(coord), intent(inout)                :: self
    real(defReal), dimension(3, 3), intent(in) :: rotationMatrix

    self % rotMat = rotationMatrix

  end subroutine setRotationMatrix

  !!
  !!
  !!
  elemental subroutine setUniIdx(self, uniIdx)
    class(coord), intent(inout)   :: self
    integer(shortInt), intent(in) :: uniIdx

    self % uniIdx = uniIdx

  end subroutine setUniIdx

  !!
  !!
  !!
  elemental subroutine setUniRootId(self, uniRootId)
    class(coord), intent(inout)   :: self
    integer(shortInt), intent(in) :: uniRootId

    self % uniRootId = uniRootId

  end subroutine setUniRootId

end module coord_class