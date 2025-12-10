module coordList_class

  use coord_class,        only : coord
  use genericProcedures,  only : fatalError, numToChar, rotateVector
  use numPrecision
  use publicObjects,      only : coordData, transportObjectStateCoordUpdateData
  use universalVariables, only : HARDCODED_MAX_NEST

  !!
  !! List of co-ordinates at diffrent level of a geometry
  !!
  !! Specifies the position of a particle in space
  !!
  !! It can exist in the following states:
  !!  ABOVE GEOMETRY     -> Nesting = 1. matIdx & uniqueId are < 0 (unassigned). Coordinates at level 1 are reliable.
  !!  PLACED IN GEOMETRY -> Nesting >=1. matIdx is assigned. Coordinates up to nesting are realible.
  !!  UNINITIALISED      -> Is neither PLACED nor ABOVE.
  !!
  !!  NOTE:
  !!   moveGlobal resets regionId & matIdx to 0
  !!   moveLocal  leaves regionId & matIdx unchanged
  !!
  !! Public Members:
  !!   nesting  -> Maximum currently occupied nexting level
  !!   lvl      -> Coordinates at each level
  !!   matIdx   -> Material index at the current position
  !!   uniqueId -> Unique cell Id at the current position
  !!
  !! Interface:
  !!   init              -> Initialise and place ABOVE GEOMETRY, given position and normalised direction
  !!   kill              -> Returns to uninitialised state
  !!   isPlaced          -> True if co-ordinates are PLACED IN GEOMETRY
  !!   isAbove           -> True of co-ordinates are ABOVE GEOMETRY
  !!   isUninitialised   -> True if co-ordinates are UNINITIALISED
  !!   addLevel          -> Increment number of occupied levels by 1.
  !!   decreaseLevel     -> Decrease nesting to a lower level n.
  !!   takeAboveGeometry -> Change state to ABOVE GEOMETRY
  !!   moveGlobal        -> Move point along direction ABOVE GEOMETRY
  !!   moveLocal         -> Move point along direction above and including level n
  !!   rotate            -> Rotate by the cosine of polar deflection mu, and azimuthal angle phi
  !!   cell              -> Return cellIdx at the lowest level
  !!   assignPosition    -> Set position and take ABOVE GEOMETRY
  !!   assignDirection   -> Set direction and do not change coordList state
  !!
  type, public :: coordList
    private
    integer(shortInt)                          :: geometryIdx = 0, materialIdx = -3, nesting = 0, uniqueId = -3
    type(coord), dimension(HARDCODED_MAX_NEST) :: lvl
  contains
    ! Build procedures
    generic            :: init => init_fromPositionAndDirection, init_fromTransportObjectState
    procedure, private :: init_fromPositionAndDirection
    procedure, private :: init_fromTransportObjectState
    procedure          :: kill
    ! State enquiry procedures
    procedure          :: isPlaced
    procedure          :: isAbove
    procedure          :: isUninitialised
    procedure          :: isValid
    ! Interface procedures
    procedure          :: addLevel
    procedure          :: assignDirection
    procedure          :: assignPosition
    procedure          :: decreaseLevel
    procedure          :: getCellIdx
    procedure          :: getCoordinatesData
    procedure          :: getDirection
    procedure          :: getElementIdx
    procedure          :: getGeometryIdx
    procedure          :: getLocalId
    procedure          :: getLowestCellIdx
    procedure          :: getLowestElementIdx
    procedure          :: getLowestMeshIdx
    procedure          :: getMaterialIdx
    procedure          :: getMeshIdx
    procedure          :: getNesting
    procedure          :: getPosition
    procedure          :: getTransportObjectStateUpdateData
    procedure          :: getUniqueId
    procedure          :: getUniverseIdx
    procedure          :: getUniverseRootId
    procedure          :: moveGlobal
    procedure          :: moveLocal
    procedure          :: rotate
    procedure          :: setCellIdx
    procedure          :: setCoordinates
    procedure          :: setDirection
    procedure          :: setElementIdx
    procedure          :: setGeometryIdx
    procedure          :: setIsRotated
    procedure          :: setLocalId
    procedure          :: setLowestElementIdx
    procedure          :: setLowestMeshIdx
    procedure          :: setMaterialIdx
    procedure          :: setMeshIdx
    procedure          :: setNesting
    procedure          :: setPosition
    procedure          :: setPositionAndDirection
    procedure          :: setRotationMatrix
    procedure          :: setUniqueId
    procedure          :: setUniverseIdx
    procedure          :: setUniverseRootId
    procedure          :: takeAboveGeom
    procedure          :: updateCoordinatesFromData
  end type coordList

contains
  !! Add another level of co-ordinates
  !!
  !! Simply increments nesting counter
  !!
  !! Args:
  !!   None
  !!
  pure subroutine addLevel(self)
    class(coordList), intent(inout) :: self

    self % nesting = self % nesting + 1

  end subroutine addLevel

  !!
  !! Assign new direction
  !!
  !! Does not change the state of co-ordinates
  !!
  !! Args:
  !!   u [in] -> New normalised direction at level 1
  !!
  !! NOTE:
  !!   Does not check if u is normalised!
  !!
  pure subroutine assignDirection(self, u)
    class(coordList), intent(inout)         :: self
    real(defReal), dimension(3), intent(in) :: u
    integer(shortInt)                       :: i
    real(defReal), dimension(3)             :: uNew

    ! Initialise uNew = u then propagate changes to all levels.
    uNew = u
    do i = 1, self % nesting
      if (1 < i .and. self % lvl(i) % getIsRotated()) then
        uNew = matmul(self % lvl(i) % getRotationMatrix(), uNew)

      end if
      call self % lvl(i) % setDirection(uNew)

    end do

  end subroutine assignDirection

  !!
  !! Take co-ordinates ABOVE GEOMETRY and assign new position
  !!
  !! Args:
  !!   r [in] -> New position at level 1
  !!
  pure subroutine assignPosition(self, r)
    class(coordList), intent(inout)         :: self
    real(defReal), dimension(3), intent(in) :: r

    call self % takeAboveGeom()
    call self % lvl(1) % setPosition(r)

  end subroutine assignPosition

  !!
  !! Decrease nesting to level n
  !!
  !! Args:
  !!   n [in] -> New nesting level
  !!
  !! Errors:
  !!   fatalError if n is -ve or larger then current nesting
  !!
  subroutine decreaseLevel(self, n)
    class(coordList), intent(inout) :: self
    integer(shortInt), intent(in)   :: n
    character(*), parameter         :: Here = 'decreaseLevel (coordList_class.f90)'

    if (n < 1 .or. self % nesting < n) &
    call fatalError(Here, 'New nesting: '//numToChar(n)//' is invalid. Current nesting is: '//numToChar(self % nesting)//'.')
    self % nesting = n

  end subroutine decreaseLevel

  !!
  !!
  !!
  elemental function getCellIdx(self, lvl) result(cellIdx)
    class(coordList), intent(in)  :: self
    integer(shortInt), intent(in) :: lvl
    integer(shortInt)             :: cellIdx

    cellIdx = self % lvl(lvl) % getCellIdx()

  end function getCellIdx

  !!
  !!
  !!
  elemental function getCoordinatesData(self, lvl) result(data)
    class(coordList), intent(in)  :: self
    integer(shortInt), intent(in) :: lvl
    type(coordData)               :: data

    data = self % lvl(lvl) % getData()

  end function getCoordinatesData

  !!
  !!
  !!
  pure function getDirection(self, lvl) result(u)
    class(coordList), intent(in)  :: self
    integer(shortInt), intent(in) :: lvl
    real(defReal), dimension(3)   :: u

    u = self % lvl(lvl) % getDirection()

  end function getDirection

  !!
  !!
  !!
  elemental function getElementIdx(self, lvl) result(elementIdx)
    class(coordList), intent(in)  :: self
    integer(shortInt), intent(in) :: lvl
    integer(shortInt)             :: elementIdx

    elementIdx = self % lvl(lvl) % getElementIdx()

  end function getElementIdx

  !!
  !!
  !!
  elemental function getGeometryIdx(self) result(geometryIdx)
    class(coordList), intent(in) :: self
    integer(shortInt)            :: geometryIdx

    geometryIdx = self % geometryIdx

  end function getGeometryIdx

  !!
  !!
  !!
  elemental function getLocalId(self, lvl) result(localId)
    class(coordList), intent(in)  :: self
    integer(shortInt), intent(in) :: lvl
    integer(shortInt)             :: localId

    localId = self % lvl(lvl) % getLocalId()

  end function getLocalId

  !!
  !! Returns the index of the cell occupied at the lowest level
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   cellIdx at the lowest ocupied level
  !!
  elemental function getLowestCellIdx(self) result(cellIdx)
    class(coordList), intent(in) :: self
    integer(shortInt)            :: cellIdx

    cellIdx = self % lvl(max(self % nesting, 1)) % getCellIdx()

  end function getLowestCellIdx

  !!
  !!
  !!
  elemental function getLowestElementIdx(self) result(elementIdx)
    class(coordList), intent(in) :: self
    integer(shortInt)            :: elementIdx

    elementIdx = self % lvl(max(self % nesting, 1)) % getElementIdx()

  end function getLowestElementIdx

  !!
  !!
  !!
  elemental function getLowestMeshIdx(self) result(meshIdx)
    class(coordList), intent(in) :: self
    integer(shortInt)            :: meshIdx

    meshIdx = self % lvl(max(self % nesting, 1)) % getMeshIdx()

  end function getLowestMeshIdx

  !!
  !!
  !!
  elemental function getMaterialIdx(self) result(materialIdx)
    class(coordList), intent(in) :: self
    integer(shortInt)            :: materialIdx

    materialIdx = self % materialIdx

  end function getMaterialIdx

  !!
  !!
  !!
  elemental function getMeshIdx(self, lvl) result(meshIdx)
    class(coordList), intent(in)  :: self
    integer(shortInt), intent(in) :: lvl
    integer(shortInt)             :: meshIdx

    meshIdx = self % lvl(lvl) % getMeshIdx()

  end function getMeshIdx

  !!
  !!
  !!
  elemental function getNesting(self) result(nesting)
    class(coordList), intent(in) :: self
    integer(shortInt)            :: nesting

    nesting = self % nesting

  end function getNesting

  !!
  !!
  !!
  elemental function getTransportObjectStateUpdateData(self) result(data)
    class(coordList), intent(in)              :: self
    integer(shortInt)                         :: level
    type(transportObjectStateCoordUpdateData) :: data

    ! Copy into data.
    level = max(1, self % nesting)
    data % geometryIdx = self % geometryIdx
    data % lowestCellIdx = self % lvl(level) % getCellIdx()
    data % lowestElementIdx = self % lvl(level) % getElementIdx()
    data % lowestMeshIdx = self % lvl(level) % getMeshIdx()
    data % materialIdx = self % materialIdx
    data % uniqueId = self % uniqueId
    data % rGlobal = self % lvl(1) % getPosition()
    data % uGlobal = self % lvl(1) % getDirection()

  end function getTransportObjectStateUpdateData

  !!
  !!
  !!
  pure function getPosition(self, lvl) result(r)
    class(coordList), intent(in)  :: self
    integer(shortInt), intent(in) :: lvl
    real(defReal), dimension(3)   :: r

    r = self % lvl(lvl) % getPosition()

  end function getPosition

  !!
  !!
  !!
  elemental function getUniqueId(self) result(uniqueId)
    class(coordList), intent(in) :: self
    integer(shortInt)            :: uniqueId

    uniqueId = self % uniqueId

  end function getUniqueId

  !!
  !!
  !!
  elemental function getUniverseIdx(self, lvl) result(uniIdx)
    class(coordList), intent(in)  :: self
    integer(shortInt), intent(in) :: lvl
    integer(shortInt)             :: uniIdx

    uniIdx = self % lvl(lvl) % getUniverseIdx()

  end function getUniverseIdx

  !!
  !!
  !!
  elemental function getUniverseRootId(self, lvl) result(uniRootId)
    class(coordList), intent(in)  :: self
    integer(shortInt), intent(in) :: lvl
    integer(shortInt)             :: uniRootId

    uniRootId = self % lvl(lvl) % getUniverseRootId()

  end function getUniverseRootId

  !!
  !! Initialise coordList
  !!
  !! Change state from UNINITIALISED to ABOVE GEOMETRY
  !!
  !! Args:
  !!   r [in] -> Position in level 1
  !!   u [in] -> Normalised direction in level 1 (norm2(u) = ONE)
  !!
  !! NOTE:
  !!   Does not check if u is normalised!
  !!
  pure subroutine init_fromPositionAndDirection(self, r, u)
    class(coordList), intent(inout)         :: self
    real(defReal), dimension(3), intent(in) :: r, u

    call self % takeAboveGeom()
    call self % lvl(1) % setPosition(r)
    call self % lvl(1) % setDirection(u)
    self % nesting = 1

  end subroutine init_fromPositionAndDirection

  !!
  !!
  !!
  elemental subroutine init_fromTransportObjectState(self, data)
    class(coordList), intent(inout)                       :: self
    type(transportObjectStateCoordUpdateData), intent(in) :: data

    self % geometryIdx = data % geometryIdx
    self % materialIdx = data % materialIdx
    self % uniqueId = data % uniqueId
    call self % lvl(1) % setPosition(data % rGlobal)
    call self % lvl(1) % setDirection(data % uGlobal)
    self % nesting = 1

  end subroutine init_fromTransportObjectState

  !!
  !! Return true if co-ordinates are above geometry
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   True if co-ordinates are ABOVE GEOMETRY
  !!
  elemental function isAbove(self) result(isIt)
    class(coordList), intent(in) :: self
    logical(defBool)             :: isIt

    isIt = self % materialIdx < 0 .and. self % nesting == 1 .and. self % uniqueId < 0

  end function isAbove

  !!
  !! Return true if co-ordinates List is placed in geometry
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   True if co-ordinates are PLACED.
  !!
  elemental function isPlaced(self) result(isIt)
    class(coordList), intent(in) :: self
    logical(defBool)             :: isIt

    isIt = 0 < self % materialIdx .and. 0 < self % uniqueId .and. 1 <= self % nesting

  end function isPlaced

  !!
  !! Return true if coordinates are uninitialised
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   True if co-ordinates are UNINITIALISED
  !!
  elemental function isUninitialised(self) result(isIt)
    class(coordList), intent(in) :: self
    logical(defBool)             :: isIt

    isIt = .not. (self % isPlaced() .or. self % isAbove())

  end function isUninitialised

  !!
  !!
  !!
  elemental function isValid(self, lvl) result(isIt)
    class(coordList), intent(in)  :: self
    integer(shortInt), intent(in) :: lvl
    logical(defBool)              :: isIt

    isIt = self % lvl(lvl) % isValid()

  end function isValid

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(coordList), intent(inout) :: self

    self % geometryIdx = 0
    self % materialIdx = -3
    self % nesting = 0
    self % uniqueId = -3

    ! Kill coordinates
    call self % lvl % kill()

  end subroutine kill

  !!
  !! Move a point ABOVE the geometry
  !!
  !! Changes state to ABOVE GEOMETRY
  !!
  !! Args:
  !!   d [in] -> Distance (+ve or -ve)
  !!
  !! Errors:
  !!   If d < 0 then movment is backwards.
  !!
  elemental subroutine moveGlobal(self, d)
    class(coordList), intent(inout) :: self
    real(defReal), intent(in)       :: d

    call self % takeAboveGeom()
    call self % lvl(1) % setPosition(self % lvl(1) % getPosition() + d * self % lvl(1) % getDirection())

  end subroutine moveGlobal

  !!
  !! Move point inside the geometry
  !!
  !! Moves above and including level n
  !! Does not change matIdx nor uniqueId
  !!
  !! Args:
  !!   d [in] -> Distance (+ve or -ve)
  !!   n [in] -> Nesting level
  !!
  !! Errors:
  !!   If d < 0.0 movement is backwards
  !!
  subroutine moveLocal(self, d, n)
    class(coordList), intent(inout) :: self
    real(defReal), intent(in)       :: d
    integer(shortInt), intent(in)   :: n
    integer(shortInt)               :: i

    call self % decreaseLevel(n)
    do i = 1, n
      call self % lvl(i) % setPosition(self % lvl(i) % getPosition() + d * self % lvl(i) % getDirection())

    end do

  end subroutine moveLocal

  !!
  !! Rotate direction of the point
  !!
  !! Does not change the state of co-ordinates
  !!
  !! Args:
  !!   mu [in]  -> Cosine of polar deflection angle <-1,1>
  !!   phi [in] -> Azimuthal deflection angle <0;2*pi>
  !!
  elemental subroutine rotate(self, mu, phi)
    class(coordList), intent(inout) :: self
    real(defReal), intent(in)       :: mu
    real(defReal), intent(in)       :: phi
    integer(shortInt)               :: i

    ! Rotate directions in all nesting levels
    call self % lvl(1) % setDirection(rotateVector(self % lvl(1) % getDirection(), mu, phi))

    ! Propagate rotation to lower levels
    do i = 2, self % nesting
      if (self % lvl(i) % getIsRotated()) then
        ! Note that rotation must be performed with the matrix
        ! Deflections by mu & phi depend on coordinates
        ! Deflection by the same my & phi may be diffrent at diffrent, rotated levels!
        call self % lvl(i) % setDirection(matmul(self % lvl(i) % getRotationMatrix(), self % lvl(i - 1) % getDirection()))

      else
        call self % lvl(i) % setDirection(self % lvl(i - 1) % getDirection())

      end if
    end do

  end subroutine rotate

  !!
  !!
  !!
  elemental subroutine setCellIdx(self, cellIdx, lvl)
    class(coordList), intent(inout) :: self
    integer(shortInt), intent(in)   :: cellIdx, lvl

    call self % lvl(lvl) % setCellIdx(cellIdx)

  end subroutine setCellIdx

  !!
  !!
  !!
  elemental subroutine setCoordinates(self, coords, lvl)
    class(coordList), intent(inout) :: self
    type(coord), intent(in)         :: coords
    integer(shortInt), intent(in)   :: lvl

    self % lvl(lvl) = coords

  end subroutine setCoordinates

  !!
  !!
  !!
  pure subroutine setDirection(self, u, lvl)
    class(coordList), intent(inout)         :: self
    real(defReal), dimension(3), intent(in) :: u
    integer(shortInt), intent(in)           :: lvl

    call self % lvl(lvl) % setDirection(u)

  end subroutine setDirection

  !!
  !!
  !!
  elemental subroutine setElementIdx(self, elementIdx, lvl)
    class(coordList), intent(inout) :: self
    integer(shortInt), intent(in)   :: elementIdx, lvl

    call self % lvl(lvl) % setElementIdx(elementIdx)

  end subroutine setElementIdx

  !!
  !!
  !!
  elemental subroutine setGeometryIdx(self, geometryIdx)
    class(coordList), intent(inout) :: self
    integer(shortInt), intent(in)   :: geometryIdx

    self % geometryIdx = geometryIdx

  end subroutine setGeometryIdx

  !!
  !!
  !!
  elemental subroutine setIsRotated(self, isRotated, lvl)
    class(coordList), intent(inout) :: self
    logical(defBool), intent(in)    :: isRotated
    integer(shortInt), intent(in)   :: lvl

    call self % lvl(lvl) % setIsRotated(isRotated)

  end subroutine setIsRotated

  !!
  !!
  !!
  elemental subroutine setLocalId(self, localId, lvl)
    class(coordList), intent(inout) :: self
    integer(shortInt), intent(in)   :: localId, lvl

    call self % lvl(lvl) % setLocalId(localId)

  end subroutine setLocalId

  !!
  !!
  !!
  elemental subroutine setLowestElementIdx(self, lowestElementIdx)
    class(coordList), intent(inout) :: self
    integer(shortInt), intent(in)   :: lowestElementIdx

    call self % lvl(self % nesting) % setElementIdx(lowestElementIdx)

  end subroutine setLowestElementIdx

  !!
  !!
  !!
  elemental subroutine setLowestMeshIdx(self, lowestMeshIdx)
    class(coordList), intent(inout) :: self
    integer(shortInt), intent(in)   :: lowestMeshIdx

    call self % lvl(self % nesting) % setMeshIdx(lowestMeshIdx)

  end subroutine setLowestMeshIdx

  !!
  !!
  !!
  elemental subroutine setMaterialIdx(self, materialIdx)
    class(coordList), intent(inout) :: self
    integer(shortInt), intent(in)   :: materialIdx

    self % materialIdx = materialIdx

  end subroutine setMaterialIdx

  !!
  !!
  !!
  elemental subroutine setMeshIdx(self, meshIdx, lvl)
    class(coordList), intent(inout) :: self
    integer(shortInt), intent(in)   :: meshIdx, lvl

    call self % lvl(lvl) % setMeshIdx(meshIdx)

  end subroutine setMeshIdx

  !!
  !!
  !!
  elemental subroutine setNesting(self, nesting)
    class(coordList), intent(inout) :: self
    integer(shortInt), intent(in)   :: nesting

    self % nesting = nesting

  end subroutine setNesting

  !!
  !!
  !!
  pure subroutine setPosition(self, r, lvl)
    class(coordList), intent(inout)         :: self
    real(defReal), dimension(3), intent(in) :: r
    integer(shortInt), intent(in)           :: lvl

    call self % lvl(lvl) % setPosition(r)

  end subroutine setPosition

  !!
  !!
  !!
  pure subroutine setPositionAndDirection(self, r, u, lvl)
    class(coordList), intent(inout)         :: self
    real(defReal), dimension(3), intent(in) :: r, u
    integer(shortInt), intent(in)           :: lvl

    call self % lvl(lvl) % setPositionAndDirection(r, u)

  end subroutine setPositionAndDirection

  !!
  !!
  !!
  pure subroutine setRotationMatrix(self, rotationMatrix, lvl)
    class(coordList), intent(inout)            :: self
    real(defReal), dimension(3, 3), intent(in) :: rotationMatrix
    integer(shortInt), intent(in)              :: lvl

    call self % lvl(lvl) % setRotationMatrix(rotationMatrix)

  end subroutine setRotationMatrix

  !!
  !!
  !!
  elemental subroutine setUniqueId(self, uniqueId)
    class(coordList), intent(inout) :: self
    integer(shortInt), intent(in)   :: uniqueId

    self % uniqueId = uniqueId

  end subroutine setUniqueId

  !!
  !!
  !!
  elemental subroutine setUniverseIdx(self, universeIdx, lvl)
    class(coordList), intent(inout) :: self
    integer(shortInt), intent(in)   :: universeIdx, lvl

    call self % lvl(lvl) % setUniverseIdx(universeIdx)

  end subroutine setUniverseIdx

  !!
  !!
  !!
  elemental subroutine setUniverseRootId(self, universeRootId, lvl)
    class(coordList), intent(inout) :: self
    integer(shortInt), intent(in)   :: universeRootId, lvl

    call self % lvl(lvl) % setUniverseRootId(universeRootId)

  end subroutine setUniverseRootId

  !!
  !! Takes coordinates above the geometry
  !!
  !! State changes to ABOVE GEOMETRY
  !!
  !! Args:
  !!   None
  !!
  !! NOTE:
  !!   If called on UNINITIALISED may result in unnormalised direction at level 1!
  !!
  elemental subroutine takeAboveGeom(self)
    class(coordList), intent(inout) :: self

    self % materialIdx = -3
    self % nesting = 1
    self % uniqueId = -3

  end subroutine takeAboveGeom

  !!
  !!
  !!
  elemental subroutine updateCoordinatesFromData(self, lvl, data)
    class(coordList), intent(inout) :: self
    integer(shortInt), intent(in)   :: lvl
    type(coordData), intent(in)     :: data

    ! Update coordinates at lowest level then update metadata.
    call self % lvl(lvl) % updateFromData(data)

  end subroutine updateCoordinatesFromData

end module coordList_class