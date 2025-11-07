module geometryStd_class

  use charMap_class,     only : charMap
  use coord_class,       only : coord
  use coordList_class,   only : coordList
  use csg_class,         only : csg
  use dictionary_class,  only : dictionary
  use genericProcedures, only : fatalError, numToChar
  use geometry_inter,    only : geometry, distCache
  use materialMenu_mod,  only : nMat
  use mesh_inter,        only : mesh
  use numPrecision
  use publicObjects,     only : coordData, newCoordData
  use RNG_class,         only : RNG
  use scalarField_inter, only : getTemperatureFieldPtr, scalarField
  use universalVariables

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: geometryStd_CptrCast

  !!
  !! Standard Geometry Model
  !!
  !! Typical geometry of a MC Neutron Transport code composed of multiple nested
  !! universes.
  !!
  !! Boundary conditions in diffrent movement models are handeled:
  !!   move       -> explicitBC
  !!   moveGlobal -> explicitBC
  !!   teleport   -> Co-ordinate transfrom
  !!
  !! Sample Dictionary Input:
  !!   geometry {
  !!     type geometryStd;
  !!     <csg_class difinition>
  !!    }
  !!
  !! Public Members:
  !!   geom -> Representation of geometry by csg_class. Contains all surfaces, cells and universe
  !!     as well as geometry graph and info about root uni and boundary surface.
  !!
  !! Interface:
  !!   Geometry Interface
  !!
  type, public, extends(geometry) :: geometryStd
    private
    type(csg)                     :: geom
  contains
    ! Superclass procedures
    procedure          :: activeMats
    procedure          :: bounds
    procedure, private :: closestDist
    procedure, private :: diveToMat
    procedure          :: getCellIdx
    procedure          :: getMeshIdxByName
    procedure          :: getMeshPtr
    procedure          :: init
    procedure          :: kill
    procedure          :: move
    procedure          :: moveGlobal
    procedure          :: placeCoord
    procedure          :: sampleInitialPosition
    procedure          :: teleport
    procedure          :: whatIsAt
  end type geometryStd

contains
  !!
  !! Returns the list of active materials used in the geometry
  !!
  !! See geometry_inter for details
  !!
  pure function activeMats(self) result(matList)
    class(geometryStd), intent(in)               :: self
    integer(shortInt), dimension(:), allocatable :: matList

    matList = self % geom % getActiveMaterialIdxs()

  end function activeMats

  !!
  !! Return Axis Aligned Bounding Box encompassing the geometry
  !!
  !! See geometry_inter for details
  !!
  function bounds(self)
    class(geometryStd), intent(in) :: self
    real(defReal), dimension(6)    :: bounds
    integer(shortInt)              :: i

    ! Get boundary surface
    bounds = self % geom % getSurfaceBounds(self % geom % getBorderIdx())

    ! Change infinite dimensions to ZERO
    do i = 1, 3
      if (bounds(i) <= -INF .and. bounds(i + 3) >= INF) bounds([i, i + 3]) = ZERO

    end do

  end function bounds

  !!
  !! Return distance to the closest surface
  !!
  !! Searches through all geometry levels. In addition to distance return level
  !! and surfIdx for crossing surface
  !!
  !! Args:
  !!   coords [inout] -> Current coordinates of a particle
  !!   maxDist [in]   -> Maximum distance of travel
  !!   dist [out]     -> Value of closest distance
  !!   surfIdx [out]  -> Surface index for the crossing returned from the universe
  !!   lvl     [out]  -> Level at which crossing is closest
  !!
  subroutine closestDist(self, maxDist, coords, updateData, cache)
    class(geometryStd), intent(in)           :: self
    real(defReal), intent(in)                :: maxDist
    type(coordList), intent(in)              :: coords
    type(coordData), intent(out)             :: updateData
    type(distCache), intent(inout), optional :: cache
    integer(shortInt)                        :: l, testIdx
    logical(defBool)                         :: update
    real(defReal)                            :: testDistance
    type(coordData)                          :: levelData

    ! Loop over all geometry levels.
    do l = 1, coords % getNesting()
      ! Check if cache is present and valid.
      update = .true.
      if (present(cache)) then
        if (l <= cache % lvl) update = .false.

      end if

      if (update) then
        ! Get universe and compute distance.
        levelData = coords % getCoordinatesData(l)
        levelData % dMax = min(maxDist, testDistance)
        call self % geom % distanceUniverse(levelData)
        testDistance = levelData % d
        testIdx = levelData % surfaceIdx

        if (present(cache)) then
          ! Update cache and mark this level as valid.
          cache % dist(l) = testDistance
          cache % surf(l) = testIdx
          cache % lvl = l

        end if

      else
        testDistance = cache % dist(l)
        testIdx = cache % surf(l)

      end if

      ! Save distance, surfIdx & level coresponding to shortest distance
      ! Take FP precision into account
      if (updateData % d * FP_REL_TOL <= updateData % d - testDistance) then
        updateData = levelData
        updateData % surfaceIdx = testIdx
        updateData % updateLevel = l

      end if

    end do

  end subroutine closestDist

  !!
  !! Descend down the geometry structure until material is reached
  !!
  !! Requires starting level to be specified.
  !! It is a private procedure common to all movement types in geometry.
  !!
  !! Args:
  !!   coords [inout] -> CoordList of a particle. Assume that coords are already valid for all
  !!     levels above and including start
  !!   start [in] -> Starting level for material search
  !!
  !! Errors:
  !!   fatalError if material cell is not found until maximum nesting is reached
  !!
  subroutine diveToMat(self, coords, start)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    integer(shortInt), intent(in)  :: start
    integer(shortInt)              :: fill, uniqueId, i
    type(coordData)                :: data
    character(*), parameter        :: Here = 'diveToMat (geometryStd_class.f90)'

    do i = start, HARDCODED_MAX_NEST
      ! Find cell fill
      call self % geom % getFill(coords % getUniverseRootId(i), coords % getLocalId(i), fill, uniqueId)

      if (0 <= fill) then ! Found material cell
        call coords % setMaterialIdx(fill)
        call coords % setUniqueId(uniqueId)
        return

      end if

      ! If reached here we have a universe fill and we descend a level
      if (i == HARDCODED_MAX_NEST) exit ! If there is nested universe at the lowest level

      ! Get current universe
      data = newCoordData(coords % getPosition(i) - &
                          self % geom % getUniverseCellOffset(coords % getUniverseIdx(i), coords % getLocalId(i)), &
                          coords % getDirection(i), universeRootId = uniqueId)

      ! Enter nested universe
      call self % geom % enterUniverse(abs(fill), data)

      ! Set new % uniRootId and place into coordList.
      call coords % addLevel()
      call coords % updateCoordinatesFromData(i + 1, data)

    end do

    call fatalError(Here, 'Failed to find material cell.')

  end subroutine diveToMat

  !!
  !! Cast geometry pointer to geometryStd class pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class geometry
  !!
  !! Result:
  !!   Null if source is not of geometryStd class
  !!   Target points to source if source is geometryStd class
  !!
  pure function geometryStd_CptrCast(source) result(ptr)
    class(geometry), pointer, intent(in) :: source
    class(geometryStd), pointer          :: ptr

    select type(source)
      class is (geometryStd)
        ptr => source

      class default
        ptr => null()

    end select

  end function geometryStd_CptrCast

  !!
  !!
  !!
  function getCellIdx(self, cellId) result(cellIdx)
    class(geometryStd), intent(in) :: self
    integer(shortInt), intent(in)  :: cellId
    integer(shortInt)              :: cellIdx

    cellIdx = self % geom % getCellIdx(cellId)

  end function getCellIdx

  !!
  !!
  !!
  function getMeshIdxByName(self, name) result(idx)
    class(geometryStd), intent(in) :: self
    character(nameLen), intent(in) :: name
    integer(shortInt)              :: idx

    idx = self % geom % getMeshIdxByName(name)

  end function getMeshIdxByName

  !!
  !!
  !!
  function getMeshPtr(self, id) result(meshPtr)
    class(geometryStd), intent(in) :: self
    integer(shortInt), intent(in)  :: id
    class(mesh), pointer           :: meshPtr

    meshPtr => self % geom % getMeshPtr(id)

  end function getMeshPtr

  !!
  !! Initialise geometry
  !!
  !! See geometry_inter for details
  !!
  subroutine init(self, dict, mats, silent)
    class(geometryStd), intent(inout)      :: self
    class(dictionary), intent(in)          :: dict
    type(charMap), intent(in)              :: mats
    logical(defBool), optional, intent(in) :: silent

    ! Build the representation
    call self % geom % init(dict, mats, silent)

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(geometryStd), intent(inout) :: self

    call self % geom % kill()

  end subroutine kill

!!
  !! Given coordinates placed in the geometry move point through the geometry
  !!
  !! See geometry_inter for details
  !!
  !! Uses explicit BC
  !!
  subroutine move(self, coords, maxDist, event, cache)
    class(geometryStd), intent(in)           :: self
    type(coordList), intent(inout)           :: coords
    real(defReal), intent(inout)             :: maxDist
    integer(shortInt), intent(out)           :: event
    type(distCache), intent(inout), optional :: cache
    integer(shortInt)                        :: borderIdx
    type(coordData)                          :: updateData
    character(*), parameter                  :: Here = 'move (geometryStd_class.f90)'

    if (.not. coords % isPlaced()) call fatalError(Here, 'Coordinate list is not placed in the geometry.')

    ! Find distance to the next surface and reset cache level to 0 afterwards
    call self % closestDist(maxDist, coords, updateData, cache)
    if (present(cache)) cache % lvl = 0

    if (maxDist < updateData % d) then ! Moves within cell
      ! Move local, register event and return early
      call coords % moveLocal(maxDist, coords % getNesting())
      event = COLL_EV
      return

    end if

    ! If reached this point a boundary was hit. Update maxDist.
    maxDist = updateData % d

    borderIdx = self % geom % getBorderIdx()
    if (updateData % surfaceIdx == borderIdx .and. updateData % updateLevel == 1) then ! Hits domain boundary
      ! Move global to the boundary and register event
      call coords % moveGlobal(updateData % d)
      event = BOUNDARY_EV

      ! Get boundary surface, apply BCs and place back in geometry.
      updateData = newCoordData(coords % getPosition(1), coords % getDirection(1))
      call self % geom % explicitSurfaceBoundaryConditions(borderIdx, updateData % r, updateData % u)
      call coords % setPositionAndDirection(updateData % r, updateData % u, 1)
      call self % placeCoord(coords)

    else
      ! The particle moves locally within its cell. Move to boundary at hit level.
      call coords % moveLocal(updateData % d, updateData % updateLevel)
      updateData % r = coords % getPosition(updateData % updateLevel)

      ! Register event and update cache level and distance
      event = CROSS_EV
      if (present(cache)) then
        cache % lvl = updateData % updateLevel - 1
        cache % dist(1:cache % lvl) = cache % dist(1:cache % lvl) - updateData % d

      end if

      ! Get universe and cross to the next cell
      call self % geom % crossUniverse(updateData % universeIdx, updateData)

      ! Get material
      call coords % updateCoordinatesFromData(updateData % updateLevel, updateData)
      call self % diveToMat(coords, updateData % updateLevel)

    end if

  end subroutine move

  !!
  !! Move a particle in the top (global) level in the geometry
  !!
  !! See geometry_inter for details
  !!
  !! Uses explicit BC
  !!
  subroutine moveGlobal(self, coords, maxDist, event)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    integer(shortInt)              :: borderIdx
    type(coordData)                :: data
    real(defReal)                  :: dist

    ! Initialise event = COLL_EV and get boundary surface.
    event = COLL_EV

    ! Find distance to the boundary
    borderIdx = self % geom % getBorderIdx()
    data = newCoordData(coords % getPosition(1), coords % getDirection(1))
    dist = self % geom % distanceSurface(borderIdx, data % r, data % u)

    ! Check if dist < maxDist. If so, update maxDist and event
    if (dist < maxDist) then
      maxDist = dist
      event = BOUNDARY_EV

    end if

    ! Move global and apply boundary conditions if applicable.
    call coords % moveGlobal(maxDist)
    if (event == BOUNDARY_EV) then
      data = newCoordData(coords % getPosition(1), coords % getDirection(1))
      call self % geom % explicitSurfaceBoundaryConditions(borderIdx, data % r, data % u)
      call coords % setPositionAndDirection(data % r, data % u, 1)

    end if

    ! Return particle to geometry
    call self % placeCoord(coords)

  end subroutine moveGlobal

  !!
  !! Place coordinate list into geometry
  !!
  !! See geometry_inter for details
  !!
  subroutine placeCoord(self, coords)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    integer(shortInt)              :: nesting
    type(coordData)                :: data
    character(*), parameter        :: Here = 'placeCoord (geometryStd_class.f90)'

    ! Check that coordList is initialised.
    nesting = coords % getNesting()
    if (nesting < 1) call fatalError(Here, 'CoordList is not initialised. Nesting is: '//numToChar(nesting)//'.')

    ! Place coordinates above geometry (in case they were placed)
    call coords % takeAboveGeom()

    ! Enter root universe.
    data = newCoordData(coords % getPosition(1), coords % getDirection(1), universeRootId = 1)
    call self % geom % enterUniverse(self % geom % getRootIdx(), data)

    ! Set new coordinates in the list.
    call coords % updateCoordinatesFromData(1, data)

    ! Dive to material
    call self % diveToMat(coords, 1)

  end subroutine placeCoord

  !!
  !!
  !!
  subroutine sampleInitialPosition(self, bottom, top, rand, materialIdx, uniqueId, r, temperature)
    class(geometryStd), intent(in)           :: self
    real(defReal), dimension(3), intent(in)  :: bottom, top
    class(RNG), intent(inout)                :: rand
    integer(shortInt), intent(out)           :: materialIdx, uniqueId
    real(defReal), dimension(3), intent(out) :: r
    real(defReal), intent(out), optional     :: temperature
    real(defReal), dimension(3)              :: randomNumbers

    ! Sample position.
    call rand % generate(randomNumbers)
    r = (top - bottom) * randomNumbers + bottom

    ! Find material under position.
    call self % whatIsAt(materialIdx, uniqueId, r, temperature = temperature)

  end subroutine sampleInitialPosition

  !!
  !! Move a particle in the top level without stopping
  !!
  !! See geometry_inter for details
  !!
  !! Uses co-ordinate transform boundary XSs
  !!
  subroutine teleport(self, coords, dist)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(in)      :: dist
    type(coordData)                :: data

    ! Move the coords above the geometry
    call coords % moveGlobal(dist)

    ! Place coordinates back into geometry
    call self % placeCoord(coords)

    ! If point is outside apply boundary transformations
    if (coords % getMaterialIdx() == OUTSIDE_MAT) then
      data = newCoordData(coords % getPosition(1), coords % getDirection(1))
      call self % geom % transformSurfaceBoundaryConditions(self % geom % getBorderIdx(), data % r, data % u)

      ! Return particle to geometry.
      call coords % setPositionAndDirection(data % r, data % u, 1)
      call self % placeCoord(coords)

    end if

  end subroutine teleport

  !!
  !! Find material and unique cell at a given location
  !!
  !! See geometry_inter for details
  !!
  subroutine whatIsAt(self, matIdx, uniqueID, r, u, temperature)
    class(geometryStd), intent(in)                    :: self
    integer(shortInt), intent(out)                    :: matIdx, uniqueID
    real(defReal), dimension(3), intent(in)           :: r
    real(defReal), dimension(3), optional, intent(in) :: u
    real(defReal), intent(out), optional              :: temperature
    class(scalarField), pointer                       :: temperatureFieldPtr
    real(defReal), dimension(3)                       :: u_l
    type(coordList)                                   :: coords

    ! If a direction is supplied, update u_l
    u_l = [ONE, ZERO, ZERO]
    if (present(u)) u_l = u

    ! Initialise coordinates
    call coords % init(r, u_l)

    ! Place coordinates
    call self % placeCoord(coords)

    ! Return material & uniqueID
    matIdx = coords % getMaterialIdx()
    uniqueID = coords % getUniqueId()

    ! Get temperature if requested.
    if (present(temperature)) then
      temperature = ZERO
      temperatureFieldPtr => getTemperatureFieldPtr()
      if (associated(temperatureFieldPtr)) temperature = temperatureFieldPtr % at(coords)

    end if

  end subroutine whatIsAt

end module geometryStd_class