module geometryStd_class

  use charMap_class,     only : charMap
  use coord_class,       only : coord
  use coordList_class,   only : coordList
  use csg_class,         only : csg
  use dictionary_class,  only : dictionary
  use genericProcedures, only : fatalError, numToChar
  use geometry_inter,    only : geometry, distCache
  use materialMenu_mod,  only : nMat
  use numPrecision
  use publicObjects,     only : coordData, newCoordData
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
    procedure          :: init
    procedure          :: kill
    procedure          :: placeCoord
    procedure          :: whatIsAt
    procedure          :: bounds
    procedure          :: move
    procedure          :: moveGlobal
    procedure          :: teleport
    procedure          :: activeMats

    procedure          :: getCellIdx
    ! Private procedures
    procedure, private :: diveToMat
    procedure, private :: closestDist
  end type geometryStd

contains

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
  !! Find material and unique cell at a given location
  !!
  !! See geometry_inter for details
  !!
  subroutine whatIsAt(self, matIdx, uniqueID, r, u)
    class(geometryStd), intent(in)                    :: self
    integer(shortInt), intent(out)                    :: matIdx, uniqueID
    real(defReal), dimension(3), intent(in)           :: r
    real(defReal), dimension(3), optional, intent(in) :: u
    type(coordList)                                   :: coords
    real(defReal), dimension(3)                       :: u_l

    ! If a direction is supplied, update u_l
    u_l = [ONE, ZERO, ZERO]
    if (present(u)) u_l = u

    ! Initialise coordinates
    call coords % init(r, u_l)

    ! Place coordinates
    call self % placeCoord(coords)

    ! Return material & uniqueID
    matIdx = coords % getMatIdx()
    uniqueID = coords % getUniqueId()

  end subroutine whatIsAt

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
    integer(shortInt)                        :: borderIdx, surfIdx, level
    real(defReal)                            :: dist
    type(coordData)                          :: data
    character(*), parameter                  :: Here = 'move (geometryStd_class.f90)'

    if (.not. coords % isPlaced()) call fatalError(Here, 'Coordinate list is not placed in the geometry.')

    ! Find distance to the next surface and reset cache level to 0 afterwards
    call self % closestDist(maxDist, coords, dist, surfIdx, level, cache)
    if (present(cache)) cache % lvl = 0

    if (maxDist < dist) then ! Moves within cell
      ! Move local, register event and return early
      call coords % moveLocal(maxDist, coords % getNesting())
      event = COLL_EV
      return

    end if

    ! Update maxDist if reached this point
    maxDist = dist

    borderIdx = self % geom % getBorderIdx()
    if (surfIdx == borderIdx .and. level == 1) then ! Hits domain boundary
      ! Move global to the boundary and register event
      call coords % moveGlobal(dist)
      event = BOUNDARY_EV

      ! Get boundary surface and apply BCs
      data = newCoordData(coords % getPosition(1), coords % getDirection(1))
      call self % geom % explicitSurfaceBoundaryConditions(borderIdx, data % r, data % u)

      ! Place back in geometry and return early
      call coords % setPositionAndDirection(data % r, data % u, 1)
      call self % placeCoord(coords)
      return

    end if

    ! If reached here then particle crosses to a different local cell. Move to boundary at hit level
    call coords % moveLocal(dist, level)

    ! Register event and update cache level and distance
    event = CROSS_EV
    if (present(cache)) then
      cache % lvl = level - 1
      cache % dist(1:cache % lvl) = cache % dist(1:cache % lvl) - dist

    end if

    ! Get universe and cross to the next cell
    data = coords % getCoordinatesData(level)
    data % surfaceIdx = surfIdx
    call self % geom % crossUniverse(coords % getUniIdx(level), data)

    ! Get material
    call coords % updateCoordinatesFromData(level, data)
    call self % diveToMat(coords, level)

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
    if (coords % getMatIdx() == OUTSIDE_MAT) then
      data = newCoordData(coords % getPosition(1), coords % getDirection(1))
      call self % geom % transformSurfaceBoundaryConditions(self % geom % getBorderIdx(), data % r, data % u)

      ! Return particle to geometry.
      call coords % setPositionAndDirection(data % r, data % u, 1)
      call self % placeCoord(coords)

    end if

  end subroutine teleport

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
  !!
  !!
  function getCellIdx(self, cellId) result(cellIdx)
    class(geometryStd), intent(in) :: self
    integer(shortInt), intent(in)  :: cellId
    integer(shortInt)              :: cellIdx

    cellIdx = self % geom % getCellIdx(cellId)

  end function getCellIdx

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
      call self % geom % getFill(coords % getUniRootId(i), coords % getLocalId(i), fill, uniqueId)

      if (0 <= fill) then ! Found material cell
        call coords % setMatIdx(fill)
        call coords % setUniqueId(uniqueId)
        return

      end if

      ! If reached here we have a universe fill and we descend a level
      if (i == HARDCODED_MAX_NEST) exit ! If there is nested universe at the lowest level

      ! Get current universe
      data = newCoordData(coords % getPosition(i) - &
                          self % geom % getUniverseCellOffset(coords % getUniIdx(i), coords % getLocalId(i)), &
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
  subroutine closestDist(self, maxDist, coords, shortestDist, surfIdx, lvl, cache)
    class(geometryStd), intent(in)           :: self
    real(defReal), intent(in)                :: maxDist
    type(coordList), intent(inout)           :: coords
    real(defReal), intent(out)               :: shortestDist
    integer(shortInt), intent(out)           :: surfIdx, lvl
    type(distCache), intent(inout), optional :: cache
    integer(shortInt)                        :: l, testIdx
    logical(defBool)                         :: update
    real(defReal)                            :: testDistance
    type(coordData)                          :: data

    ! Initialise variables and loop over all geometry levels.
    shortestDist = INF
    surfIdx = 0
    lvl = 0

    do l = 1, coords % getNesting()
      ! Check if cache is present and valid.
      update = .true.
      if (present(cache)) then
        if (l <= cache % lvl) update = .false.

      end if

      if (update) then
        ! Get universe and compute distance.
        data = coords % getCoordinatesData(l)
        data % dMax = maxDist
        call self % geom % distanceUniverse(data)
        call coords % updateCoordinatesFromData(l, data)
        testDistance = data % d
        testIdx = data % surfaceIdx

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
      if ((shortestDist - testDistance) < shortestDist * FP_REL_TOL) cycle
      shortestDist = testDistance
      surfIdx = testIdx
      lvl = l

    end do

  end subroutine closestDist

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

end module geometryStd_class