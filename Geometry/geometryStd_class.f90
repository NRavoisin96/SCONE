module geometryStd_class

  use charMap_class,     only : charMap
  use coord_class,       only : coord
  use coordList_class,   only : coordList
  use csg_class,         only : csg
  use dictionary_class,  only : dictionary
  use genericProcedures, only : fatalError, numToChar
  use geometry_inter,    only : geometry, distCache
  use numPrecision
  use surface_inter,     only : surface
  use universalVariables
  use universe_inter,    only : universe

  ! Nuclear Data
  use materialMenu_mod,  only : nMat

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
    type(csg) :: geom
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: placeCoord
    procedure :: whatIsAt
    procedure :: bounds
    procedure :: move_noCache
    procedure :: move_withCache
    procedure :: moveGlobal
    procedure :: teleport
    procedure :: activeMats
    ! Private procedures
    procedure, private :: diveToMat
    procedure, private :: closestDist
    procedure, private :: closestDist_cache
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
  elemental subroutine kill(self)
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
    class(universe), pointer       :: uni
    type(coord)                    :: new
    character(100), parameter      :: Here = 'placeCoord (geometryStd_class.f90)'

    ! Check that coordList is initialised.
    nesting = coords % getNesting()
    if (nesting < 1) call fatalError(Here, 'CoordList is not initialised. Nesting is: '//numToChar(nesting)//'.')

    ! Place coordinates above geometry (in case they were placed)
    call coords % takeAboveGeom()

    ! Enter root universe.
    uni => self % geom % unis % getPtr_fast(self % geom % rootIdx)
    call uni % enter(coords % getPosition(1), coords % getDirection(1), new)

    ! Set new coordinates in the list.
    call new % setUniRootId(1)
    call coords % setCoordinates(new, 1)

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
    real(defReal), dimension(3)                       :: u_l = [ONE, ZERO, ZERO]

    ! If a direction is supplied, update u_l
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
    class(surface), pointer        :: surf
    integer(shortInt)              :: i

    ! Get boundary surface
    surf => self % geom % surfs % getPtr(self % geom % borderIdx)
    bounds = surf % getBoundingBox()

    ! Change infinite dimensions to ZERO
    do i = 1, 3
      if (bounds(i) <= -INF .and. bounds(i + 3) >= INF) then
        bounds(i) = ZERO
        bounds(i + 3) = ZERO
      end if
    end do

  end function bounds

  !!
  !! Given coordinates placed in the geometry move point through the geometry
  !!
  !! See geometry_inter for details
  !!
  !! Uses explicit BC
  !!
  subroutine move_noCache(self, coords, maxDist, event)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    integer(shortInt)              :: surfIdx, level
    real(defReal)                  :: dist
    class(surface), pointer        :: surf
    real(defReal), dimension(3)    :: r, u
    class(universe), pointer       :: uni
    type(coord)                    :: levelCoords
    character(100), parameter      :: Here = 'move (geometryStd_class.f90)'

    if (.not. coords % isPlaced()) call fatalError(Here, 'Coordinate list is not placed in the geometry.')

    ! Find distance to the next surface
    call self % closestDist(coords, maxDist, dist, surfIdx, level)

    if (maxDist < dist) then ! Moves within cell
      ! Move local, register event and return early
      call coords % moveLocal(maxDist, coords % getNesting())
      event = COLL_EV
      return

    end if

    ! Update maxDist if reached this point
    maxDist = dist

    if (surfIdx == self % geom % borderIdx .and. level == 1) then ! Hits domain boundary
      ! Move global to the boundary
      call coords % moveGlobal(dist)
      event = BOUNDARY_EV

      ! Get boundary surface and apply boundary conditions
      surf => self % geom % surfs % getPtr(self % geom % borderIdx)
      r = coords % getPosition(1)
      u = coords % getDirection(1)
      call surf % explicitBC(r, u)

      ! Place back in geometry and return early
      call coords % setPosition(r, 1)
      call coords % setDirection(u, 1)
      call self % placeCoord(coords)
      return

    end if

    ! If reached here then particle crosses to a different local cell. Move to boundary at hit level
    call coords % moveLocal(dist, level)
    event = CROSS_EV

    ! Get universe and cross to the next cell
    uni => self % geom % unis % getPtr_fast(coords % getUniIdx(level))
    levelCoords = coords % getCoordinates(level)
    call uni % cross(levelCoords, surfIdx)

    ! Get material
    call coords % setCoordinates(levelCoords, level)
    call self % diveToMat(coords, level)

  end subroutine move_noCache

  !!
  !! Given coordinates placed in the geometry move point through the geometry
  !!
  !! See geometry_inter for details
  !!
  !! Uses explicit BC
  !!
  subroutine move_withCache(self, coords, maxDist, event, cache)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(inout)   :: maxDist
    integer(shortInt), intent(out) :: event
    type(distCache), intent(inout) :: cache
    integer(shortInt)              :: surfIdx, level
    real(defReal)                  :: dist
    class(surface), pointer        :: surf
    real(defReal), dimension(3)    :: r, u
    class(universe), pointer       :: uni
    type(coord)                    :: levelCoords
    character(100), parameter      :: Here = 'move_withCache (geometryStd_class.f90)'

    if (.not. coords % isPlaced()) call fatalError(Here, 'Coordinate list is not placed in the geometry.')

    ! Find distance to the next surface and reset cache level to 0 afterwards
    call self % closestDist_cache(coords, cache, maxDist, dist, surfIdx, level)
    cache % lvl = 0

    if (maxDist < dist) then ! Moves within cell
      ! Move local, register event and return early
      call coords % moveLocal(maxDist, coords % getNesting())
      event = COLL_EV
      return

    end if

    ! Update maxDist if reached this point
    maxDist = dist

    if (surfIdx == self % geom % borderIdx .and. level == 1) then ! Hits domain boundary
      ! Move global to the boundary and register event
      call coords % moveGlobal(dist)
      event = BOUNDARY_EV

      ! Get boundary surface and apply BCs
      surf => self % geom % surfs % getPtr(self % geom % borderIdx)
      r = coords % getPosition(1)
      u = coords % getDirection(1)
      call surf % explicitBC(r, u)

      ! Place back in geometry and return early
      call coords % setPosition(r, 1)
      call coords % setDirection(u, 1)
      call self % placeCoord(coords)
      return

    end if

    ! If reached here then particle crosses to a different local cell. Move to boundary at hit level
    call coords % moveLocal(dist, level)

    ! Register event and update cache level and distance
    event = CROSS_EV
    cache % lvl = level - 1
    cache % dist(1:level - 1) = cache % dist(1:level - 1) - dist

    ! Get universe and cross to the next cell
    uni => self % geom % unis % getPtr_fast(coords % getUniIdx(level))
    levelCoords = coords % getCoordinates(level)
    call uni % cross(levelCoords, surfIdx)

    ! Get material
    call coords % setCoordinates(levelCoords, level)
    call self % diveToMat(coords, level)

  end subroutine move_withCache

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
    class(surface), pointer        :: surf
    real(defReal), dimension(3)    :: r, u
    real(defReal)                  :: dist

    ! Initialise event = COLL_EV and get boundary surface.
    event = COLL_EV
    surf => self % geom % surfs % getPtr(self % geom % borderIdx)

    ! Find distance to the boundary
    dist = surf % distance(coords % getPosition(1), coords % getDirection(1))

    ! Check if dist < maxDist. If so, update maxDist and event
    if (dist < maxDist) then
      maxDist = dist
      event = BOUNDARY_EV

    end if

    ! Move global and apply boundary conditions if applicable.
    call coords % moveGlobal(maxDist)
    if (event == BOUNDARY_EV) then
      r = coords % getPosition(1)
      u = coords % getDirection(1)
      call surf % explicitBC(r, u)
      call coords % setPosition(r, 1)
      call coords % setDirection(u, 1)

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
    class(surface), pointer        :: surf
    real(defReal), dimension(3)    :: r, u

    ! Move the coords above the geometry
    call coords % moveGlobal(dist)

    ! Place coordinates back into geometry
    call self % placeCoord(coords)

    ! If point is outside apply boundary transformations
    if (coords % getMatIdx() == OUTSIDE_MAT) then
      surf => self % geom % surfs % getPtr(self % geom % borderIdx)
      r = coords % getPosition(1)
      u = coords % getDirection(1)
      call surf % transformBC(r, u)

      ! Return particle to geometry.
      call coords % setPosition(r, 1)
      call coords % setDirection(u, 1)
      call self % placeCoord(coords)
    end if

  end subroutine teleport

  !!
  !! Returns the list of active materials used in the geometry
  !!
  !! See geometry_inter for details
  !!
  !! NOTE: This function uses VOID_MAT and UNDEF_MAT from universalVariables
  !!
  function activeMats(self) result(matList)
    class(geometryStd), intent(in)               :: self
    integer(shortInt), dimension(:), allocatable :: matList
    integer(shortInt)                            :: N, lastIdx

    ! Takes the list of materials present in the geometry from geomGraph
    N = size(self % geom % graph % usedMats)
    lastIdx = self % geom % graph % usedMats(N)

    ! Check if the last entry of the list is an actual material or void
    if (lastIdx == VOID_MAT) then
      N = N - 1
      lastIdx = self % geom % graph % usedMats(N)
    end if
  
    ! Check if the last entry of the list is an undefined material and if so update N
    if (lastIdx == UNDEF_MAT) N = N - 1
    matList = self % geom % graph % usedMats(1:N)

  end function activeMats

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
    class(universe), pointer       :: uni
    type(coord)                    :: new
    real(defReal), dimension(3)    :: offset
    character(100), parameter      :: Here = 'diveToMat (geometryStd_class.f90)'

    do i = start, HARDCODED_MAX_NEST
      ! Find cell fill
      call self % geom % graph % getFill(coords % getUniRootId(i), coords % getLocalId(i), fill, uniqueId)

      if (fill >= 0) then ! Found material cell
        call coords % setMatIdx(fill)
        call coords % setUniqueId(uniqueId)
        return

      end if

      ! If reached here we have a universe fill and we descend a level
      if (i == HARDCODED_MAX_NEST) exit ! If there is nested universe at the lowest level
      fill = abs(fill)

      ! Get current universe
      uni => self % geom % unis % getPtr_fast(coords % getUniIdx(i))

      ! Get cell offset
      offset = uni % cellOffset(coords % getCoordinates(i))

      ! Get nested universe
      uni => self % geom % unis % getPtr_fast(fill)

      ! Enter nested universe
      call coords % addLevel()
      call uni % enter(coords % getPosition(i) - offset, coords % getDirection(i), new)

      ! Set new % uniRootId and place into coordList.
      call new % setUniRootId(uniqueId)
      call coords % setCoordinates(new, i + 1)

    end do

    call fatalError(Here, 'Failed to find material cell. Should not happen after &
                    &geometry checks during build...')

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
  subroutine closestDist(self, coords, maxDist, dist, surfIdx, lvl)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    real(defReal), intent(in)      :: maxDist
    real(defReal), intent(out)     :: dist
    integer(shortInt), intent(out) :: surfIdx, lvl
    integer(shortInt)              :: l, testIdx
    real(defReal)                  :: testDist
    class(universe), pointer       :: uni
    type(coord)                    :: levelCoords

    ! Initialise variables and loop over all geometry levels.
    dist = INF
    surfIdx = 0
    lvl = 0
    do l = 1, coords % getNesting()
      ! Get universe and compute distance.
      levelCoords = coords % getCoordinates(l)
      call levelCoords % setEndPosition(levelCoords % getPosition() + maxDist * levelCoords % getDirection())
      uni => self % geom % unis % getPtr_fast(levelCoords % getUniIdx())
      call uni % distance(levelCoords, testDist, testIdx)
      call coords % setCoordinates(levelCoords, l)

      ! Save distance, surfIdx & level coresponding to shortest distance
      ! Take FP precision into account
      if ((dist - testDist) < dist * FP_REL_TOL) cycle
      dist = testDist
      surfIdx = testIdx
      lvl = l

    end do

  end subroutine closestDist

  !!
  !! Return distance to the closest surface
  !!
  !! Searches through all geometry levels. In addition to distance return level
  !! and surfIdx for crossing surface
  !!
  !! Args:
  !!   coords [inout] -> Current coordinates of a particle
  !!   cache [inout]  -> Distance cache. Use valid distances from cache. Put calculated
  !!                     distances on the cache.
  !!   maxDist        -> Maximum distance of travel.
  !!   dist [out]     -> Value of closest distance
  !!   surfIdx [out]  -> Surface index for the crossing returned from the universe
  !!   lvl [out]      -> Level at which crossing is closest
  !!
  subroutine closestDist_cache(self, coords, cache, maxDist, dist, surfIdx, lvl)
    class(geometryStd), intent(in) :: self
    type(coordList), intent(inout) :: coords
    type(distCache), intent(inout) :: cache
    real(defReal), intent(in)      :: maxDist
    real(defReal), intent(out)     :: dist
    integer(shortInt), intent(out) :: surfIdx, lvl
    integer(shortInt)              :: l, testIdx
    real(defReal)                  :: testDist
    class(universe), pointer       :: uni
    type(coord)                    :: levelCoords

    ! Initialise variables and loop over all geometry levels.
    dist = INF
    surfIdx = 0
    lvl = 0
    do l = 1, coords % getNesting()
      ! Update Cache if distance is not valid
      if (cache % lvl < l) then
        ! Get universe
        levelCoords = coords % getCoordinates(l)
        uni => self % geom % unis % getPtr_fast(levelCoords % getUniIdx())

        ! Find distance
        call levelCoords % setEndPosition(levelCoords % getPosition() + maxDist * levelCoords % getDirection())
        call uni % distance(levelCoords, cache % dist(l), cache % surf(l))
        call coords % setCoordinates(levelCoords, l)
        cache % lvl = cache % lvl + 1

      end if

      ! Read distance and crossing memento from cache
      testDist = cache % dist(l)
      testIdx  = cache % surf(l)

      ! Save distance, surfIdx & level coresponding to shortest distance
      ! Take FP precision into account
      if ((dist - testDist) < dist * FP_REL_TOL) cycle
      dist = testDist
      surfIdx = testIdx
      lvl = l

    end do

  end subroutine closestDist_cache

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