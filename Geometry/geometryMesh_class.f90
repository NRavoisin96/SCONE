module geometryMesh_class

  use charMap_class,     only : charMap
  use coordList_class,   only : coordList
  use dictionary_class,  only : dictionary
  use genericProcedures, only : ceilingBinarySearch, fatalError, numToChar, removeDuplicates
  use geometry_inter,    only : distCache, geometry
  use meshShelf_class,   only : meshBox, meshShelf
  use mesh_inter,        only : mesh
  use numPrecision
  use publicObjects,     only : coordData, newCoordData
  use RNG_class,         only : RNG
  use scalarField_inter, only : getTemperatureFieldPtr, scalarField
  use universalVariables

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(geometry) :: geometryMesh
    private
    integer(shortInt), dimension(:), allocatable    :: elementIdOffsets, fills, localIdOffsets
    integer(shortInt), dimension(:, :), allocatable :: cumulativeToActiveElementsMap
    real(defReal)                                   :: totalVolume = ZERO
    real(defReal), dimension(:), allocatable        :: cumulativeVolumes
    type(meshShelf)                                 :: meshes
  contains
    procedure :: activeMats
    procedure :: bounds
    procedure :: getMeshIdxByName
    procedure :: getMeshPtr
    procedure :: init
    procedure :: kill
    procedure :: move
    procedure :: moveGlobal
    procedure :: placeCoord
    procedure :: sampleInitialPosition
    procedure :: teleport
    procedure :: whatIsAt
  end type geometryMesh

contains
  !!
  !!
  !!
  pure function activeMats(self) result(matList)
    class(geometryMesh), intent(in)              :: self
    integer(shortInt), dimension(:), allocatable :: matList

    matList = removeDuplicates(self % fills)

  end function activeMats

  !!
  !!
  !!
  function bounds(self)
    class(geometryMesh), intent(in) :: self
    real(defReal), dimension(6)     :: bounds

    bounds = reshape(self % meshes % getOverallBoundingBoxBounds(), [6])

  end function bounds

  !!
  !!
  !!
  function getMeshIdxByName(self, name) result(idx)
    class(geometryMesh), intent(in) :: self
    character(nameLen), intent(in)  :: name
    integer(shortInt)               :: idx

    idx = self % meshes % getMeshIdxByName(name)

  end function getMeshIdxByName

  !!
  !!
  !!
  function getMeshPtr(self, id) result(meshPtr)
    class(geometryMesh), intent(in) :: self
    integer(shortInt), intent(in)   :: id
    class(mesh), pointer            :: meshPtr

    meshPtr => self % meshes % getMeshPtr(self % meshes % getMeshIdx(id))

  end function getMeshPtr

  !!
  !!
  !!
  subroutine init(self, dict, mats, silent)
    class(geometryMesh), intent(inout)            :: self
    class(dictionary), intent(in)                 :: dict
    type(charMap), intent(in)                     :: mats
    logical(defBool), intent(in), optional        :: silent
    character(nameLen), dimension(:), allocatable :: meshFills
    class(dictionary), pointer                    :: meshDict, meshesDict
    integer(shortInt)                             :: cumulativeVolumeIdx, i, iLessOne, j, fill, fillIdx, &
                                                     nActiveElements, nFills, nLocalIds, nLocalIdsInMesh, nMeshes
    logical(defBool)                              :: loud
    real(defReal)                                 :: volume
    type(meshBox)                                 :: box
    character(*), parameter                       :: here = 'init (geometryMesh_class.f90)'

    ! Choose whether to display messages
    loud = .true.
    if (present(silent)) loud = .not. silent

    ! Print beginning
    if (loud) then
      print *, repeat('<>', MAX_COL / 2)
      print *, "/\/\ READING GEOMETRY /\/\"
      print *, "Importing Meshes"

    end if

    ! Build meshes.
    meshesDict => dict % getDictPtr('meshes')
    call self % meshes % init(meshesDict, mats)
    if (loud) print *, "DONE!"

    ! Loop through all meshes and count the number of localIds in each.
    nMeshes = self % meshes % getSize()
    allocate(self % elementIdOffsets(nMeshes), self % localIdOffsets(nMeshes))
    nActiveElements = 0
    nLocalIds = 0
    self % elementIdOffsets(1) = 0
    self % localIdOffsets(1) = 0
    do i = 1, nMeshes
      box = self % meshes % getMeshBox(i)
      nActiveElements = nActiveElements + box % ptr % getElementsNumber(.true.)
      nLocalIdsInMesh = box % ptr % getLocalIdsNumber()
      nLocalIds = nLocalIds + nLocalIdsInMesh
      
      if (1 < i) then
        iLessOne = i - 1
        self % elementIdOffsets(i) = self % elementIdOffsets(iLessOne) + box % ptr % getElementsNumber()
        self % localIdOffsets(i) = self % localIdOffsets(iLessOne) + nLocalIdsInMesh

      end if

    end do
    allocate(self % cumulativeToActiveElementsMap(2, nActiveElements), self % cumulativeVolumes(nActiveElements), &
             self % fills(nLocalIds))

    ! Loop through all mesh subdictionaries and build fills and boundary conditions.
    if (loud) print *, 'Building fills.'
    cumulativeVolumeIdx = 0
    fillIdx = 0
    do i = 1, nMeshes
      box = self % meshes % getMeshBox(i)

      ! Get subdictionary corresponding to the current mesh.
      meshDict => meshesDict % getDictPtr(box % name)

      ! Get fills from dictionary.
      call meshDict % get(meshFills, 'fills')

      ! Check that the number of fills matches the number of localIds in the current mesh.
      nFills = size(meshFills)
      if (nFills /= box % ptr % getLocalIdsNumber()) call fatalError(here, &
      'The number of fills does not match the number of local ids for mesh with id: '//numToChar(box % ptr % getId())//'.')

      ! Convert all fills to materials.
      do j = 1, nFills
        fill = mats % getOrDefault(meshFills(j), NOT_FOUND)
        if (fill == NOT_FOUND) call fatalError(here, 'Unknown material: '//trim(meshFills(j))//'.')

        fillIdx = fillIdx + 1
        self % fills(fillIdx) = fill

      end do

      ! Append cumulative volumes and total volume.
      do j = 1, box % ptr % getElementsNumber()
        if (.not. box % ptr % getElementIsActive(j)) cycle
        volume = box % ptr % getElementVolume(j)
        self % totalVolume = self % totalVolume + volume
        cumulativeVolumeIdx = cumulativeVolumeIdx + 1
        self % cumulativeVolumes(cumulativeVolumeIdx) = self % totalVolume
        self % cumulativeToActiveElementsMap(:, cumulativeVolumeIdx) = [i, j]

      end do

    end do

    if (loud) print *, 'DONE!'

    ! Print geometry information and end
    if (loud) then
      print *, "GEOMETRY INFORMATION "
      print '(2X, 2A)', "Number of Meshes: ", numToChar(nMeshes)
      print *, "\/\/ FINISHED READING GEOMETRY \/\/"
      print *, repeat('<>', MAX_COL / 2)

    end if

  end subroutine init

  !!
  !!
  !!
  subroutine kill(self)
    class(geometryMesh), intent(inout) :: self

    ! Local.
    if (allocated(self % elementIdOffsets)) deallocate(self % elementIdOffsets)
    if (allocated(self % fills)) deallocate(self % fills)
    if (allocated(self % localIdOffsets)) deallocate(self % localIdOffsets)
    self % totalVolume = ZERO
    if (allocated(self % cumulativeVolumes)) deallocate(self % cumulativeVolumes)
    call self % meshes % kill()

  end subroutine kill

  !!
  !!
  !!
  subroutine move(self, coords, maxDist, event, cache)
    class(geometryMesh), intent(in)          :: self
    type(coordList), intent(inout)           :: coords
    real(defReal), intent(inout)             :: maxDist
    integer(shortInt), intent(out)           :: event
    type(distCache), intent(inout), optional :: cache
    class(mesh), pointer                     :: meshPtr
    integer(shortInt)                        :: faceIdx
    integer(shortInt), dimension(N_BC_TYPES) :: boundaryConditions
    type(coordData)                          :: updateData
    character(*), parameter                  :: here = 'move (geometryMesh_class.f90)'

    if (.not. coords % isPlaced()) call fatalError(Here, 'Coordinate list is not placed in the geometry.')

    ! Find distance to the next mesh face. First retrieve the pointer of the mesh occupied by the coordinates.
    updateData = coords % getCoordinatesData(1)
    updateData % dMax = maxDist
    meshPtr => self % meshes % getMeshPtr(updateData % meshIdx)

    ! Find the distance to the next mesh face intersection.
    call meshPtr % distanceToNextFace(updateData)

    if (maxDist < updateData % d) then ! Moves within cell
      ! Move local, register event and return early
      call coords % moveLocal(maxDist, coords % getNesting())
      event = COLL_EV
      return

    end if

    ! If reached this point a face was hit. Update maxDist and move coordinates to hit position.
    maxDist = updateData % d
    call coords % moveLocal(maxDist, coords % getNesting())

    ! Check if the face that was hit is a boundary face.
    if (meshPtr % getFaceIsBoundary(updateData % faceIdx)) then
      ! We have hit a boundary. Check if boundary condition is void and return early if so. Else, apply boundary conditions.
      event = BOUNDARY_EV
      boundaryConditions = meshPtr % getFaceBoundaryConditions(updateData % faceIdx)
      if (boundaryConditions(TRANSPORT_BCs) == VACUUM_BC) then
        call coords % setMatIdx(OUTSIDE_FILL)
        call coords % setUniqueId(0)

      else
        updateData % r = coords % getPosition(1)
        call meshPtr % explicitBoundaryConditions(updateData % faceIdx, TRANSPORT_BCs, updateData)
        call coords % updateCoordinatesFromData(1, updateData)

      end if

    else
      ! Particle moves from element to element. Update coordinates, then find materialIdx and uniqueId.
      event = CROSS_EV
      updateData % r = coords % getPosition(1)
      call coords % updateCoordinatesFromData(1, updateData)
      call coords % setMatIdx(self % fills(self % localIdOffsets(updateData % meshIdx) + updateData % localId))
      call coords % setUniqueId(self % elementIdOffsets(updateData % meshIdx) + updateData % elementIdx)

    end if

  end subroutine move

  !!
  !!
  !!
  subroutine moveGlobal(self, coords, maxDist, event)
    class(geometryMesh), intent(in) :: self
    type(coordList), intent(inout)  :: coords
    real(defReal), intent(inout)    :: maxDist
    integer(shortInt), intent(out)  :: event
    character(*), parameter         :: here = 'moveGlobal (geometryMesh_class.f90)'

    call fatalError(here, 'Unsupported procedure.')

  end subroutine moveGlobal

  !!
  !!
  !!
  subroutine placeCoord(self, coords)
    class(geometryMesh), intent(in) :: self
    type(coordList), intent(inout)  :: coords
    integer(shortInt)               :: nesting
    type(coordData)                 :: data
    character(*), parameter         :: here = 'placeCoord (geometryMesh_class.f90)'

    ! Check that coordList is initialised.
    nesting = coords % getNesting()
    if (nesting < 1) call fatalError(here, 'CoordList is not initialised. Nesting is: '//numToChar(nesting)//'.')

    ! Place coordinates above geometry (in case they were placed).
    call coords % takeAboveGeom()

    ! Find host mesh and its element.
    data = newCoordData(coords % getPosition(1), coords % getDirection(1))
    call self % meshes % findHostMesh(data)

    ! Set new coordinates in the list.
    call coords % updateCoordinatesFromData(1, data)

    ! Get material corresponding to localId in mesh.
    if (0 < data % meshIdx) then
      call coords % setMatIdx(self % fills(self % localIdOffsets(data % meshIdx) + data % localId))
      call coords % setUniqueId(self % elementIdOffsets(data % meshIdx) + data % elementIdx)

    end if

  end subroutine placeCoord

  !!
  !!
  !!
  subroutine sampleInitialPosition(self, bottom, top, rand, materialIdx, uniqueId, r, temperature)
    class(geometryMesh), intent(in)          :: self
    real(defReal), dimension(3), intent(in)  :: bottom, top
    class(RNG), intent(inout)                :: rand
    integer(shortInt), intent(out)           :: materialIdx, uniqueId
    real(defReal), dimension(3), intent(out) :: r
    real(defReal), intent(out), optional     :: temperature
    class(mesh), pointer                     :: meshPtr
    class(scalarField), pointer              :: temperatureFieldPtr
    integer(shortInt)                        :: cumulativeIdx, elementIdx, localId, meshIdx
    real(defReal)                            :: randomNumber
    type(coordList)                          :: coords

    ! First sample an element at random. Use cumulative volumes for sampling.
    call rand % generate(randomNumber, self % totalVolume)
    cumulativeIdx = ceilingBinarySearch(self % cumulativeVolumes, randomNumber)

    ! Get mesh corresponding to cumulativeIdx.
    meshIdx = self % cumulativeToActiveElementsMap(1, cumulativeIdx)
    meshPtr => self % meshes % getMeshPtr(meshIdx)

    ! Now sample a position at random within the element in the mesh.
    elementIdx = self % cumulativeToActiveElementsMap(2, cumulativeIdx)
    call meshPtr % sampleInitialPosition(elementIdx, rand, localId, r)

    ! Return materialIdx and uniqueId.
    materialIdx = self % fills(self % localIdOffsets(meshIdx) + localId)
    uniqueId = self % elementIdOffsets(meshIdx) + elementIdx

    ! Get temperature if requested.
    if (present(temperature)) then
      temperature = ZERO
      temperatureFieldPtr => getTemperatureFieldPtr()
      if (associated(temperatureFieldPtr)) then
        call coords % setElementIdx(elementIdx, 1)
        call coords % setMeshIdx(meshIdx, 1)
        temperature = temperatureFieldPtr % at(coords)

      end if

    end if

  end subroutine sampleInitialPosition

  !!
  !!
  !!
  subroutine teleport(self, coords, dist)
    class(geometryMesh), intent(in) :: self
    type(coordList), intent(inout)  :: coords
    real(defReal), intent(in)       :: dist
    character(*), parameter         :: here = 'teleport (geometryMesh_class.f90)'

    call fatalError(here, 'Unsupported procedure.')

  end subroutine teleport

  !!
  !!
  !!
  subroutine whatIsAt(self, matIdx, uniqueId, r, u, temperature)
    class(geometryMesh), intent(in)                   :: self
    integer(shortInt), intent(out)                    :: matIdx, uniqueID
    real(defReal), dimension(3), intent(in)           :: r
    real(defReal), dimension(3), optional, intent(in) :: u
    real(defReal), intent(out), optional              :: temperature
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

end module geometryMesh_class