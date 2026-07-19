module hostElementDeterminationPackage_class

  use coord_class,            only : coord
  use dictionary_class,       only : dictionary
  use errors_mod,             only : fatalError
  use genericProcedures,      only : numToChar
  use geometry_inter,         only : geometry
  use geometryStd_class,      only : geometryStd
  use geometryFactory_func,   only : new_geometry
  use geometryReg_mod,        only : geomIdx, geomPtr
  use meshUniverse_class,     only : meshUniverse
  use nuclearDataReg_mod,     only : ndReg_init => init
  use numPrecision
  use outputFile_class,       only : outputFile
  use physicsPackage_inter,   only : physicsPackage
  use rng_class,              only : rng
  use timer_mod,              only : registerTimer
  use universe_inter,         only : universe
  use unstructuredMesh_inter, only : unstructuredMesh

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(physicsPackage) :: hostElementDeterminationPackage
    private
    character(nameLen)                           :: outputFormat = '', patchType = ''
    character(pathLen)                           :: outputFile = ''
    class(geometry), pointer                     :: geom => null()
    integer(shortInt)                            :: octreeDepth = 0, octreeNMaxFaces = 0, pop = 0
    integer(shortInt), dimension(:), allocatable :: depths, seeds
    logical(defBool)                             :: singleFaceShortcut = .true.
    real(defReal)                                :: octreeStorageSize = ZERO
    real(defReal), dimension(:), allocatable     :: averageHostTimes_patch, averageInitialisationTimes_patch, hostTimes_other, &
                                                    initialisationTimes_other, patchSearchStorageSizes
    real(defReal), dimension(:, :), allocatable  :: hostTimes_patch, initialisationTimes_patch
  contains
    procedure :: collectResults
    procedure :: init
    procedure :: run
    procedure :: runSingle_other
    procedure :: runSingle_patch
  end type hostElementDeterminationPackage

contains
  !!
  !!
  !!
  subroutine collectResults(self)
    class(hostElementDeterminationPackage), intent(in) :: self
    character(nameLen)                                 :: name
    integer(shortInt)                                  :: i, j
    type(outputFile)                                   :: out

    call out % init(self % outputFormat, filename = self % outputFile)

    name = 'pop'
    call out % printValue(self % pop, name)

    name = 'patchType'
    call out % printValue(self % patchType, name)

    name = 'singleFaceShortcut'
    call out % printValue(merge(1, 0, self % singleFaceShortcut), name)

    name = 'octreeDepth'
    call out % printValue(self % octreeDepth, name)

    name = 'octreeNMaxFaces'
    call out % printValue(self % octreeNMaxFaces, name)

    name = 'seeds'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, [size(self % seeds)])
    do i = 1, size(self % seeds)
      call out % addValue(self % seeds(i))

    end do
    call out % endArray()
    call out % endBlock()

    ! Print values for octree.
    name = 'octreeStorageSize'
    call out % printValue(self % octreeStorageSize, name)

    name = 'rawOctreeInitTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, [size(self % seeds)])
    do i = 1, size(self % seeds)
      call out % addValue(self % initialisationTimes_other(i))

    end do
    call out % endArray()
    call out % endBlock()

    name = 'rawOctreeHostTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, [size(self % seeds)])
    do i = 1, size(self % seeds)
      call out % addValue(self % hostTimes_other(i))

    end do
    call out % endArray()
    call out % endBlock()

    ! Print values for Patch-Search.
    do i = 1, size(self % depths)
      name = 'patchStorageSize_D'//numToChar(self % depths(i))
      call out % printValue(self % patchSearchStorageSizes(i), name)

      name = 'rawPatchInitTimes_D'//numToChar(self % depths(i))
      call out % startBlock(name)
      name = 'Res'
      call out % startArray(name, [size(self % seeds)])
      do j = 1, size(self % seeds)
        call out % addValue(self % initialisationTimes_patch(i, j))

      end do
      call out % endArray()
      call out % endBlock()

      name = 'rawPatchHostTimes_D'//numToChar(self % depths(i))
      call out % startBlock(name)
      name = 'Res'
      call out % startArray(name, [size(self % seeds)])
      do j = 1, size(self % seeds)
        call out % addValue(self % hostTimes_patch(i, j))

      end do
      call out % endArray()
      call out % endBlock()

    end do

  end subroutine collectResults

  !!
  !!
  !!
  subroutine init(self, dict)
    class(hostElementDeterminationPackage), intent(inout) :: self
    class(dictionary), intent(inout)                      :: dict
    character(nameLen)                                    :: geometryName
    integer(shortInt)                                     :: nDepths, nRuns
    type(outputFile)                                      :: testOutput
    character(*), parameter                               :: HERE = 'init (hostElementDeterminationPackage_class.f90)'

    ! Retrieve patch type.
    call dict % getOrDefault(self % patchType, 'patchType', 'patchSearchAcceleration')

    ! Retrieve singleFaceShortcut from dictionary.
    call dict % getOrDefault(self % singleFaceShortcut, 'singleFaceShortcut', .true.)

    ! Read outputfile path
    call dict % getOrDefault(self % outputFile, 'outputFile', './output')

    ! Get output format and verify
    ! Initialise output file before calculation (so mistake in format will be caught early)
    call dict % getOrDefault(self % outputFormat, 'outputFormat', 'asciiMATLAB')
    call testOutput % init(self % outputFormat)

    ! Read parameters.
    call dict % get(self % octreeDepth, 'octreeDepth')
    call dict % get(self % octreeNMaxFaces, 'octreeNMaxFaces')
    call dict % get(self % depths, 'depths')
    call dict % get(self % seeds, 'seeds')
    call dict % get(self % pop, 'pop')

    ! Allocate memory.
    nDepths = size(self % depths)
    nRuns = size(self % seeds)
    allocate(self % hostTimes_patch(nDepths, nRuns), self % hostTimes_other(nRuns), &
             self % initialisationTimes_patch(nDepths, nRuns), self % initialisationTimes_other(nRuns), &
             self % patchSearchStorageSizes(nDepths))

    ! Build Nuclear Data
    call ndReg_init(dict % getDictPtr("nuclearData"))

    ! Build geometry
    geometryName = 'testGeometry'
    call new_geometry(dict % getDictPtr('geometry'), geometryName)
    self % geom => geomPtr(geomIdx(geometryName))

  end subroutine init

  !!
  !!
  !!
  subroutine run(self)
    class(hostElementDeterminationPackage), intent(inout) :: self
    class(universe), pointer                              :: universePtr
    class(unstructuredMesh), pointer                      :: unstructuredMeshPtr
    integer(shortInt)                                     :: i, j, nRuns
    real(defReal), dimension(3)                           :: u
    real(defReal), dimension(6)                           :: bounds
    type(coord)                                           :: coords
    type(geometryStd), pointer                            :: geometryStdPtr
    type(meshUniverse), pointer                           :: meshUniversePtr

    ! Downcast.
    select type(temp => self % geom)
      type is(geometryStd)
        geometryStdPtr => temp

    end select

    ! Get pointer to mesh universe then downcast.
    universePtr => geometryStdPtr % geom % unis % getPtr_fast(2)
    select type(temp => universePtr)
      type is(meshUniverse)
        meshUniversePtr => temp

    end select

    ! Get pointer to the mesh in the mesh universe then downcast.
    select type(temp => meshUniversePtr % mesh % ptr)
      class is(unstructuredMesh)
        unstructuredMeshPtr => temp

    end select

    u = [ONE, ZERO, ZERO]
    call coords % setDirection(u)
    bounds = self % geom % bounds()
    nRuns = size(self % seeds)

    ! Do octree first.
    do i = 1, nRuns
      call self % runSingle_other(i, bounds, unstructuredMeshPtr, coords)

    end do

    ! Now do Patch-Search.
    do i = 1, size(self % depths)
      do j = 1, nRuns
        call self % runSingle_patch(i, j, bounds, unstructuredMeshPtr, coords)

      end do

    end do

    call self % collectResults()

  end subroutine run

  !!
  !!
  !!
  subroutine runSingle_other(self, runNumber, bounds, uMesh, coords)
    class(hostElementDeterminationPackage), intent(inout) :: self
    integer(shortInt), intent(in)                         :: runNumber
    real(defReal), dimension(6), intent(in)               :: bounds
    class(unstructuredMesh), intent(inout)                :: uMesh
    type(coord), intent(inout)                            :: coords
    integer(longInt)                                      :: seed
    integer(shortInt)                                     :: i
    real(defReal)                                         :: t1, t2
    real(defReal), dimension(3)                           :: boundsDifference, randomNumbers
    type(dictionary)                                      :: dict
    type(RNG)                                             :: pRNG

    ! Initialise dictionary then fill it.
    call dict % init(3)
    call dict % store('type', 'octreeAcceleration')
    call dict % store('depth', self % octreeDepth)
    call dict % store('nMaxFaces', self % octreeNMaxFaces)

    ! Now initialise acceleration structure in the mesh.
    call cpu_time(t1)
    call uMesh % initAccelerationStructure(dict)
    call cpu_time(t2)
    self % initialisationTimes_other(runNumber) = t2 - t1
    self % octreeStorageSize = uMesh % getAccelerationStructureStorageSize() / 1.0e6_defReal

    boundsDifference = bounds(4:6) - bounds(1:3)
    seed = self % seeds(runNumber)
    call pRNG % init(seed)
    call cpu_time(t1)
    do i = 1, self % pop
      ! Sample Position
      call pRNG % generate(randomNumbers)
      call coords % setPosition(boundsDifference * randomNumbers + bounds(1:3))

      ! Find element occupied by coordinates.
      call uMesh % findHostElement(coords)

    end do
    call cpu_time(t2)
    self % hostTimes_other(runNumber) = t2 - t1
    call uMesh % killAccelerationStructure()

  end subroutine runSingle_other

  !!
  !!
  !!
  subroutine runSingle_patch(self, depthNumber, runNumber, bounds, uMesh, coords)
    class(hostElementDeterminationPackage), intent(inout) :: self
    integer(shortInt), intent(in)                         :: depthNumber, runNumber
    real(defReal), dimension(6), intent(in)               :: bounds
    class(unstructuredMesh), intent(inout)                :: uMesh
    type(coord), intent(inout)                            :: coords
    integer(longInt)                                      :: seed
    integer(shortInt)                                     :: i
    real(defReal)                                         :: t1, t2
    real(defReal), dimension(3)                           :: boundsDifference, randomNumbers
    type(dictionary)                                      :: dict
    type(RNG)                                             :: pRNG

    ! Initialise dictionary then fill it.
    call dict % init(3)
    call dict % store('type', self % patchType)
    call dict % store('depth', self % depths(depthNumber))

    call dict % store('singleFaceShortcut', merge(1, 0, self % singleFaceShortcut))

    ! Now initialise acceleration structure in the mesh.
    call cpu_time(t1)
    call uMesh % initAccelerationStructure(dict)
    call cpu_time(t2)
    self % initialisationTimes_patch(depthNumber, runNumber) = t2 - t1
    self % patchSearchStorageSizes(depthNumber) = uMesh % getAccelerationStructureStorageSize() / 1.0e6_defReal

    boundsDifference = bounds(4:6) - bounds(1:3)
    seed = self % seeds(runNumber)
    call pRNG % init(seed)
    call cpu_time(t1)
    do i = 1, self % pop
      ! Sample Position
      call pRNG % generate(randomNumbers)
      call coords % setPosition(boundsDifference * randomNumbers + bounds(1:3))

      ! Find element occupied by coordinates.
      call uMesh % findHostElement(coords)

    end do
    call cpu_time(t2)
    self % hostTimes_patch(depthNumber, runNumber) = t2 - t1
    call uMesh % killAccelerationStructure()

  end subroutine runSingle_patch

end module hostElementDeterminationPackage_class