module octreeOptimisationPackage_class

  use coord_class,            only : coord
  use dictionary_class,       only : dictionary
  use genericProcedures,      only : numToChar
  use geometry_inter,         only : geometry
  use geometryFactory_func,   only : new_geometry
  use geometryReg_mod,        only : geomIdx, geomPtr
  use geometryStd_class,      only : geometryStd
  use hashFunctions_func,     only : FNV_1
  use meshUniverse_class,     only : meshUniverse
  use nuclearDataReg_mod,     only : ndReg_init => init
  use numPrecision
  use outputFile_class,       only : outputFile
  use physicsPackage_inter,   only : physicsPackage
  use rng_class,              only : RNG
  use universe_inter,         only : universe
  use unstructuredMesh_inter, only : unstructuredMesh

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(physicsPackage) :: octreeOptimisationPackage
    private
    character(pathLen)                             :: outputFile = '', outputFormat = ''
    class(geometry), pointer                       :: geom => null()
    integer(longInt)                               :: seed = 0_longInt
    integer(shortInt)                              :: nRuns = 0, population = 0
    integer(shortInt), dimension(:), allocatable   :: depths, nMaxFaces
    real(defReal), dimension(:, :), allocatable    :: storageSizes
    real(defReal), dimension(:, :, :), allocatable :: hostTimes, initialisationTimes
  contains
    procedure :: collectResults
    procedure :: init
    procedure :: run
    procedure :: runSingle
  end type octreeOptimisationPackage

contains
  !!
  !!
  !!
  subroutine collectResults(self)
    class(octreeOptimisationPackage), intent(in) :: self
    character(nameLen)                           :: name, suffix
    integer(shortInt)                            :: i, j, k
    type(outputFile)                             :: out

    call out % init(self % outputFormat, filename = self % outputFile)

    name = 'seed'
    call out % printValue(self % seed, name)

    name = 'pop'
    call out % printValue(self % population, name)

    name = 'nRuns'
    call out % printValue(self % nRuns, name)

    ! Print values.
    name = 'depths'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, [size(self % depths)])
    do i = 1, size(self % depths)
      call out % addValue(self % depths(i))

    end do
    call out % endArray()
    call out % endBlock()

    name = 'nMaxFaces'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, [size(self % nMaxFaces)])
    do i = 1, size(self % nMaxFaces)
      call out % addValue(self % nMaxFaces(i))

    end do
    call out % endArray()
    call out % endBlock()

    do i = 1, size(self % depths)
      do j = 1, size(self % nMaxFaces)
        suffix = '_D'//numToChar(self % depths(i))//'_N'//numToChar(self % nMaxFaces(j))

        name = 'storageSize'//trim(suffix)
        call out % printValue(self % storageSizes(i, j), name)

        name = 'rawInitTimes'//trim(suffix)
        call out % startBlock(name)
        name = 'Res'
        call out % startArray(name, [self % nRuns])
        do k = 1, self % nRuns
          call out % addValue(self % initialisationTimes(i, j, k))

        end do
        call out % endArray()
        call out % endBlock()

        name = 'rawHostTimes'//trim(suffix)
        call out % startBlock(name)
        name = 'Res'
        call out % startArray(name, [self % nRuns])
        do k = 1, self % nRuns
          call out % addValue(self % hostTimes(i, j, k))

        end do
        call out % endArray()
        call out % endBlock()

      end do

    end do

  end subroutine collectResults

  !!
  !!
  !!
  subroutine init(self, dict)
    class(octreeOptimisationPackage), intent(inout) :: self
    class(dictionary), intent(inout)                :: dict
    character(8)                                    :: date
    character(10)                                   :: time
    character(nameLen)                              :: geometryName
    character(:), allocatable                       :: string
    integer(shortInt)                               :: seed_temp, sizeDepths, sizeNMaxFaces
    type(outputFile)                                :: testOutput

    ! Read outputfile path
    call dict % getOrDefault(self % outputFile, 'outputFile', './output')

    ! Get output format and verify
    ! Initialise output file before calculation (so mistake in format will be caught early)
    call dict % getOrDefault(self % outputFormat, 'outputFormat', 'asciiMATLAB')
    call testOutput % init(self % outputFormat)

    ! Read parameters.
    call dict % get(self % depths, 'depths')
    call dict % get(self % nMaxFaces, 'nMaxFaces')
    call dict % getOrDefault(self % nRuns, 'nRuns', 20)

    ! Allocate memory.
    sizeDepths = size(self % depths)
    sizeNMaxFaces = size(self % nMaxFaces)
    allocate(self % storageSizes(sizeDepths, sizeNMaxFaces), &
             self % hostTimes(sizeDepths, sizeNMaxFaces, self % nRuns), &
             self % initialisationTimes(sizeDepths, sizeNMaxFaces, self % nRuns))

    ! Read calculation settings
    call dict % get(self % population, 'pop')

    ! *** It is a bit silly but dictionary cannot store longInt for now
    !     so seeds are limited to 32 bits (can be -ve)
    if (dict % isPresent('seed')) then
      call dict % get(seed_temp, 'seed')

    else
      ! Obtain time string and hash it to obtain random seed
      call date_and_time(date, time)
      string = date//time
      call FNV_1(string, seed_temp)

    end if
    self % seed = seed_temp

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
    class(octreeOptimisationPackage), intent(inout) :: self
    class(universe), pointer                        :: universePtr
    class(unstructuredMesh), pointer                :: unstructuredMeshPtr
    integer(shortInt)                               :: i, j, k
    real(defReal), dimension(3)                     :: u
    real(defReal), dimension(6)                     :: bounds
    type(coord)                                     :: coords
    type(geometryStd), pointer                      :: geometryStdPtr
    type(meshUniverse), pointer                     :: meshUniversePtr

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

    ! Loop over all runs, then all nMaxFaces, then all depths.
    do i = 1, size(self % depths)
      do j = 1, size(self % nMaxFaces)
        do k = 1, self % nRuns
          call self % runSingle(i, j, k, bounds, unstructuredMeshPtr, coords)

        end do

      end do

    end do

    call self % collectResults()

  end subroutine run

  !!
  !!
  !!
  subroutine runSingle(self, depthNumber, nMaxFacesNumber, runNumber, bounds, uMesh, coords)
    class(octreeOptimisationPackage), intent(inout) :: self
    integer(shortInt), intent(in)                   :: depthNumber, nMaxFacesNumber, runNumber
    real(defReal)                                   :: t1, t2
    real(defReal), dimension(6), intent(in)         :: bounds
    class(unstructuredMesh), intent(inout)          :: uMesh
    type(coord), intent(inout)                      :: coords
    integer(shortInt)                               :: i
    real(defReal), dimension(3)                     :: boundsDifference, randomNumbers
    type(dictionary)                                :: dict
    type(RNG)                                       :: pRNG

    ! Initialise dictionary then fill it.
    call dict % init(3)
    call dict % store('type', 'octreeAcceleration')
    call dict % store('depth', self % depths(depthNumber))
    call dict % store('nMaxFaces', self % nMaxFaces(nMaxFacesNumber))

    ! Now initialise acceleration structure in the mesh.
    call cpu_time(t1)
    call uMesh % initAccelerationStructure(dict)
    call cpu_time(t2)
    self % initialisationTimes(depthNumber, nMaxFacesNumber, runNumber) = t2 - t1

    ! Get storage size.
    self % storageSizes(depthNumber, nMaxFacesNumber) = uMesh % getAccelerationStructureStorageSize() / 1.0e6_defReal

    boundsDifference = bounds(4:6) - bounds(1:3)
    call pRNG % init(self % seed)
    call cpu_time(t1)
    do i = 1, self % population
      ! Sample Position
      call pRNG % generate(randomNumbers)
      call coords % setPosition(boundsDifference * randomNumbers + bounds(1:3))

      ! Find element occupied by coordinates.
      call uMesh % findHostElement(coords)

    end do
    call cpu_time(t2)
    self % hostTimes(depthNumber, nMaxFacesNumber, runNumber) = t2 - t1

    call uMesh % killAccelerationStructure()

  end subroutine runSingle

end module octreeOptimisationPackage_class