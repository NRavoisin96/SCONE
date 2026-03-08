module octreeOptimisationPackage_class

  use coord_class,            only : coord
  use dictionary_class,       only : dictionary
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
    real(defReal), dimension(:, :), allocatable    :: averageHostTimes, averageInitialisationTimes, averageTotalTimes, &
                                                      stdHostTimes, stdInitialisationTimes, stdTotalTimes
    real(defReal), dimension(:, :, :), allocatable :: hostTimes, initialisationTimes, totalTimes
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
    character(nameLen)                           :: name
    integer(shortInt)                            :: i, j
    integer(shortInt), dimension(:), allocatable :: arrayShape
    type(outputFile)                             :: out

    call out % init(self % outputFormat, filename = self % outputFile)

    name = 'seed'
    call out % printValue(self % seed, name)

    name = 'pop'
    call out % printValue(self % population, name)

    ! Print values.
    arrayShape = [size(self % depths), size(self % nMaxFaces)]
    name = 'averageInitialisationTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do j = 1, size(self % nMaxFaces)
      do i = 1, size(self % depths)
        call out % addValue(self % averageInitialisationTimes(i, j))

      end do

    end do
    call out % endArray()
    call out % endBlock()

    name = 'stdInitialisationTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do j = 1, size(self % nMaxFaces)
      do i = 1, size(self % depths)
        call out % addValue(self % stdInitialisationTimes(i, j))

      end do

    end do
    call out % endArray()
    call out % endBlock()

    name = 'averageHostTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do j = 1, size(self % nMaxFaces)
      do i = 1, size(self % depths)
        call out % addValue(self % averageHostTimes(i, j))

      end do

    end do
    call out % endArray()
    call out % endBlock()

    name = 'stdHostTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do j = 1, size(self % nMaxFaces)
      do i = 1, size(self % depths)
        call out % addValue(self % stdHostTimes(i, j))

      end do

    end do
    call out % endArray()
    call out % endBlock()

    name = 'averageTotalTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do j = 1, size(self % nMaxFaces)
      do i = 1, size(self % depths)
        call out % addValue(self % averageTotalTimes(i, j))

      end do

    end do
    call out % endArray()
    call out % endBlock()

    name = 'stdTotalTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do j = 1, size(self % nMaxFaces)
      do i = 1, size(self % depths)
        call out % addValue(self % stdTotalTimes(i, j))

      end do

    end do
    call out % endArray()
    call out % endBlock()

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
    allocate(self % averageHostTimes(sizeDepths, sizeNMaxFaces), self % averageInitialisationTimes(sizeDepths, sizeNMaxFaces), &
             self % averageTotalTimes(sizeDepths, sizeNMaxFaces), self % hostTimes(sizeDepths, sizeNMaxFaces, self % nRuns), &
             self % stdHostTimes(sizeDepths, sizeNMaxFaces), self % stdInitialisationTimes(sizeDepths, sizeNMaxFaces), &
             self % stdTotalTimes(sizeDepths, sizeNMaxFaces), self % initialisationTimes(sizeDepths, sizeNMaxFaces, self % nRuns), &
             self % totalTimes(sizeDepths, sizeNMaxFaces, self % nRuns))

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
    real(defReal)                                   :: sumInit, sumHost, sumTot, sumofSquaresInit, sumOfSquaresHost, &
                                                       sumOfSquaresTot, t_init, t_host, t_tot
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
    self % totalTimes = self % initialisationTimes + self % hostTimes

    ! Compute means and variances.
    do i = 1, size(self % depths)
      do j = 1, size(self % nMaxFaces)
        sumInit = ZERO
        sumHost = ZERO
        sumTot = ZERO
        sumOfSquaresInit = ZERO
        sumOfSquaresHost = ZERO
        sumOfSquaresTot = ZERO
        do k = 1, self % nRuns
          t_init = self % initialisationTimes(i, j, k)
          sumInit = sumInit + t_init
          sumofSquaresInit = sumofSquaresInit + t_init * t_init

          t_host = self % hostTimes(i, j, k)
          sumHost = sumHost + t_host
          sumOfSquaresHost = sumOfSquaresHost + t_host * t_host

          t_tot = self % totalTimes(i, j, k)
          sumTot = sumTot + t_tot
          sumOfSquaresTot = sumOfSquaresTot + t_tot * t_tot

        end do
        self % averageInitialisationTimes(i, j) = sumInit / self % nRuns
        self % stdInitialisationTimes(i, j) = sqrt((sumofSquaresInit - sumInit * sumInit / self % nRuns) / (self % nRuns - 1))

        self % averageHostTimes(i, j) = sumHost / self % nRuns
        self % stdHostTimes(i, j) = sqrt((sumofSquaresHost - sumHost * sumHost / self % nRuns) / (self % nRuns - 1))

        self % averageTotalTimes(i, j) = sumTot / self % nRuns
        self % stdTotalTimes(i, j) = sqrt((sumofSquaresTot - sumTot * sumTot / self % nRuns) / (self % nRuns - 1))

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

    print *, 'Depth: ', self % depths(depthNumber)
    print *, 'nMaxFaces: ', self % nMaxFaces(nMaxFacesNumber)
    print *, 'Run: ', runNumber

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
    print*, "-------------------------------------------------------------"
    print*, "/\/\ Initialisation procedure time /\/\"
    print*, "CPU time: ", self % initialisationTimes(depthNumber, nMaxFacesNumber, runNumber), " seconds"
    print*, "-------------------------------------------------------------"

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
    print*, "-------------------------------------------------------------"
    print*, "/\/\ Host element determination time /\/\"
    print*, "CPU time: ", self % hostTimes(depthNumber, nMaxFacesNumber, runNumber), " seconds"
    print*, "-------------------------------------------------------------"

    call uMesh % killAccelerationStructure()

  end subroutine runSingle

end module octreeOptimisationPackage_class