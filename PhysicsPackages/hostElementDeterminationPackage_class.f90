module hostElementDeterminationPackage_class

  use coord_class,            only : coord
  use dictionary_class,       only : dictionary
  use errors_mod,             only : fatalError
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
    character(nameLen)                           :: outputFile = '', outputFormat = ''
    class(geometry), pointer                     :: geom => null()
    integer(shortInt)                            :: octreeDepth = 0, octreeNMaxFaces = 0, pop = 0
    integer(shortInt), dimension(:), allocatable :: depths, seeds
    real(defReal)                                :: averageHostTimes_other = ZERO, averageInitialisationTimes_other = ZERO, &
                                                    averageTotalTimes_other = ZERO, stdHostTimes_other = ZERO, &
                                                    stdInitialisationTimes_other = ZERO, stdTotalTimes_other = ZERO
    real(defReal), dimension(:), allocatable     :: averageHostTimes_patch, averageInitialisationTimes_patch, &
                                                    averageTotalTimes_patch, averageSavingsHostTimes_patch, &
                                                    averageSavingsInitialisationTimes_patch, averageSavingsTotalTimes_patch, &
                                                    stdHostTimes_patch, stdInitialisationTimes_patch, stdTotalTimes_patch, &
                                                    stdSavingsHostTimes_patch, stdSavingsInitialisationTimes_patch, &
                                                    stdSavingsTotalTimes_patch, hostTimes_other, initialisationTimes_other, &
                                                    totalTimes_other, averageSavingsInitialisationTimes_patch_1, &
                                                    stdSavingsInitialisationTimes_patch_1
    real(defReal), dimension(:, :), allocatable  :: hostTimes_patch, initialisationTimes_patch, savingsHostTimes_patch, &
                                                    savingsInitialisationTimes_patch, savingsTotalTimes_patch, totalTimes_patch
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
    integer(shortInt)                                  :: i
    integer(shortInt), dimension(:), allocatable       :: arrayShape
    type(outputFile)                                   :: out

    call out % init(self % outputFormat, filename = self % outputFile)

    name = 'pop'
    call out % printValue(self % pop, name)

    ! Print values for octree.
    name = 'octreeMeanInitialisationTime'
    call out % printValue(self % averageInitialisationTimes_other, name)

    name = 'octreeStdInitialisationTime'
    call out % printValue(self % stdInitialisationTimes_other, name)

    name = 'octreeMeanHostTime'
    call out % printValue(self % averageHostTimes_other, name)

    name = 'octreeStdHostTime'
    call out % printValue(self % stdHostTimes_other, name)

    name = 'octreeTotalTime'
    call out % printValue(self % averageTotalTimes_other, name)

    name = 'octreeStdTotalTime'
    call out % printValue(self % stdTotalTimes_other, name)

    ! Print value for Patch-Search.
    arrayShape = [size(self % depths)]
    name = 'averageInitialisationTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do i = 1, size(self % depths)
      call out % addValue(self % averageInitialisationTimes_patch(i))

    end do

    call out % endArray()
    call out % endBlock()

    name = 'stdInitialisationTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do i = 1, size(self % depths)
      call out % addValue(self % stdInitialisationTimes_patch(i))

    end do
    call out % endArray()
    call out % endBlock()

    name = 'averageHostTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do i = 1, size(self % depths)
      call out % addValue(self % averageHostTimes_patch(i))

    end do
    call out % endArray()
    call out % endBlock()

    name = 'stdHostTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do i = 1, size(self % depths)
      call out % addValue(self % stdHostTimes_patch(i))

    end do
    call out % endArray()
    call out % endBlock()

    name = 'averageTotalTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do i = 1, size(self % depths)
      call out % addValue(self % averageTotalTimes_patch(i))

    end do
    call out % endArray()
    call out % endBlock()

    name = 'stdTotalTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do i = 1, size(self % depths)
      call out % addValue(self % stdTotalTimes_patch(i))

    end do
    call out % endArray()
    call out % endBlock()

    name = 'savingsInitialisationTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do i = 1, size(self % depths)
      call out % addValue(self % averageSavingsInitialisationTimes_patch(i) * 100)

    end do

    call out % endArray()
    call out % endBlock()

    name = 'stdSavingsInitialisationTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do i = 1, size(self % depths)
      call out % addValue(self % stdSavingsInitialisationTimes_patch(i) * 100)

    end do
    call out % endArray()
    call out % endBlock()

    name = 'savingsHostTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do i = 1, size(self % depths)
      call out % addValue(self % averageSavingsHostTimes_patch(i) * 100)

    end do
    call out % endArray()
    call out % endBlock()

    name = 'stdSavingsHostTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do i = 1, size(self % depths)
      call out % addValue(self % stdSavingsHostTimes_patch(i) * 100)

    end do
    call out % endArray()
    call out % endBlock()

    name = 'savingsTotalTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do i = 1, size(self % depths)
      call out % addValue(self % averageSavingsTotalTimes_patch(i) * 100)

    end do
    call out % endArray()
    call out % endBlock()

    name = 'stdSavingsTotalTimes'
    call out % startBlock(name)
    name = 'Res'
    call out % startArray(name, arrayShape)
    do i = 1, size(self % depths)
      call out % addValue(self % stdSavingsTotalTimes_patch(i) * 100)

    end do
    call out % endArray()
    call out % endBlock()

    if(any(self % depths == 1)) then
      arrayShape = [size(self % depths) - 1]
      name = 'savingsTotalTimes_1'
      call out % startBlock(name)
      name = 'Res'
      call out % startArray(name, arrayShape)
      do i = 1, size(self % depths) - 1
        call out % addValue(self % averageSavingsInitialisationTimes_patch_1(i) * 100)

      end do
      call out % endArray()
      call out % endBlock()

      name = 'stdSavingsTotalTimes_1'
      call out % startBlock(name)
      name = 'Res'
      call out % startArray(name, arrayShape)
      do i = 1, size(self % depths) - 1
        call out % addValue(self % stdSavingsInitialisationTimes_patch_1(i) * 100)

      end do
      call out % endArray()
      call out % endBlock()

    end if

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
    allocate(self % averageHostTimes_patch(nDepths), self % averageInitialisationTimes_patch(nDepths), &
             self % averageTotalTimes_patch(nDepths), self % averageSavingsHostTimes_patch(nDepths), &
             self % averageSavingsInitialisationTimes_patch(nDepths), self % averageSavingsTotalTimes_patch(nDepths), &
             self % stdHostTimes_patch(nDepths), self % stdInitialisationTimes_patch(nDepths), &
             self % stdTotalTimes_patch(nDepths), self % stdSavingsHostTimes_patch(nDepths), &
             self % stdSavingsInitialisationTimes_patch(nDepths), self % stdSavingsTotalTimes_patch(nDepths), &
             self % hostTimes_patch(nDepths, nRuns), self % initialisationTimes_patch(nDepths, nRuns), &
             self % savingsHostTimes_patch(nDepths, nRuns), self % savingsInitialisationTimes_patch(nDepths, nRuns), &
             self % savingsTotalTimes_patch(nDepths, nRuns), self % totalTimes_patch(nDepths, nRuns), &
             self % hostTimes_other(nRuns), self % initialisationTimes_other(nRuns), &
             self % totalTimes_other(nRuns), self % averageSavingsInitialisationTimes_patch_1(nDepths - 1), &
             self % stdSavingsInitialisationTimes_patch_1(nDepths - 1))

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
    real(defReal)                                         :: sumInit, sumHost, sumTot, sumofSquaresInit, sumOfSquaresHost, &
                                                             sumOfSquaresTot, sumSavingsInit, sumSavingsInit_1, sumSavingsHost, &
                                                             sumSavingsTot, sumOfSquaresSavingsInit, sumOfSquaresSavingsInit_1, &
                                                             sumOfSquaresSavingsHost, sumOfSquaresSavingsTot, t_init, t_host, &
                                                             t_tot, t_savings_init, t_savings_init_1, t_savings_host, &
                                                             t_savings_tot, inverseNRuns, inverseNRunsLessOne
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
    inverseNRuns = ONE / nRuns
    inverseNRunsLessOne = ONE / (nRuns - 1)

    ! Do octree first.
    do i = 1, nRuns
      call self % runSingle_other(i, bounds, unstructuredMeshPtr, coords)

    end do
    self % totalTimes_other = self % initialisationTimes_other + self % hostTimes_other

    ! Compute statistics for octree.
    sumInit = ZERO
    sumHost = ZERO
    sumTot = ZERO
    sumofSquaresInit = ZERO
    sumOfSquaresHost = ZERO
    sumOfSquaresTot = ZERO
    do i = 1, nRuns
      t_init = self % initialisationTimes_other(i)
      sumInit = sumInit + t_init
      sumofSquaresInit = sumofSquaresInit + t_init * t_init

      t_host = self % hostTimes_other(i)
      sumHost = sumHost + t_host
      sumOfSquaresHost = sumOfSquaresHost + t_host * t_host

      t_tot = self % totalTimes_other(i)
      sumTot = sumTot + t_tot
      sumOfSquaresTot = sumOfSquaresTot + t_tot * t_tot

    end do
    self % averageInitialisationTimes_other = sumInit * inverseNRuns
    self % stdInitialisationTimes_other = sqrt((sumofSquaresInit - sumInit * sumInit * inverseNRuns) * inverseNRunsLessOne)
    
    self % averageHostTimes_other = sumHost * inverseNRuns
    self % stdHostTimes_other = sqrt((sumofSquaresHost - sumHost * sumHost * inverseNRuns) * inverseNRunsLessOne)
    
    self % averageTotalTimes_other = sumTot * inverseNRuns
    self % stdTotalTimes_other = sqrt((sumofSquaresTot - sumTot * sumTot * inverseNRuns) * inverseNRunsLessOne)

    ! Now do Patch-Search.
    do i = 1, size(self % depths)
      do j = 1, nRuns
        call self % runSingle_patch(i, j, bounds, unstructuredMeshPtr, coords)

      end do

    end do
    self % totalTimes_patch = self % initialisationTimes_patch + self % hostTimes_patch

    ! Compute statistics for Patch-Search.
    do i = 1, size(self % depths)
      sumInit = ZERO
      sumHost = ZERO
      sumTot = ZERO
      sumSavingsInit = ZERO
      sumSavingsInit_1 = ZERO
      sumSavingsHost = ZERO
      sumSavingsTot = ZERO
      
      sumOfSquaresInit = ZERO
      sumOfSquaresSavingsInit_1 = ZERO
      sumOfSquaresHost = ZERO
      sumOfSquaresTot = ZERO
      sumOfSquaresSavingsInit = ZERO
      sumOfSquaresSavingsHost = ZERO
      sumOfSquaresSavingsTot = ZERO
      do j = 1, nRuns
        t_init = self % initialisationTimes_patch(i, j)
        sumInit = sumInit + t_init
        sumOfSquaresInit = sumOfSquaresInit + t_init * t_init

        t_host = self % hostTimes_patch(i, j)
        sumHost = sumHost + t_host
        sumOfSquaresHost = sumOfSquaresHost + t_host * t_host

        t_tot = self % totalTimes_patch(i, j)
        sumTot = sumTot + t_tot
        sumOfSquaresTot = sumOfSquaresTot + t_tot * t_tot

        t_savings_init = (t_init - self % initialisationTimes_other(j)) / self % initialisationTimes_other(j)
        sumSavingsInit = sumSavingsInit + t_savings_init
        sumOfSquaresSavingsInit = sumOfSquaresSavingsInit + t_savings_init * t_savings_init

        t_savings_host = (t_host - self % hostTimes_other(j)) / self % hostTimes_other(j)
        sumSavingsHost = sumSavingsHost + t_savings_host
        sumOfSquaresSavingsHost = sumOfSquaresSavingsHost + t_savings_host * t_savings_host

        t_savings_tot = (t_tot - self % totalTimes_other(j)) / self % totalTimes_other(j)
        sumSavingsTot = sumSavingsTot + t_savings_tot
        sumOfSquaresSavingsTot = sumOfSquaresSavingsTot + t_savings_tot * t_savings_tot

        if(1 < i .and. any(self % depths == 1)) then
          t_savings_init_1 = (t_init - self % initialisationTimes_patch(1, j)) / self % initialisationTimes_patch(1, j)
          sumSavingsInit_1 = sumSavingsInit_1 + t_savings_init_1
          sumOfSquaresSavingsInit_1 = sumOfSquaresSavingsInit_1 + t_savings_init_1 * t_savings_init_1

        end if

      end do
      self % averageInitialisationTimes_patch(i) = sumInit * inverseNRuns
      self % stdInitialisationTimes_patch(i) = sqrt((sumofSquaresInit - sumInit * sumInit * inverseNRuns) * inverseNRunsLessOne)

      self % averageHostTimes_patch(i) = sumHost * inverseNRuns
      self % stdHostTimes_patch(i) = sqrt((sumofSquaresHost - sumHost * sumHost * inverseNRuns) * inverseNRunsLessOne)

      self % averageTotalTimes_patch(i) = sumTot * inverseNRuns
      self % stdTotalTimes_patch(i) = sqrt((sumofSquaresTot - sumTot * sumTot * inverseNRuns) * inverseNRunsLessOne)

      self % averageSavingsInitialisationTimes_patch(i) = sumSavingsInit * inverseNRuns
      self % stdSavingsInitialisationTimes_patch(i) = &
      sqrt((sumOfSquaresSavingsInit - sumSavingsInit * sumSavingsInit * inverseNRuns) * inverseNRunsLessOne)

      self % averageSavingsHostTimes_patch(i) = sumSavingsHost * inverseNRuns
      self % stdSavingsHostTimes_patch(i) = &
      sqrt((sumOfSquaresSavingsHost - sumSavingsHost * sumSavingsHost * inverseNRuns) * inverseNRunsLessOne)

      self % averageSavingsTotalTimes_patch(i) = sumSavingsTot * inverseNRuns
      self % stdSavingsTotalTimes_patch(i) = &
      sqrt((sumOfSquaresSavingsTot - sumSavingsTot * sumSavingsTot * inverseNRuns) * inverseNRunsLessOne)

      if(1 < i .and. any(self % depths == 1)) then
        self % averageSavingsInitialisationTimes_patch_1(i - 1) = sumSavingsInit_1 * inverseNRuns
        self % stdSavingsInitialisationTimes_patch_1(i - 1) = &
        sqrt((sumOfSquaresSavingsInit_1 - sumSavingsInit_1 * sumSavingsInit_1 * inverseNRuns) * inverseNRunsLessOne)

      end if

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

    print *, 'Octree: '
    print *, 'Seed: ', self % seeds(runNumber)

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
    print*, "-------------------------------------------------------------"
    print*, "/\/\ Initialisation procedure time /\/\"
    print*, "CPU time: ", self % initialisationTimes_other(runNumber), " seconds"
    print *, "Acceleration structure storage size: ", uMesh % getAccelerationStructureStorageSize() / 1.0e6_defReal, " MBs"
    print*, "-------------------------------------------------------------"

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
    print*, "-------------------------------------------------------------"
    print*, "/\/\ Host element determination time /\/\"
    print*, "CPU time: ", self % hostTimes_other(runNumber), " seconds"
    print*, "-------------------------------------------------------------"

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

    print *, 'Patch-Search: '
    print *, 'Depth: ', self % depths(depthNumber)
    print *, 'Seed: ', self % seeds(runNumber)

    ! Initialise dictionary then fill it.
    call dict % init(2)
    call dict % store('type', 'patchSearchAcceleration')
    call dict % store('depth', self % depths(depthNumber))

    ! Now initialise acceleration structure in the mesh.
    call cpu_time(t1)
    call uMesh % initAccelerationStructure(dict)
    call cpu_time(t2)
    self % initialisationTimes_patch(depthNumber, runNumber) = t2 - t1
    print*, "-------------------------------------------------------------"
    print*, "/\/\ Initialisation procedure time /\/\"
    print*, "CPU time: ", self % initialisationTimes_patch(depthNumber, runNumber), " seconds"
    print *, "Acceleration structure storage size: ", uMesh % getAccelerationStructureStorageSize() / 1.0e6_defReal, " MBs"
    print*, "-------------------------------------------------------------"

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
    print*, "-------------------------------------------------------------"
    print*, "/\/\ Host element determination time /\/\"
    print*, "CPU time: ", self % hostTimes_patch(depthNumber, runNumber), " seconds"
    print*, "-------------------------------------------------------------"

    call uMesh % killAccelerationStructure()

  end subroutine runSingle_patch

end module hostElementDeterminationPackage_class