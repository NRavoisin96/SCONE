module NTHPackage_class

  use collisionOperator_class,          only : collisionOperator
  use dictionary_class,                 only : dictionary
  use eigenPhysicsPackage_class,        only : eigenPhysicsPackage
  use errors_mod,                       only : fatalError
  use fixedSourcePhysicsPackage_class,  only : fixedSourcePhysicsPackage
  use genericProcedures,                only : numToChar, printFishLineR
  use geometryReg_mod,                  only : fieldPtrByName
  use heatTransferPhysicsPackage_class, only : heatTransferPhysicsPackage
  use numPrecision
  use outputFile_class,                 only : outputFile
  use particleDungeon_class,            only : particleDungeon
  use particlePhysicsPackage_inter,     only : initParticlePhysicsPackagePayload, particlePhysicsPackage
  use physicalParticle_inter,           only : physicalParticle
  use physicsPackage_inter,             only : copyPayload, init_super => init, initPhysicsPackagePayload, physicsPackage
  use RNG_class,                        only : RNG
  use scalarField_inter,                only : castScalarFieldPtr, scalarField
  use tallyAdmin_class,                 only : tallyAdmin
  use tallyResult_class,                only : castTallyResultArraysPtr, tallyResult, tallyResultArrays
  use timer_mod,                        only : secToChar, timerReset, timerStart, timerStop, timerTime
  use transportOperator_inter,          only : transportOperator
  use universalVariables,               only : nameDensity, nameHeatSource, nameTemperature

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(physicsPackage) :: NTHPackage
    private
    class(particlePhysicsPackage), allocatable :: neutronicsPackage
    real(defReal), dimension(:), allocatable   :: densityMeans, densityVariances
    type(heatTransferPhysicsPackage)           :: heatTransferPackage
  contains
    procedure          :: collectSpecificResults
    procedure          :: displayCycleProgress
    procedure          :: init
    procedure          :: kill
    procedure          :: run
    procedure, private :: runCycles
  end type NTHPackage

contains
  !!
  !!
  !!
  subroutine collectSpecificResults(self, out)
    class(NTHPackage), intent(in)   :: self
    type(outputFile), intent(inout) :: out

  end subroutine collectSpecificResults

  !!
  !!
  !!
  subroutine displayCycleProgress(self, cycleNumber, nInitialParticles, nFinalParticles, elapsedTime, endTime, timeToEnd)
    class(NTHPackage), intent(in) :: self
    integer(shortInt), intent(in) :: cycleNumber, nInitialParticles, nFinalParticles
    real(defReal), intent(in)     :: elapsedTime, endTime, timeToEnd

    ! Display progress
    call printFishLineR(cycleNumber)
    print *
    if(.not. self % neutronicsPackage % getCyclesActive()) then
      print *, 'Cycle: ', numToChar(cycleNumber), ' of ', numToChar(0)

    else
      print *, 'Cycle: ', numToChar(cycleNumber), ' of ', numToChar(self % neutronicsPackage % getCyclesNumber())

    end if
    print *, 'Pop: ', numToChar(nInitialParticles) , ' -> ', numToChar(nFinalParticles)
    print *, 'Elapsed time: ', trim(secToChar(elapsedTime))
    print *, 'End time:     ', trim(secToChar(endTime))
    print *, 'Time to end:  ', trim(secToChar(timeToEnd))

  end subroutine displayCycleProgress

  !!
  !!
  !!
  subroutine init(self, payload)
    class(NTHPackage), intent(inout)             :: self
    class(initPhysicsPackagePayload), intent(in) :: payload
    character(nameLen)                           :: type
    class(dictionary), pointer                   :: activeTalliesDict, clerksDict, currentPackageDict, &
                                                    inactiveTalliesDict, packagesDict
    class(scalarField), pointer                  :: heatSourceFieldPtr, temperatureFieldPtr
    type(initParticlePhysicsPackagePayload)      :: neutronicPackagePayload
    type(initPhysicsPackagePayload)              :: heatTransgerPackagePayload
    character(*), parameter                      :: here = 'init (NTHPackage_class.f90)'

    ! Initialise superclass.
    call init_super(self, payload)

    ! Check that there are fields for temperature and fission power.
    heatSourceFieldPtr => castScalarFieldPtr(fieldPtrByName(nameHeatSource))
    temperatureFieldPtr => castScalarFieldPtr(fieldPtrByName(nameTemperature))

    ! Get dictionary to physics package definitions.
    if (.not. payload % dict % isPresent('packages')) call fatalError(here, 'Missing "packages" subdictionary.')
    packagesDict => payload % dict % getDictPtr('packages')

    ! Verify inputs for neutron tranport package.
    if (.not. packagesDict % isPresent('neutronics')) call fatalError(here, 'Missing "neutronics" package subdictionary.')
    currentPackageDict => packagesDict % getDictPtr('neutronics')

    ! Check that there are inactive and active tallies.
    if (.not. currentPackageDict % isPresent('inactiveTally')) &
    call fatalError(here, 'Missing "inactiveTally" subdictionary for neutronics package.')

    if (.not. currentPackageDict % isPresent('activeTally')) &
    call fatalError(here, 'Missing "activeTally" subdictionary for neutronics package.')

    ! Get pointers to inactive and active tally subdictionaries.
    inactiveTalliesDict => currentPackageDict % getDictPtr('inactiveTally')
    activeTalliesDict => currentPackageDict % getDictPtr('activeTally')

    ! Check that there are clerks defined for inactive and active tally subdictionaries.
    if (.not. inactiveTalliesDict % isPresent('clerks')) &
    call fatalError(here, 'Missing "clerks" subdictionary in inactive tallies for neutronics package.')
    if (.not. activeTalliesDict % isPresent('clerks')) &
    call fatalError(here, 'Missing "clerks" subdictionary in active tallies for neutronics package.')

    ! Check that there is a 'fissionPower' entry in both inactive and active tallies.
    clerksDict => inactiveTalliesDict % getDictPtr('clerks')
    if (.not. clerksDict % isPresent('fissionPower')) &
    call fatalError(here, 'Missing "fissionPower" clerk definition in inactive tallies for neutronics package.')

    clerksDict => activeTalliesDict % getDictPtr('clerks')
    if (.not. clerksDict % isPresent('fissionPower')) &
    call fatalError(here, 'Missing "fissionPower" clerk definition in active tallies for neutronics package.')

    ! Verify inputs for heat transfer package.
    if (.not. packagesDict % isPresent('heatTransfer')) call fatalError(here, 'Missing "heatTransfer" package subdictionary.')
    currentPackageDict => packagesDict % getDictPtr('heatTransfer')

    ! Now initialise both physics package.
    currentPackageDict => packagesDict % getDictPtr('neutronics')
    call currentPackageDict % get(type, 'type')
    select case(type)
      case('eigenPhysicsPackage')
        allocate(eigenPhysicsPackage :: self % neutronicsPackage)

      case('fixedSourcePhysicsPackage')
        allocate(fixedSourcePhysicsPackage :: self % neutronicsPackage)

      case default
        call fatalError(here, 'Invalid neutronics physics package.')

    end select
    call copyPayload(payload, neutronicPackagePayload)
    neutronicPackagePayload % dict => currentPackageDict
    call self % neutronicsPackage % init(neutronicPackagePayload)
    
    call copyPayload(payload, heatTransgerPackagePayload)
    heatTransgerPackagePayload % dict => packagesDict % getDictPtr('heatTransfer')
    call self % heatTransferPackage % init(heatTransgerPackagePayload)

    ! Allocate memory (hardcode here.)
    allocate(self % densityMeans(10), self % densityVariances(10))

  end subroutine init

  !!
  !!
  !!
  subroutine kill(self)
    class(NTHPackage), intent(inout) :: self

    if (allocated(self % neutronicsPackage)) then
      call self % neutronicsPackage % kill()
      deallocate(self % neutronicsPackage)

    end if
    if(allocated(self % densityMeans)) deallocate(self % densityMeans)
    if(allocated(self % densityVariances)) deallocate(self % densityVariances)
    call self % heatTransferPackage % kill()

  end subroutine kill

  !!
  !!
  !!
  subroutine run(self)
    class(NTHPackage), intent(inout) :: self
    class(scalarField), pointer      :: densityFieldPtr, heatSourceFieldPtr, temperatureFieldPtr
    integer(shortInt)                :: nInactiveCycles, timerMain
    logical(defBool)                 :: hasTimerStarted     
    type(tallyAdmin), pointer        :: tallyAdminPtr
    character(*), parameter          :: HERE = 'run (NTHPackage_class.f90)'

    ! Initialise hasTimerStarted
    hasTimerStarted = .false.

    ! Generate initial state for neutronics package.
    if (.not. allocated(self % neutronicsPackage)) call fatalError(HERE, 'Neutronics physics package is not allocated.')
    call self % neutronicsPackage % generateInitialState()

    ! Get pointers to heat source and temperature fields.
    densityFieldPtr => castScalarFieldPtr(fieldPtrByName(nameDensity))
    heatSourceFieldPtr => castScalarFieldPtr(fieldPtrByName(nameHeatSource))
    temperatureFieldPtr => castScalarFieldPtr(fieldPtrByName(nameTemperature))

    timerMain = self % getTimerMain()

    ! Get number of inactive cycles.
    nInactiveCycles = self % neutronicsPackage % getInactiveCyclesNumber()
    if (0 < nInactiveCycles) then
      ! Reset and start timer.
      call timerReset(timerMain)
      call timerStart(timerMain)
      hasTimerStarted = .true.
      
      ! Get pointer to tallyAdmin.
      tallyAdminPtr => self % neutronicsPackage % getTallyAdminPtr()
      call self % runCycles(nInactiveCycles, densityFieldPtr, heatSourceFieldPtr, temperatureFieldPtr, tallyAdminPtr, .true.)
      call self % neutronicsPackage % setCyclesActive()

    end if

    ! Run active cycles.
    if(.not. hasTimerStarted) then
      call timerReset(timerMain)
      call timerStart(timerMain)

    end if
    tallyAdminPtr => self % neutronicsPackage % getTallyAdminPtr()
    call self % runCycles(self % neutronicsPackage % getCyclesNumber(), densityFieldPtr, heatSourceFieldPtr, &
                          temperatureFieldPtr, tallyAdminPtr)

  end subroutine run

  !!
  !!
  !!
  subroutine runCycles(self, nCycles, densityFieldPtr, heatSourceFieldPtr, temperatureFieldPtr, tallyAdminPtr, flush)
    class(NTHPackage), intent(inout)           :: self
    integer(shortInt), intent(in)              :: nCycles
    class(scalarField), pointer, intent(inout) :: densityFieldPtr, heatSourceFieldPtr, temperatureFieldPtr
    type(tallyAdmin), pointer, intent(inout)   :: tallyAdminPtr
    logical(defBool), intent(in), optional     :: flush
    class(physicalParticle), allocatable       :: p
    class(tallyResult), allocatable            :: tallyResults
    class(transportOperator), allocatable      :: transOp
    integer(shortInt)                          :: i, j, k, geometryIdx, nInitialParticles, timerMain
    logical(defBool)                           :: flushResults
    real(defReal)                              :: endTime, elapsedTime
    real(defReal), dimension(10)               :: densities, sum, sumOfSquares
    real(defReal), dimension(:), allocatable   :: meanTemperatures
    type(collisionOperator)                    :: collOp
    type(particleDungeon)                      :: buffer
    type(RNG)                                  :: pRNG
    type(tallyResultArrays), pointer           :: tallyResultArraysPtr
    character(*), parameter                    :: HERE = 'runCycles (NTHPackage_class.f90)'

    flushResults = .false.
    if (present(flush)) flushResults = flush

    ! Initialise shared variables.
    geometryIdx = 0
    nInitialParticles = 0
    sum = ZERO
    sumOfSquares = ZERO

    ! Create parallel region here.
    !$omp parallel default(shared) private(buffer, collOp, i, j, p, pRNG, tallyResultArraysPtr, transOp)

    ! Initialise buffer then get collision and transport operators from neutronics package.
    call buffer % init(self % neutronicsPackage % getBufferSize())
    collOp = self % neutronicsPackage % getCollisionOperator()
    allocate(transOp, source = self % neutronicsPackage % getTransportOperator())
    !$omp barrier

    ! Loop through inactive cycles first (if any).
    do i = 1, nCycles
      ! Run neutronics simulation and get fission power results from tallyAdminPtr.
      call self % neutronicsPackage % runCycle(i, nCycles, p, transOp, geometryIdx, nInitialParticles, collOp, buffer, &
                                               pRNG, tallyAdminPtr, .false.)
      !$omp barrier
      
      !$omp master
      call tallyAdminPtr % getResult(tallyResults, 'fissionPower')

      ! Downcast tallyResults to correct type and update heat source field.
      tallyResultArraysPtr => castTallyResultArraysPtr(tallyResults)
      if (.not. allocated(tallyResultArraysPtr % results)) call fatalError(HERE, 'Empty tally results.')
      do j = 1, size(tallyResultArraysPtr % results)
        if (tallyResultArraysPtr % results(j) % clerkName == 'fissionPower') then
          call heatSourceFieldPtr % setValues(tallyResultArraysPtr % results(j) % values)

          do k = 1, 10
            print *, 'Fission power in element '//numToChar(k)//': ', tallyResultArraysPtr % results(j) % values(k), '+/-', &
            tallyResultArraysPtr % results(j) % standardDeviations(k)

          end do

        end if

      end do
      !$omp end master
      !$omp barrier

      ! Run heat transfer simulation and update temperature field. This is ugly for now (NR).
      call self % heatTransferPackage % runWalkers()
      !$omp barrier
      
      !$omp master
      meanTemperatures = self % heatTransferPackage % getMeans()
      call temperatureFieldPtr % setValues(meanTemperatures)
      densities = 1.933346e4_defReal - 7.9647e-1_defReal * meanTemperatures
      call densityFieldPtr % setValues(densities)

      ! Update neutronics package nuclear data using the new temperature field and flush tallies.
      call self % neutronicsPackage % updateNuclearData()
      if (flushResults) then
        call self % heatTransferPackage % flushResults()
        call tallyAdminPtr % flush('fissionPower')

      else
        sum = sum + densities
        sumOfSquares = sumOfSquares + densities * densities

      end if

      ! Display progress so far.
      timerMain = self % getTimerMain()
      call timerStop(timerMain)
      elapsedTime = timerTime(timerMain)
      endTime = self % neutronicsPackage % getTotalCyclesNumber() * elapsedTime / &
                self % neutronicsPackage % getCurrentCycleNumber(i)
      call self % displayCycleProgress(i, nInitialParticles, 0, elapsedTime, endTime, &
                                       max(ZERO, endTime - elapsedTime))
      call tallyAdminPtr % display()
      !$omp end master
      !$omp barrier

    end do

    !$omp end parallel

    if (.not. flushResults) then
      ! Normalise densities.
      do i = 1, 10
        if(nCycles == 1) then
          self % densityMeans(i) = sum(i)
          self % densityVariances(i) = ZERO

        else
          self % densityMeans(i) = sum(i) / nCycles
          self % densityVariances(i) = (sumOfSquares(i) - sum(i) * sum(i) / nCycles) / (nCycles - 1)

        end if
        print *, 'Density of element '//numToChar(i)//': ', self % densityMeans(i), '+/-', sqrt(self % densityVariances(i))

      end do

    end if

  end subroutine runCycles

end module NTHPackage_class