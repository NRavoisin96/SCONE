module NTHPackage_class

  use collisionOperator_class,          only : collisionOperator
  use dictionary_class,                 only : dictionary
  use eigenPhysicsPackage_class,        only : eigenPhysicsPackage
  use errors_mod,                       only : fatalError
  use fixedSourcePhysicsPackage_class,  only : fixedSourcePhysicsPackage
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
  use timer_mod,                        only : timerReset, timerStart
  use transportOperator_inter,          only : transportOperator
  use universalVariables,               only : nameHeatSource, nameTemperature

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(physicsPackage) :: NTHPackage
    private
    class(particlePhysicsPackage), allocatable :: neutronicsPackage
    type(heatTransferPhysicsPackage)           :: heatTransferPackage
  contains
    procedure          :: collectSpecificResults
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
    call self % heatTransferPackage % kill()

  end subroutine kill

  !!
  !!
  !!
  subroutine run(self)
    class(NTHPackage), intent(inout)      :: self
    class(scalarField), pointer           :: heatSourceFieldPtr, temperatureFieldPtr
    integer(shortInt)                     :: nInactiveCycles, timerMain
    type(tallyAdmin), pointer             :: tallyAdminPtr
    character(*), parameter               :: HERE = 'run (NTHPackage_class.f90)'

    ! Generate initial state for neutronics package.
    if (.not. allocated(self % neutronicsPackage)) call fatalError(HERE, 'Neutronics physics package is not allocated.')
    call self % neutronicsPackage % generateInitialState()

    ! Get pointers to heat source and temperature fields.
    heatSourceFieldPtr => castScalarFieldPtr(fieldPtrByName(nameHeatSource))
    temperatureFieldPtr => castScalarFieldPtr(fieldPtrByName(nameTemperature))

    timerMain = self % getTimerMain()

    ! Get number of inactive cycles.
    nInactiveCycles = self % neutronicsPackage % getInactiveCyclesNumber()
    if (0 < nInactiveCycles) then
      ! Reset and start timer.
      call timerReset(timerMain)
      call timerStart(timerMain)
      
      ! Get pointer to tallyAdmin.
      tallyAdminPtr => self % neutronicsPackage % getTallyAdminPtr()
      call self % runCycles(nInactiveCycles, heatSourceFieldPtr, temperatureFieldPtr, tallyAdminPtr, .true.)
      call self % neutronicsPackage % setCyclesActive()

    end if

    ! Run active cycles.
    call timerReset(timerMain)
    call timerStart(timerMain)
    tallyAdminPtr => self % neutronicsPackage % getTallyAdminPtr()
    call self % runCycles(self % neutronicsPackage % getCyclesNumber(), heatSourceFieldPtr, temperatureFieldPtr, tallyAdminPtr)

  end subroutine run

  !!
  !!
  !!
  subroutine runCycles(self, nCycles, heatSourceFieldPtr, temperatureFieldPtr, tallyAdminPtr, flush)
    class(NTHPackage), intent(inout)           :: self
    integer(shortInt), intent(in)              :: nCycles
    class(scalarField), pointer, intent(inout) :: heatSourceFieldPtr, temperatureFieldPtr
    type(tallyAdmin), pointer, intent(inout)   :: tallyAdminPtr
    logical(defBool), intent(in), optional     :: flush
    class(physicalParticle), allocatable       :: p
    class(tallyResult), allocatable            :: tallyResults
    class(transportOperator), allocatable      :: transOp
    integer(shortInt)                          :: i, j, geometryIdx, nInitialParticles
    logical(defBool)                           :: flushResults
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
                                               pRNG, tallyAdminPtr)
      !$omp barrier
      
      !$omp master
      call tallyAdminPtr % getResult(tallyResults, 'fissionPower')

      ! Downcast tallyResults to correct type and update heat source field.
      tallyResultArraysPtr => castTallyResultArraysPtr(tallyResults)
      if (.not. allocated(tallyResultArraysPtr % results)) call fatalError(HERE, 'Empty tally results.')
      do j = 1, size(tallyResultArraysPtr % results)
        if (tallyResultArraysPtr % results(j) % clerkName == 'fissionPower') then
          call heatSourceFieldPtr % setValues(tallyResultArraysPtr % results(j) % values)

        end if

      end do
      !$omp end master
      !$omp barrier

      ! Run heat transfer simulation and update temperature field. This is ugly for now (NR).
      call self % heatTransferPackage % runWalkers()
      !$omp barrier
      
      !$omp master
      call temperatureFieldPtr % setValues(self % heatTransferPackage % getMeans())

      ! Update neutronics package nuclear data using the new temperature field and flush tallies.
      call self % neutronicsPackage % updateNuclearData()
      if (flushResults) then
        call self % heatTransferPackage % flushResults()
        call tallyAdminPtr % flush('fissionPower')

      end if
      !$omp end master
      !$omp barrier

    end do

    !$omp end parallel

  end subroutine runCycles

end module NTHPackage_class