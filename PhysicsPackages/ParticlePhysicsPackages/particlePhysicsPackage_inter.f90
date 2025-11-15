module particlePhysicsPackage_inter

  use collisionOperator_class,       only : collisionOperator
  use dictionary_class,              only : dictionary
  use errors_mod,                    only : fatalError
  use fieldFactory_func,             only : new_field
  use genericProcedures,             only : numToChar
  use nuclearDatabase_inter,         only : nuclearDatabase
  use nuclearDataReg_mod,            only : activateNuclearDataRegistry => activate, getNuclearDataRegistry => get
  use numPrecision
  use outputFile_class,              only : outputFile
  use particleDungeon_class,         only : particleDungeon
  use physicalParticle_inter,        only : physicalParticle
  use RNG_class,                     only : RNG
  use source_inter,                  only : source
  use sourceFactory_func,            only : new_source
  use tallyAdmin_class,              only : tallyAdmin
  use timer_mod,                     only : timerReset, timerStart, timerStop, timerTime
  use transportOperator_inter,       only : transportOperator
  use transportOperatorFactory_func, only : new_transportOperator
  use physicsPackage_inter,          only : init_super => init, initPhysicsPackagePayload, kill_super => kill, physicsPackage
  use universalVariables,            only : nameWW, P_NEUTRON_CE, P_NEUTRON_MG

  implicit none
  private

  ! Public procedures.
  public :: collectSpecificResults, init, initCycle, kill

  !!
  !!
  !!
  type, public, abstract, extends(physicsPackage) :: particlePhysicsPackage
    private
    integer(shortInt)                             :: bufferSize = 0, particleType = 0
    real(defReal)                                 :: time_transport = ZERO
    logical(defBool)                              :: printSource = .false.
    class(nuclearDatabase), pointer               :: nucData => null()
    class(RNG), pointer                           :: pRNG => null()
    class(source), allocatable                    :: particleSource
    class(transportOperator), allocatable         :: transOp
    type(collisionOperator)                       :: collOp
    type(particleDungeon), pointer                :: currentCycle => null()
  contains
    procedure                                    :: collectSpecificResults
    procedure(displayCycleProgress), deferred    :: displayCycleProgress
    procedure                                    :: generateInitialState
    procedure                                    :: generateSource
    procedure                                    :: getBufferSize
    procedure                                    :: getCollisionOperator
    procedure                                    :: getCurrentCyclePtr
    procedure(getCycleParticlesNumber), deferred :: getCycleParticlesNumber
    procedure                                    :: getInactiveCyclesNumber
    procedure                                    :: getParticleType
    procedure                                    :: getPrintSource
    procedure                                    :: getRNGPtr
    procedure(getTallyAdminPtr), deferred        :: getTallyAdminPtr
    procedure                                    :: getTransportOperator
    procedure                                    :: init
    procedure                                    :: initCycle
    procedure                                    :: kill
    procedure(processEndOfCycle), deferred       :: processEndOfCycle
    procedure                                    :: runCycle
    procedure                                    :: runCycles
    procedure                                    :: setCurrentCyclePtr
    procedure                                    :: setCyclesActive
    procedure(trackParticleHistory), deferred    :: trackParticleHistory
    procedure                                    :: updateNuclearData
  end type particlePhysicsPackage

  abstract interface
    !!
    !!
    !!
    subroutine displayCycleProgress(self, cycleNumber, nInitialParticles, nFinalParticles, elapsedTime, endTime, timeToEnd)
      import                                    :: defReal, particlePhysicsPackage, shortInt
      class(particlePhysicsPackage), intent(in) :: self
      integer(shortInt), intent(in)             :: cycleNumber, nInitialParticles, nFinalParticles
      real(defReal), intent(in)                 :: elapsedTime, endTime, timeToEnd
    end subroutine displayCycleProgress

    !!
    !!
    !!
    function getCycleParticlesNumber(self) result(nParticles)
      import                                    :: particlePhysicsPackage, shortInt
      class(particlePhysicsPackage), intent(in) :: self
      integer(shortInt)                         :: nParticles
    end function getCycleParticlesNumber

    !!
    !!
    !!
    function getTallyAdminPtr(self) result(tallyAdminPtr)
      import                                    :: particlePhysicsPackage, tallyAdmin
      class(particlePhysicsPackage), intent(in) :: self
      type(tallyAdmin), pointer                 :: tallyAdminPtr
    end function getTallyAdminPtr

    !!
    !!
    !!
    subroutine processEndOfCycle(self, nFinalParticles)
      import                                       :: particlePhysicsPackage, shortInt
      class(particlePhysicsPackage), intent(inout) :: self
      integer(shortInt), intent(out)               :: nFinalParticles
    end subroutine processEndOfCycle

    !!
    !!
    !!
    subroutine trackParticleHistory(self, transOp, collOp, p, buffer, tally)
      import :: collisionOperator, particleDungeon, particlePhysicsPackage, physicalParticle, tallyAdmin, transportOperator
      class(particlePhysicsPackage), intent(in) :: self
      class(transportOperator), intent(inout)   :: transOp
      type(collisionOperator), intent(inout)    :: collOp
      class(physicalParticle), intent(inout)    :: p
      type(particleDungeon), intent(inout)      :: buffer
      type(tallyAdmin), intent(inout)           :: tally
    end subroutine trackParticleHistory

  end interface

  !!
  !!
  !!
  type, public, extends(initPhysicsPackagePayload) :: initParticlePhysicsPackagePayload
    integer(shortInt) :: defaultBufferSize = 0, currentCycleSizeMultiplier = 1
    logical(defBool)  :: isSourceRequired = .false.
  end type initParticlePhysicsPackagePayload

contains
  !!
  !!
  !!
  subroutine collectSpecificResults(self, out)
    class(particlePhysicsPackage), intent(in) :: self
    type(outputFile), intent(inout)           :: out
    character(nameLen)                        :: name

    name = 'Transport_time'
    call out % printValue(self % time_transport, name)

  end subroutine collectSpecificResults

  !!
  !!
  !!
  subroutine generateInitialState(self)
    class(particlePhysicsPackage), intent(inout) :: self

    ! Do nothing.

  end subroutine generateInitialState

  !!
  !!
  !!
  subroutine generateSource(self)
    class(particlePhysicsPackage), intent(inout) :: self

    call self % particleSource % generate(self % currentCycle, self % getParticlesNumber(), self % pRNG)

  end subroutine generateSource

  !!
  !!
  !!
  elemental function getBufferSize(self) result(bufferSize)
    class(particlePhysicsPackage), intent(in) :: self
    integer(shortInt)                         :: bufferSize

    bufferSize = self % bufferSize

  end function getBufferSize

  !!
  !!
  !!
  function getCollisionOperator(self) result(collOp)
    class(particlePhysicsPackage), intent(in) :: self
    type(collisionOperator)                   :: collOp

    collOp = self % collOp

  end function getCollisionOperator

  !!
  !!
  !!
  function getCurrentCyclePtr(self) result(currentCyclePtr)
    class(particlePhysicsPackage), intent(in) :: self
    type(particleDungeon), pointer            :: currentCyclePtr

    currentCyclePtr => self % currentCycle

  end function getCurrentCyclePtr

  !!
  !!
  !!
  elemental function getInactiveCyclesNumber(self) result(nInactiveCycles)
    class(particlePhysicsPackage), intent(in) :: self
    integer(shortInt)                         :: nInactiveCycles

    nInactiveCycles = 0

  end function getInactiveCyclesNumber

  !!
  !!
  !!
  elemental function getParticleType(self) result(particleType)
    class(particlePhysicsPackage), intent(in) :: self
    integer(shortInt)                         :: particleType

    particleType = self % particleType

  end function getParticleType

  !!
  !!
  !!
  elemental function getPrintSource(self) result(printSource)
    class(particlePhysicsPackage), intent(in) :: self
    logical(defBool)                          :: printSource

    printSource = self % printSource

  end function getPrintSource

  !!
  !!
  !!
  function getRNGPtr(self) result(pRNGPtr)
    class(particlePhysicsPackage), intent(in) :: self
    type(RNG), pointer                        :: pRNGPtr

    pRNGPtr => self % pRNG

  end function getRNGPtr

  !!
  !!
  !!
  function getTransportOperator(self) result(transOp)
    class(particlePhysicsPackage), intent(in) :: self
    class(transportOperator), allocatable     :: transOp

    allocate(transOp, source = self % transOp)

  end function getTransportOperator

  !!
  !!
  !!
  subroutine init(self, payload)
    class(particlePhysicsPackage), intent(inout)     :: self
    class(initPhysicsPackagePayload), intent(in)     :: payload
    type(initParticlePhysicsPackagePayload), pointer :: payloadPtr
    character(nameLen)                               :: energy, nucData
    type(dictionary)                                 :: sourceDict
    character(*), parameter                          :: here = 'init (particlePhysicsPackage_inter.f90)'

    ! Downcast payload to correct type.
    select type(ptr => payload)
      type is(initParticlePhysicsPackagePayload)
        payloadPtr => ptr

      class default
        call fatalError(here, 'Invalid payload type.')

    end select

    ! Initialise superclass.
    call init_super(self, payloadPtr)

    ! Load energy from dictionary.
    call payloadPtr % dict % get(energy, 'dataType')
    
    ! Process type of data.
    select case(energy)
      case('mg')
        self % particleType = P_NEUTRON_MG

      case('ce')
        self % particleType = P_NEUTRON_CE

      case default
        call fatalError(here, "dataType must be 'mg' or 'ce'.")

    end select

    ! Load nuclear data, parallel buffer size, and whether to print particle source per cycle from dictionary.
    call payloadPtr % dict % get(nucData, 'XSdata')
    call payloadPtr % dict % getOrDefault(self % bufferSize, 'buffer', payloadPtr % defaultBufferSize)
    call payloadPtr % dict % getOrDefault(self % printSource, 'printSource', .false.)

    ! Initialise RNG.
    allocate(self % pRNG)
    call self % pRNG % init(self % getInitialSeed())

    ! Activate Nuclear Data. Note: all materials are active.
    call activateNuclearDataRegistry(self % particleType, nucData, payloadPtr % geometry % activeMats())
    self % nucData => getNuclearDataRegistry(self % particleType)

    ! Call visualisation.
    if (payloadPtr % dict % isPresent('viz')) call self % buildVisualisation(payloadPtr % dict % getDictPtr('viz'))

    ! Build collision operator.
    call self % collOp % init(payloadPtr % dict % getDictPtr('collisionOperator'))

    ! Build transport operator.
    call new_transportOperator(self % transOp, payloadPtr % dict % getDictPtr('transportOperator'))

    ! Read variance reduction option as a geometry field.
    if (payloadPtr % dict % isPresent('varianceReduction')) &
    call new_field(payloadPtr % dict % getDictPtr('varianceReduction'), nameWW)

    ! Read source.
    if (payloadPtr % dict % isPresent('source')) then
      call new_source(self % particleSource, payloadPtr % dict % getDictPtr('source'), payloadPtr % geometry)

    else
      if (payloadPtr % isSourceRequired) call fatalError(here, 'Missing "source" dictionary.')
      ! Build source.
      call sourceDict % init(1)
      call sourceDict % store('type', merge('CEFissionSource', 'MGFissionSource', self % particleType == P_NEUTRON_CE))
      call new_source(self % particleSource, sourceDict, payloadPtr % geometry)
      call sourceDict % kill()

    end if

    ! Allocate currentCycle and initialise it to correct size.
    allocate(self % currentCycle)
    call self % currentCycle % init(payloadPtr % currentCycleSizeMultiplier * self % getParticlesNumber())

  end subroutine init

  !!
  !!
  !!
  subroutine initCycle(self)
    class(particlePhysicsPackage), intent(inout) :: self

    ! Do nothing.

  end subroutine initCycle

  !!
  !!
  !!
  subroutine kill(self)
    class(particlePhysicsPackage), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % bufferSize = 0
    self % particleType = 0
    self % time_transport = ZERO
    self % printSource = .false.
    self % nucData => null()
    if (associated(self % pRNG)) deallocate(self % pRNG)
    if (allocated(self % particleSource)) then
      call self % particleSource % kill()
      deallocate(self % particleSource)

    end if
    if (allocated(self % transOp)) then
      call self % transOp % kill()
      deallocate(self % transOp)

    end if
    call self % collOp % kill()
    if (associated(self % currentCycle)) then
      call self % currentCycle % kill()
      deallocate(self % currentCycle)

    end if

  end subroutine kill

  !!
  !!
  !!
  subroutine runCycle(self, cycleNumber, nCycles, p, transOp, geometryIdx, nInitialParticles, collOp, buffer, pRNG, tally)
    class(particlePhysicsPackage), intent(inout)        :: self
    integer(shortInt), intent(in)                       :: cycleNumber, nCycles
    class(physicalParticle), allocatable, intent(inout) :: p
    class(transportOperator), intent(inout)             :: transOp
    integer(shortInt), intent(inout)                    :: geometryIdx, nInitialParticles
    type(collisionOperator), intent(inout)              :: collOp
    type(particleDungeon), intent(inout)                :: buffer
    type(RNG), intent(inout)                            :: pRNG
    type(tallyAdmin), pointer, intent(inout)            :: tally
    integer(shortInt)                                   :: i, nFinalParticles, timerMain
    real(defReal)                                       :: elapsedTime, endTime

    !$omp master
    ! Prepare current cycle.
    call self % initCycle()
    if (self % printSource) call self % currentCycle % printToFile(trim(self % getOutputFile())//'_source'//numToChar(cycleNumber))
    nInitialParticles = self % getCycleParticlesNumber()
    call tally % reportCycleStart(self % currentCycle)
    !$omp end master

    ! Wait for master thread before launching parallel execution.
    !$omp barrier

    geometryIdx = self % getGeometryIdx()
  
    !$omp do schedule(dynamic)
    do i = 1, nInitialParticles
      ! Create RNG which can be thread private.
      pRNG = self % pRNG

      ! Generate a particle from the dungeon, prepate it and track its history.
      p = self % currentCycle % copy(i)
      call p % setGeometryIdx(geometryIdx)
      call p % setRNGPtr(pRNG)
      call p % strideRNG(i)

      call self % trackParticleHistory(transOp, collOp, p, buffer, tally)

    end do
    !$omp end do

    !$omp master
    ! Process end of cycle results.
    call self % processEndOfCycle(nFinalParticles)
    
    ! Stop timer and display progress so far.
    timerMain = self % getTimerMain()
    call timerStop(timerMain)
    elapsedTime = timerTime(timerMain)
    self % time_transport = elapsedTime
    endTime = self % getTotalCyclesNumber() * elapsedTime / self % getCurrentCycleNumber(cycleNumber)
    call self % displayCycleProgress(cycleNumber, nInitialParticles, nFinalParticles, elapsedTime, endTime, &
                                     max(ZERO, endTime - elapsedTime))
    call tally % display()
    !$omp end master

  end subroutine runCycle

  !!
  !!
  !!
  subroutine runCycles(self, nCycles, reset)
    class(particlePhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                :: nCycles
    logical(defBool), intent(in), optional       :: reset
    class(physicalParticle), allocatable         :: p
    class(transportOperator), allocatable        :: transOp
    integer(shortInt)                            :: geometryIdx, i, nInitialParticles, timerMain
    logical(defBool)                             :: resetTimer
    type(tallyAdmin), pointer                    :: tallyAdminPtr
    type(particleDungeon)                        :: buffer
    type(collisionOperator)                      :: collOp
    type(RNG)                                    :: pRNG

    ! Reset and start timer.
    resetTimer = .true.
    if (present(reset)) resetTimer = reset

    if (resetTimer) then
      timerMain = self % getTimerMain()
      call timerReset(timerMain)
      call timerStart(timerMain)

    end if

    geometryIdx = 0
    nInitialParticles = 0

    ! Create parallel region once outside the main loop for performance.
    !$omp parallel private(buffer, collOp, i, p, pRNG, transOp) &
    !$omp shared(geometryIdx, nCycles, nInitialParticles, self, tallyAdminPtr)

    ! Create particle buffer and a transport operator which can be made thread private
    call buffer % init(self % bufferSize)
    collOp = self % collOp
    allocate(transOp, source = self % transOp)

    !$omp master
    tallyAdminPtr => self % getTallyAdminPtr()
    !$omp end master
    !$omp barrier

    ! Loop through all cycles.
    do i = 1, nCycles
      call self % runCycle(i, nCycles, p, transOp, geometryIdx, nInitialParticles, collOp, buffer, pRNG, tallyAdminPtr)

    end do
    !$omp end parallel

  end subroutine runCycles

  !!
  !!
  !!
  subroutine setCurrentCyclePtr(self, currentCyclePtr)
    class(particlePhysicsPackage), intent(inout) :: self
    type(particleDungeon), pointer, intent(in)   :: currentCyclePtr

    self % currentCycle => currentCyclePtr

  end subroutine setCurrentCyclePtr

  !!
  !!
  !!
  elemental subroutine setCyclesActive(self)
    class(particlePhysicsPackage), intent(inout) :: self

    ! Do nothing by default.

  end subroutine setCyclesActive

  !!
  !!
  !!
  subroutine updateNuclearData(self)
    class(particlePhysicsPackage), intent(inout) :: self
    character(*), parameter                      :: here = 'updateNuclearData (particlePhysicsPackage_inter.f90)'

    if (.not. associated(self % nucData)) call fatalError(here, 'Attempting to update unassociated nuclear data.')
    call self % nucData % updateMaterialsProperties()

  end subroutine updateNuclearData

end module particlePhysicsPackage_inter