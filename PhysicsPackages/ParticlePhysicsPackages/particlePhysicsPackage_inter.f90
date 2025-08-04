module particlePhysicsPackage_inter

  use collisionOperator_class,       only : collisionOperator
  use dictionary_class,              only : dictionary
  use fieldFactory_func,             only : new_field
  use genericProcedures,             only : fatalError, numToChar
  use nuclearDatabase_inter,         only : nuclearDatabase
  use nuclearDataReg_mod,            only : activateNuclearDataRegistry => activate, getNuclearDataRegistry => get
  use numPrecision
  use outputFile_class,              only : outputFile
  use particle_class,                only : particle
  use particleDungeon_class,         only : particleDungeon
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
    procedure                                    :: generateSource
    procedure                                    :: getBufferSize
    procedure                                    :: getCurrentCyclePtr
    procedure(getCycleParticlesNumber), deferred :: getCycleParticlesNumber
    procedure                                    :: getParticleType
    procedure                                    :: getPrintSource
    procedure                                    :: getRNGPtr
    procedure(getTallyAdminPtr), deferred        :: getTallyAdminPtr
    procedure                                    :: init
    procedure                                    :: initCycle
    procedure                                    :: kill
    procedure(processEndOfCycle), deferred       :: processEndOfCycle
    procedure                                    :: runCycles
    procedure                                    :: setCurrentCyclePtr
    procedure(trackParticleHistory), deferred    :: trackParticleHistory
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
      import :: collisionOperator, particle, particleDungeon, particlePhysicsPackage, tallyAdmin, transportOperator
      class(particlePhysicsPackage), intent(in) :: self
      class(transportOperator), intent(inout)   :: transOp
      type(collisionOperator), intent(inout)    :: collOp
      type(particle), intent(inout)             :: p
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
  function getCurrentCyclePtr(self) result(currentCyclePtr)
    class(particlePhysicsPackage), intent(in) :: self
    type(particleDungeon), pointer            :: currentCyclePtr

    currentCyclePtr => self % currentCycle

  end function getCurrentCyclePtr

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
      call sourceDict % init(3)
      call sourceDict % store('type', 'fissionSource')
      call sourceDict % store('data', trim(energy))
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
  subroutine runCycles(self, nCycles)
    class(particlePhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                :: nCycles
    integer(shortInt)                            :: i, j, nFinalParticles, nInitialParticles, timerMain
    type(tallyAdmin), pointer                    :: tallyAdminPtr
    type(particleDungeon)                        :: buffer
    type(particle)                               :: p
    type(collisionOperator)                      :: collOp
    class(transportOperator), allocatable        :: transOp
    type(RNG), target                            :: pRNG
    real(defReal)                                :: elapsedTime, endTime

    ! Reset and start timer.
    timerMain = self % getTimerMain()
    call timerReset(timerMain)
    call timerStart(timerMain)

    ! Create parallel region once outside the main loop for performance.
    !$omp parallel private(p, buffer, pRNG, collOp, transOp, j) &
    !$omp shared(nFinalParticles, nInitialParticles, tallyAdminPtr)

    ! Create particle buffer and a transport operator which can be made thread private
    call buffer % init(self % bufferSize)
    collOp = self % collOp
    allocate(transOp, source = self % transOp)

    ! Loop through all cycles.
    do i = 1, nCycles
      !$omp master
      ! Prepare current cycle.
      call self % initCycle()
      if (self % printSource) call self % currentCycle % printToFile(trim(self % getOutputFile())//'_source'//numToChar(i))
      nInitialParticles = self % getCycleParticlesNumber()
      tallyAdminPtr => self % getTallyAdminPtr()
      call tallyAdminPtr % reportCycleStart(self % currentCycle)
      !$omp end master

      ! Wait for master thread before launching parallel execution.
      !$omp barrier

      ! Initialise particle.
      p % geomIdx = self % getGeometryIdx()
    
      !$omp do schedule(dynamic)
      do j = 1, nInitialParticles
        ! Create RNG which can be thread private
        pRNG = self % pRNG
        p % pRNG => pRNG
        call p % pRNG % stride(j)

        ! Obtain particle current cycle dungeon and prepare particle.
        call self % currentCycle % copy(p, j)
        call self % trackParticleHistory(transOp, collOp, p, buffer, tallyAdminPtr)

      end do
      !$omp end do

      !$omp master
      ! Process end of cycle results.
      call self % processEndOfCycle(nFinalParticles)

      ! Stop timer and display progress so far.
      call timerStop(timerMain)
      elapsedTime = timerTime(timerMain)
      endTime = nCycles * elapsedTime / i
      call self % displayCycleProgress(i, nInitialParticles, nFinalParticles, elapsedTime, endTime, &
                                       max(ZERO, endTime - elapsedTime))
      call tallyAdminPtr % display()
      self % time_transport = self % time_transport + elapsedTime
      !$omp end master

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

end module particlePhysicsPackage_inter