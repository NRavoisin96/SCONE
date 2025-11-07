module fixedSourcePhysicsPackage_class

  use collisionOperator_class,      only : collisionOperator
  use dictionary_class,             only : dictionary
  use errors_mod,                   only : fatalError
  use genericProcedures,            only : numToChar, printFishLineR
  use numPrecision
  use outputFile_class,             only : outputFile
  use particleDungeon_class,        only : particleDungeon
  use particlePhysicsPackage_inter, only : collectSpecificResults_super => collectSpecificResults, &
                                           init_super => init, initParticlePhysicsPackagePayload, &
                                           particlePhysicsPackage, kill_super => kill
  use physicalParticle_inter,       only : physicalParticle
  use physicsPackage_inter,         only : copyPayload, initPhysicsPackagePayload
  use RNG_class,                    only : RNG
  use tallyAdmin_class,             only : tallyAdmin
  use timer_mod,                    only : secToChar
  use transportOperator_inter,      only : transportOperator
  use universalVariables

  implicit none
  private

  !!
  !! Physics Package for fixed source calculations
  !!
  type, public, extends(particlePhysicsPackage) :: fixedSourcePhysicsPackage
    private
    ! Building blocks
    type(tallyAdmin), pointer                   :: tally => null()

    ! Settings
    integer(shortInt)                           :: bufferShift = 0

    ! Calculation components
    type(particleDungeon), pointer              :: commonBuffer => null()
  contains
    procedure :: collectSpecificResults
    procedure :: displayCycleProgress
    procedure :: getCycleParticlesNumber
    procedure :: getTallyAdminPtr
    procedure :: init
    procedure :: initCycle
    procedure :: kill
    procedure :: printSettings
    procedure :: processEndOfCycle
    procedure :: run
    procedure :: trackParticleHistory
  end type fixedSourcePhysicsPackage

contains
  !!
  !! Print calculation results to file
  !!
  subroutine collectSpecificResults(self, out)
    class(fixedSourcePhysicsPackage), intent(in) :: self
    type(outputFile), intent(inout)              :: out
    character(nameLen)                           :: name

    ! Call superclass.
    call collectSpecificResults_super(self, out)

    name = 'Source_batches'
    call out % printValue(self % getCyclesNumber(), name)

    ! Print tally
    call self % tally % print(out)

  end subroutine collectSpecificResults

  !!
  !!
  !!
  subroutine displayCycleProgress(self, cycleNumber, nInitialParticles, nFinalParticles, elapsedTime, endTime, timeToEnd)
    class(fixedSourcePhysicsPackage), intent(in) :: self
    integer(shortInt), intent(in)                :: cycleNumber, nInitialParticles, nFinalParticles
    real(defReal), intent(in)                    :: elapsedTime, endTime, timeToEnd

    ! Display progress
    call printFishLineR(cycleNumber)
    print *
    print *, 'Source batch: ', numToChar(cycleNumber), ' of ', numToChar(self % getCyclesNumber())
    print *, 'Pop:          ', numToChar(nInitialParticles)
    print *, 'Elapsed time: ', trim(secToChar(elapsedTime))
    print *, 'End time:     ', trim(secToChar(endTime))
    print *, 'Time to end:  ', trim(secToChar(timeToEnd))

  end subroutine displayCycleProgress

  !!
  !!
  !!
  function getCycleParticlesNumber(self) result(nParticles)
    class(fixedSourcePhysicsPackage), intent(in) :: self
    integer(shortInt)                            :: nParticles

    nParticles = self % getParticlesNumber()

  end function getCycleParticlesNumber

  !!
  !!
  !!
  function getTallyAdminPtr(self) result(tallyAdminPtr)
    class(fixedSourcePhysicsPackage), intent(in) :: self
    type(tallyAdmin), pointer                    :: tallyAdminPtr

    tallyAdminPtr => self % tally

  end function getTallyAdminPtr

  !!
  !! Initialise from individual components and dictionaries for source and tally
  !!
  subroutine init(self, payload)
    class(fixedSourcePhysicsPackage), intent(inout) :: self
    class(initPhysicsPackagePayload), intent(in)    :: payload
    type(initParticlePhysicsPackagePayload)         :: initPayload
    integer(shortInt)                               :: commonBufferSize
    character(*), parameter                         :: Here = 'init (fixedSourcePhysicsPackage_class.f90)'

    ! Create payload.
    call copyPayload(payload, initPayload)
    initPayload % defaultBufferSize = 50
    initPayload % isSourceRequired = .true.

    ! Initialise superclass.
    call init_super(self, initPayload)

    ! Initialise tally Admin
    allocate(self % tally)
    call self % tally % init(payload % dict % getDictPtr('tally'))

    ! Is the common buffer turned on? Set the size if so
    if (payload % dict % isPresent('commonBufferSize')) then
      call payload % dict % get(commonBufferSize,'commonBufferSize')
      allocate(self % commonBuffer)
      call self % commonBuffer % init(commonBufferSize)

      ! Set threshold at which to shift particles from private buffer
      ! to common buffer
      call payload % dict % getOrDefault(self % bufferShift, 'bufferShift', 10)
      if (self % getBufferSize() < self % bufferShift) &
      call fatalError(Here, 'Buffer size should be greater than the shift threshold')

    end if

    call self % printSettings()

  end subroutine init

  !!
  !!
  !!
  subroutine initCycle(self)
    class(fixedSourcePhysicsPackage), intent(inout) :: self

    ! Prepare source for current cycle.
    call self % generateSource()

  end subroutine initCycle

  !!
  !! Deallocate memory
  !!
  subroutine kill(self)
    class(fixedSourcePhysicsPackage), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    if (associated(self % tally)) then
      call self % tally % kill()
      deallocate(self % tally)

    end if
    self % bufferShift = 0

  end subroutine kill

  !!
  !! Print settings of the physics package
  !!
  subroutine printSettings(self)
    class(fixedSourcePhysicsPackage), intent(in) :: self

    print *, repeat("<>", 50)
    print *, "/\/\ FIXED SOURCE CALCULATION /\/\"
    print *, "Source batches:       ", numToChar(self % getCyclesNumber())
    print *, "Population per batch: ", numToChar(self % getParticlesNumber())
    print *, "Initial RNG Seed:     ", numToChar(self % getInitialSeed())
    print *
    print *, repeat("<>", 50)

  end subroutine printSettings

  !!
  !!
  !!
  subroutine processEndOfCycle(self, nFinalParticles)
    class(fixedSourcePhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(out)                  :: nFinalParticles
    type(RNG), pointer                              :: RNGPtr
    type(particleDungeon), pointer                  :: currentCyclePtr

    ! Update RNG.
    nFinalParticles = self % getParticlesNumber()
    RNGPtr => self % getRNGPtr()
    call RNGPtr % stride(nFinalParticles)

    ! Send end of cycle report.
    currentCyclePtr => self % getCurrentCyclePtr()
    call self % tally % reportCycleEnd(currentCyclePtr)

  end subroutine processEndOfCycle

  !!
  !!
  !!
  subroutine run(self)
    class(fixedSourcePhysicsPackage), intent(inout) :: self

    print *, repeat("<>", 50)
    print *, "/\/\ FIXED SOURCE CALCULATION /\/\"

    call self % runCycles(self % getCyclesNumber())
    call self % collectResults()

    print *
    print *, "\/\/ END OF FIXED SOURCE CALCULATION \/\/"
    print *

  end subroutine run

  !!
  !!
  !!
  subroutine trackParticleHistory(self, transOp, collOp, p, buffer, tally)
    class(fixedSourcePhysicsPackage), intent(in)        :: self
    class(transportOperator), intent(inout)             :: transOp
    type(collisionOperator), intent(inout)              :: collOp
    class(physicalParticle), allocatable, intent(inout) :: p
    type(particleDungeon), intent(inout)                :: buffer
    type(tallyAdmin), intent(inout)                     :: tally
    class(physicalParticle), allocatable                :: transfer
    integer(shortInt)                                   :: bufferExtra, i

    bufferLoop: do
      call p % setKEff(ONE)
      call self % placeCoord(p % getCoordsPtr())

      ! Save state
      call p % savePreHistoryState()

      ! Transport particle until its death
      history: do
        call transOp % transport(p, tally)
        if (p % getIsDead()) exit history

        call collOp % collide(p, tally, buffer, buffer)
        if (p % getIsDead()) exit history

      end do history

      ! If buffer is quite full, shift some particles to the commonBuffer
      if (associated(self % commonBuffer) .and. self % bufferShift < buffer % popSize()) then
        bufferExtra = buffer % popSize() - self % bufferShift
        do i = 1, bufferExtra
          call buffer % release(transfer)
          call self % commonBuffer % detainCritical(transfer)

        end do

      end if

      ! Clear out buffer
      if (.not. buffer % isEmpty()) then
        call buffer % release(p)

      elseif (associated(self % commonBuffer)) then
        ! Clear out common queue
        ! Note the apparently redundant critical sections (one here in PP, one in the dungeon).
        ! This is to prevent the situation where two threads both enter the conditional and compete
        ! for the final particle in the dungeon. The first thread would pop the particle while the
        ! second would try to pop from an empty dungeon.
        !$omp critical
        if (.not. self % commonBuffer % isEmpty()) call self % commonBuffer % releaseCritical(p)
        !$omp end critical
        if (p % getIsDead()) exit bufferLoop

      else
        exit bufferLoop

      end if

    end do bufferLoop

  end subroutine trackParticleHistory

end module fixedSourcePhysicsPackage_class