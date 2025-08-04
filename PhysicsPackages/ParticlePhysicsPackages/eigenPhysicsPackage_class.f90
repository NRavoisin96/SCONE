module eigenPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,            only : fatalError, printFishLineR, numToChar
  use dictionary_class,             only : dictionary
  use outputFile_class,             only : outputFile

  ! Timers
  use timer_mod,                    only : secToChar

  ! Particle classes and Random number generator
  use particle_class,               only : particle
  use particleDungeon_class,        only : particleDungeon
  use RNG_class,                    only : RNG

  ! Physics package interface
  use particlePhysicsPackage_inter, only : collectSpecificResults_super => collectSpecificResults, &
                                           init_super => init, initParticlePhysicsPackagePayload, &
                                           particlePhysicsPackage, kill_super => kill
  use physicsPackage_inter,         only : copyPayload, initPhysicsPackagePayload

  ! Geometry
  use geometryReg_mod,              only : gr_fieldIdx => fieldIdx, gr_fieldPtr => fieldPtr

  ! Fields
  use field_inter,                  only : field
  use uniFissSitesField_class,      only : uniFissSitesField, uniFissSitesField_TptrCast
  use fieldFactory_func,            only : new_field

  ! Operators
  use collisionOperator_class,      only : collisionOperator
  use transportOperator_inter,      only : transportOperator

  ! Tallies
  use tallyAdmin_class,             only : tallyAdmin
  use tallyResult_class,            only : tallyResult
  use keffAnalogClerk_class,        only : keffResult

  implicit none
  private

  !!
  !! Physics Package for eigenvalue calculations
  !!
  type, public, extends(particlePhysicsPackage) :: eigenPhysicsPackage
    private
    ! Building blocks
    type(tallyAdmin), pointer                   :: inactiveTally => null()
    type(tallyAdmin), pointer                   :: activeTally => null()
    type(tallyAdmin), pointer                   :: inactiveAtch => null()
    type(tallyAdmin), pointer                   :: activeAtch => null()
    class(uniFissSitesField), pointer           :: ufsField => null()

    ! Settings
    integer(shortInt)                           :: N_inactive = 0
    real(defReal)                               :: k_eff = ONE
    logical(defBool)                            :: inactiveCycles = .true., UFS = .false.

    ! Calculation components
    type(particleDungeon), pointer              :: nextCycle => null()
  contains
    procedure :: collectSpecificResults
    procedure :: displayCycleProgress
    procedure :: generateInitialState
    procedure :: getCycleParticlesNumber
    procedure :: getTallyAdminPtr
    procedure :: init
    procedure :: kill
    procedure :: printSettings
    procedure :: processEndOfCycle
    procedure :: run
    procedure :: trackParticleHistory
  end type eigenPhysicsPackage

contains
  !!
  !! Print calculation results to file
  !!
  subroutine collectSpecificResults(self, out)
    class(eigenPhysicsPackage), intent(in) :: self
    type(outputFile), intent(inout)        :: out
    character(nameLen)                     :: name

    ! Call superclass.
    call collectSpecificResults_super(self, out)

    name = 'Inactive_Cycles'
    call out % printValue(self % N_inactive,name)

    name = 'Active_Cycles'
    call out % printValue(self % getCyclesNumber(), name)

    ! Print Inactive tally
    name = 'inactive'
    call out % startBlock(name)
    call self % inactiveTally % print(out)
    call out % endBlock()

    ! Print Active attachment
    ! Is printed into the root block
    call self % activeAtch % print(out)

    name = 'active'
    call out % startBlock(name)
    call self % activeTally % print(out)
    call out % endBlock()

  end subroutine collectSpecificResults

  !!
  !!
  !!
  subroutine displayCycleProgress(self, cycleNumber, nInitialParticles, nFinalParticles, elapsedTime, endTime, timeToEnd)
    class(eigenPhysicsPackage), intent(in) :: self
    integer(shortInt), intent(in)          :: cycleNumber, nInitialParticles, nFinalParticles
    real(defReal), intent(in)              :: elapsedTime, endTime, timeToEnd

    ! Display progress
    call printFishLineR(cycleNumber)
    print *
    print *, 'Cycle: ', numToChar(cycleNumber), ' of ', &
    merge(numToChar(self % N_inactive), numToChar(self % getCyclesNumber()), self % inactiveCycles)
    print *, 'Pop: ', numToChar(nInitialParticles) , ' -> ', numToChar(nFinalParticles)
    print *, 'Elapsed time: ', trim(secToChar(elapsedTime))
    print *, 'End time:     ', trim(secToChar(endTime))
    print *, 'Time to end:  ', trim(secToChar(timeToEnd))

  end subroutine displayCycleProgress

  !!
  !!
  !!
  subroutine generateInitialState(self)
    class(eigenPhysicsPackage), intent(inout) :: self

    ! Allocate and initialise nextCycle dungeon.
    allocate(self % nextCycle)
    call self % nextCycle % init(3 * self % getParticlesNumber())

    ! Generate initial source
    print *, "GENERATING INITIAL FISSION SOURCE"
    call self % generateSource()
    print *, "DONE!"

  end subroutine generateInitialState

  !!
  !!
  !!
  function getCycleParticlesNumber(self) result(nParticles)
    class(eigenPhysicsPackage), intent(in) :: self
    integer(shortInt)                      :: nParticles
    type(particleDungeon), pointer         :: currentCyclePtr

    currentCyclePtr => self % getCurrentCyclePtr()
    nParticles = currentCyclePtr % popSize()

  end function getCycleParticlesNumber

  !!
  !!
  !!
  function getTallyAdminPtr(self) result(tallyAdminPtr)
    class(eigenPhysicsPackage), intent(in) :: self
    type(tallyAdmin), pointer              :: tallyAdminPtr

    if (self % inactiveCycles) then
      tallyAdminPtr => self % inactiveTally

    else
      tallyAdminPtr => self % activeTally

    end if

  end function getTallyAdminPtr


  !!
  !! Initialise from individual components and dictionaries for inactive and active tally
  !!
  subroutine init(self, payload)
    class(eigenPhysicsPackage), intent(inout)    :: self
    class(initPhysicsPackagePayload), intent(in) :: payload
    type(initParticlePhysicsPackagePayload)      :: initPayload
    type(dictionary)                             :: locDict1, locDict2
    class(field), pointer                        :: field
    character(*), parameter                      :: Here = 'init (eigenPhysicsPackage_class.f90)'

    ! Create payload from input.
    call copyPayload(payload, initPayload)
    initPayload % defaultBufferSize = 1000
    initPayload % currentCycleSizeMultiplier = 3

    ! Initialise superclass.
    call init_super(self, initPayload)

    ! Read calculation settings
    call payload % dict % get(self % N_inactive, 'inactive')

    ! Initial k_effective guess
    call payload % dict % getOrDefault(self % k_eff, 'keff_0', ONE)

    ! Read uniform fission site option as a geometry field
    if (payload % dict % isPresent('uniformFissionSites')) then
      self % ufs = .true.
      ! Build and initialise
      call new_field(payload % dict % getDictPtr('uniformFissionSites'), nameUFS)
      ! Save UFS field
      field => gr_fieldPtr(gr_fieldIdx(nameUFS))
      self % ufsField => uniFissSitesField_TptrCast(field)
      ! Initialise
      call self % ufsField % estimateVol(payload % geometry, self % getParticleType(), self % getRNGPtr())

    end if

    ! Initialise active & inactive tally Admins
    allocate(self % inactiveTally)
    call self % inactiveTally % init(payload % dict % getDictPtr('inactiveTally'))

    allocate(self % activeTally)
    call self % activeTally % init(payload % dict % getDictPtr('activeTally'))

    ! Initialise active and inactive tally attachments
    ! Inactive tally attachment
    call locDict1 % init(2)
    call locDict2 % init(2)

    call locDict2 % store('type','keffAnalogClerk')
    call locDict1 % store('keff', locDict2)
    call locDict1 % store('display',['keff'])

    allocate(self % inactiveAtch)
    call self % inactiveAtch % init(locDict1)

    call locDict2 % kill()
    call locDict1 % kill()

    ! Active tally attachment
    call locDict1 % init(2)
    call locDict2 % init(2)

    call locDict2 % store('type','keffImplicitClerk')
    call locDict1 % store('keff', locDict2)
    call locDict1 % store('display',['keff'])

    allocate(self % activeAtch)
    call self % activeAtch % init(locDict1)

    call locDict2 % kill()
    call locDict1 % kill()

    ! Attach attachments to result tallies
    call self % inactiveTally % push(self % inactiveAtch)
    call self % activeTally % push(self % activeAtch)

    call self % printSettings()

  end subroutine init

  !!
  !! Deallocate memory
  !!
  subroutine kill(self)
    class(eigenPhysicsPackage), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local-
    if (associated(self % inactiveTally)) then
      call self % inactiveTally % kill()
      deallocate(self % inactiveTally)

    end if

    if (associated(self % activeTally)) then
      call self % activeTally % kill()
      deallocate(self % activeTally)

    end if

    if (associated(self % inactiveAtch)) then
      call self % inactiveAtch % kill()
      deallocate(self % inactiveAtch)

    end if

    if (associated(self % activeAtch)) then
      call self % activeAtch % kill()
      deallocate(self % activeAtch)
 
    end if

    if (associated(self % ufsField)) then
      call self % ufsField % kill()
      deallocate(self % ufsField)

    end if

    if (associated(self % nextCycle)) then
      call self % nextCycle % kill()
      deallocate(self % nextCycle)

    end if

    ! Reset local values.
    self % N_inactive = 0
    self % k_eff = ONE
    self % inactiveCycles = .true.
    self % UFS = .false.

  end subroutine kill

  !!
  !! Print settings of the physics package
  !!
  subroutine printSettings(self)
    class(eigenPhysicsPackage), intent(in) :: self

    print *, repeat("<>", 50)
    print *, "/\/\ EIGENVALUE CALCULATION WITH POWER ITERATION METHOD /\/\"
    print *, "Inactive Cycles:     ", numToChar(self % N_inactive)
    print *, "Active Cycles:       ", numToChar(self % getCyclesNumber())
    print *, "Particle Population: ", numToChar(self % getParticlesNumber())
    print *, "Initial RNG Seed:    ", numToChar(self % getInitialSeed())
    print *
    print *, repeat("<>", 50)

  end subroutine printSettings

  !!
  !!
  !!
  subroutine processEndOfCycle(self, nFinalParticles)
    class(eigenPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(out)            :: nFinalParticles
    type(tallyAdmin), pointer                 :: attachmentPtr, tallyAdminPtr
    type(particleDungeon), pointer            :: currentCyclePtr
    type(RNG), pointer                        :: RNGPtr
    integer(shortInt)                         :: nParticles
    type(particleDungeon), pointer            :: tempCycle
    class(tallyResult), allocatable           :: result
    character(*), parameter                   :: here = 'processEndOfCycle (eigenPhysicsPackage_class.f90)'

    ! Get correct tally pointers.
    if (self % inactiveCycles) then
      attachmentPtr => self % inactiveAtch
      tallyAdminPtr => self % inactiveTally

    else
      attachmentPtr => self % activeAtch
      tallyAdminPtr => self % activeTally

    end if

    ! Clean up source bank from cycle that just finished.
    currentCyclePtr => self % getCurrentCyclePtr()
    call currentCyclePtr % cleanPop()

    ! Update main RNG stream to be ready for next cycle.
    RNGPtr => self % getRNGPtr()
    nParticles = self % getParticlesNumber()
    call RNGPtr % stride(nParticles + 1)

    ! Send end of cycle report to tally.
    nFinalParticles = self % nextCycle % popSize()
    call tallyAdminPtr % reportCycleEnd(self % nextCycle)

    ! Update UFS if used.
    if (self % UFS) call self % ufsField % updateMap()

    ! Normalise population then flip cycles.
    call self % nextCycle % normSize(nParticles, RNGPtr)
    tempCycle => self % nextCycle
    self % nextCycle => currentCyclePtr
    call self % setCurrentCyclePtr(tempCycle)

    ! Get new k_eff.
    call attachmentPtr % getResult(result, 'keff')

    ! Downcast result to correct type.
    select type(result)
      type is(keffResult)
        self % k_eff = result % keff(1)

      class default
        call fatalError(here, 'Invalid result type.')

    end select
    self % nextCycle % k_eff = self % k_eff

  end subroutine processEndOfCycle

  !!
  !!
  !!
  subroutine run(self)
    class(eigenPhysicsPackage), intent(inout) :: self

    print *, repeat("<>", 50)
    print *, "/\/\ EIGENVALUE CALCULATION /\/\"

    call self % generateInitialState()
    call self % runCycles(self % N_inactive)
    self % inactiveCycles = .false.
    call self % runCycles(self % getCyclesNumber())
    call self % collectResults()

    print *
    print *, "\/\/ END OF EIGENVALUE CALCULATION \/\/"
    print *

  end subroutine run

  !!
  !!
  !!
  subroutine trackParticleHistory(self, transOp, collOp, p, buffer, tally)
    class(eigenPhysicsPackage), intent(in)  :: self
    class(transportOperator), intent(inout) :: transOp
    type(collisionOperator), intent(inout)  :: collOp
    type(particle), intent(inout)           :: p
    type(particleDungeon), intent(inout)    :: buffer
    type(tallyAdmin), intent(inout)         :: tally

    bufferLoop: do
      ! Initialize the particle's state for this history.
      p % k_eff = self % k_eff
      call self % placeCoord(p % coords)
      call p % savePreHistory()

      ! Transport the particle until it dies
      history: do
        call transOp % transport(p, tally)
        if (p % isDead) exit history

        ! CRITICAL: Secondaries are sent to self % nextCycle
        call collOp % collide(p, tally, buffer, self % nextCycle)
        if (p % isDead) exit history

      end do history

      ! Check the local buffer for a secondary particle to continue the history
      if (buffer % isEmpty()) then
        exit bufferLoop ! The entire family history is complete

      else
        call buffer % release(p) ! Get the next particle from the buffer

      end if

    end do bufferLoop

  end subroutine trackParticleHistory

end module eigenPhysicsPackage_class