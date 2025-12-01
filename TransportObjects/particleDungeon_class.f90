module particleDungeon_class

  use CEParticleState_class,        only : CEParticleState
  use errors_mod,                   only : fatalError
  use genericProcedures,            only : numToChar
  use MGParticleState_class,        only : MGParticleState
  use numPrecision
  use physicalParticle_inter,       only : physicalParticle
  use physicalParticleFactory_func, only : new_physicalParticle
  use physicalParticleState_class,  only : castPhysicalParticleStatePtr, physicalParticleState, physicalParticleStateBox
  use RNG_class,                    only : RNG
  use transportObjectState_class,   only : transportObjectState

  implicit none
  private

  !!
  !! particleDungeon stores particle phase-space
  !! Used in eigenvalue calculation to store fission sites generated in a cycle
  !! Similar structures are refered to as:
  !! Store: MONK and Serpent(?)
  !! Fission Bank: OpenMC and MCNP(?)
  !!
  !! NOTE INCONSISTENT DEFINITIONS
  !! ****
  !! For convenience it allows storing the value of k-eff that can be retrieved to adjust fission site
  !! generation rate during a calculation. It is not currently used during normalisation but
  !! is used by analog k-eff tally. It is necessary to clarify behaviour.
  !! ****
  !!
  !! Technically not the whole particles are stored but only the key data defined in phaseCoord
  !!
  !! Dungeon can work like stacks or arrays. Stack-like behaviour is not really thread safe
  !! so it can be utilised when collecting and processing secondary particles in history
  !! that should be processed during the course of one cycle. Alternatively, one can use the
  !! critical variations of the stack-like procedures.
  !! Array-like behaviour allows to easily distribute particles among threads. As long as indices
  !! assigned to different threads do not overlap, reading is thread-safe (I hope-MAK).
  !!
  !!
  !! INTERFACE:
  !!   Stack-like interface:
  !!     detain(particle)          -> add a particle to the top
  !!     detainCritical(particle)  -> add a particle to the top with a critical operation
  !!     release(particle)         -> removes a particle from the top. Sets p % isDead = .false.
  !!     releaseCritical(particle) -> removes a particle from the top with a critical operation.
  !!                                   Sets p % isDead = .false.
  !!
  !!   Array-like interface:
  !!     replace(particle, i) -> overwrite prisoner data at index i
  !!     copy(particle, i)    -> copy prisoner at index i into particle. Sets p % isDead = .false.
  !!     get(i)               -> function returns particle state at index i
  !!
  !!   Misc procedures:
  !!     isEmpty()         -> returns .true. if there are no more particles
  !!     cleanPop()        -> kill or prisoners
  !!     normWeight(totWgt)-> normalise dungeon population so its total weight is totWgt
  !!     normSize(N)       -> normalise dungeon population so it contains N particles
  !!                          does not take nonuniform weight of particles into account
  !!     setSize(n)        -> sizes dungeon to have n dummy particles for ease of overwriting
  !!     printToFile(name) -> prints population in ASCII format to file "name"
  !!
  !!   Build procedures:
  !!     init(maxSize)     -> allocate space to store maximum of maxSize particles
  !!     kill()            -> return to uninitialised state
  !!
  type, public :: particleDungeon
    private
    integer(shortInt) :: pop = 0     ! Current population size of the dungeon
    real(defReal)     :: k_eff = ONE ! k-eff for fission site generation rate normalisation

    ! Storage space
    type(physicalParticleStateBox), dimension(:), allocatable :: prisoners
  contains
    generic            :: allocatePrisoner => allocatePrisoner_particle, allocatePrisoner_particleState
    procedure, private :: allocatePrisoner_particle
    procedure, private :: allocatePrisoner_particleState
    procedure, private :: checkIdx
    procedure, private :: checkPopulation
    procedure          :: cleanPop
    procedure          :: copy
    generic            :: detain => detain_particle, detain_particleState
    procedure, private :: detain_particle
    procedure, private :: detain_particleState
    generic            :: detainCritical => detainCritical_particle, detainCritical_particleState
    procedure, private :: detainCritical_particle
    procedure, private :: detainCritical_particleState
    procedure          :: get
    procedure          :: getKEff
    procedure          :: init
    procedure          :: isEmpty
    procedure          :: kill
    generic            :: killPrisoner => killPrisoner_shortInt, killPrisoner_shortIntArray
    procedure, private :: killPrisoner_shortInt
    procedure, private :: killPrisoner_shortIntArray
    procedure          :: normSize
    procedure          :: normWeight
    procedure          :: popSize
    procedure          :: popWeight
    procedure          :: printToFile
    procedure          :: release
    procedure          :: releaseCritical
    generic            :: replace => replace_particle, replace_particleState, replace_prisoner
    procedure, private :: replace_particle
    procedure, private :: replace_particleState
    procedure, private :: replace_prisoner
    procedure          :: setKEff
    procedure          :: setSize
    procedure          :: sortByBroodID
  end type particleDungeon

contains
  !!
  !!
  !!
  subroutine allocatePrisoner_particle(self, p, idx, replace)
    class(particleDungeon), intent(inout)    :: self
    class(physicalParticle), intent(in)      :: p
    integer(shortInt), intent(in)            :: idx
    logical(defBool), intent(in)             :: replace
    class(transportObjectState), allocatable :: stateCopy

    call p % copyCurrentState(stateCopy)
    call self % allocatePrisoner_particleState(stateCopy, idx, replace)

  end subroutine allocatePrisoner_particle

  !!
  !!
  !!
  subroutine allocatePrisoner_particleState(self, state, idx, replace)
    class(particleDungeon), intent(inout)   :: self
    class(transportObjectState), intent(in) :: state
    integer(shortInt), intent(in)           :: idx
    logical(defBool), intent(in)            :: replace

    if (replace) then
      call self % checkIdx(idx)

    else
      call self % checkPopulation(idx)

    end if

    call self % killPrisoner(idx)
    allocate(self % prisoners(idx) % ptr, source = castPhysicalParticleStatePtr(state))

  end subroutine allocatePrisoner_particleState

  !!
  !!
  !!
  subroutine checkIdx(self, idx)
    class(particleDungeon), intent(in) :: self
    integer(shortInt), intent(in)      :: idx
    character(*), parameter            :: here = 'checkIdx (particleDungeon_class.f90)'

    ! Protect against out-of-bounds access.
    if (idx < 1 .or. self % pop < idx) &
    call fatalError(here, 'Out of bounds access with index: '//numToChar(idx)//&
                          ' with particle population of: '//numToChar(self % pop)//'.')

  end subroutine checkIdx

  !!
  !!
  !!
  subroutine checkPopulation(self, population)
    class(particleDungeon), intent(in) :: self
    integer(shortInt), intent(in)      :: population
    integer(shortInt)                  :: nPrisoners
    character(*), parameter            :: here = 'checkPopulation (particleDungeon_class.f90)'

    ! Check for population overflow.
    nPrisoners = size(self % prisoners)
    if (nPrisoners < population) &
    call fatalError(here, 'Ran out of space for particles. Maximum size: '//numToChar(nPrisoners)//'.'&
                          'Current population: '//numToChar(self % pop)//'.')

  end subroutine checkPopulation

  !!
  !! Kill or particles in the dungeon
  !!
  pure subroutine cleanPop(self)
    class(particleDungeon), intent(inout) :: self

    self % pop = 0

  end subroutine cleanPop

  !!
  !! Copy particle from a location inside the dungeon
  !!
  !! Makes particle alive at exit. Also sets the broodID of the particle
  !! making it ready to be transported.
  !!
  !! Args:
  !!  p   [inout] -> Particle to be filled with data
  !!  idx [in]    -> Index of the particle to be copied
  !!
  !! Errors:
  !!  fatalError if requested index is 0, -ve or above current population
  !!
  function copy(self, idx) result(p)
    class(particleDungeon), intent(in)   :: self
    integer(shortInt), intent(in)        :: idx
    class(physicalParticle), allocatable :: p

    ! Generate new particle from state in dungeon then set its broodId.
    call self % checkIdx(idx)
    p = new_physicalParticle(self % prisoners(idx) % ptr)
    call p % setBroodId(idx)

  end function copy

  !!
  !! Store particle in the dungeon
  !!
  subroutine detain_particle(self, p)
    class(particleDungeon), intent(inout)    :: self
    class(physicalParticle), intent(in)      :: p

    !$omp atomic update
    ! Increase population and weight
    self % pop = self % pop + 1
    !$omp end atomic

    ! Load new state.
    call self % allocatePrisoner(p, self % pop, .false.)

  end subroutine detain_particle

  !!
  !! Store particle in the dungeon with a critical operation
  !!
  subroutine detainCritical_particle(self, p)
    class(particleDungeon), intent(inout)    :: self
    class(physicalParticle), intent(in)      :: p

    !$omp critical (dungeon)
    ! Increase population and weight
    self % pop = self % pop + 1
    
    ! Load new state.
    call self % allocatePrisoner(p, self % pop, .false.)
    !$omp end critical (dungeon)

  end subroutine detainCritical_particle

  !!
  !! Store phaseCoord in the dungeon
  !!
  subroutine detain_particleState(self, state)
    class(particleDungeon), intent(inout)   :: self
    class(transportObjectState), intent(in) :: state

    ! Increase population
    !$omp atomic update
    self % pop = self % pop + 1
    !$omp end atomic

    ! Load new state.
    call self % allocatePrisoner(state, self % pop, .false.)

  end subroutine detain_particleState

  !!
  !! Store phaseCoord in the dungeon with a critical operation
  !!
  subroutine detainCritical_particleState(self, state)
    class(particleDungeon), intent(inout)   :: self
    class(transportObjectState), intent(in) :: state

    ! Increase population
    !$omp critical (dungeon)
    self % pop = self % pop + 1

    ! Load new state.
    call self % allocatePrisoner(state, self % pop, .false.)
    !$omp end critical (dungeon)

  end subroutine detainCritical_particleState

  !!
  !! Return particleState from a location inside the dungeon
  !! Gives fatalError if requested index is 0, -ve or above current population
  !!
  function get(self, idx) result(statePtr)
    class(particleDungeon), intent(in)    :: self
    integer(shortInt), intent(in)         :: idx
    class(physicalParticleState), pointer :: statePtr

    call self % checkIdx(idx)
    statePtr => self % prisoners(idx) % ptr

  end function get

  !!
  !!
  !!
  elemental function getKEff(self) result(k_eff)
    class(particleDungeon), intent(in) :: self
    real(defReal)                      :: k_eff

    k_eff = self % k_eff

  end function getKEff

  !!
  !! Allocate space for the particles
  !!
  subroutine init(self,maxSize)
    class(particleDungeon), intent(inout) :: self
    integer(shortInt), intent(in)         :: maxSize

    ! Kill dungeon for safety then allocate space.
    call self % kill()
    allocate(self % prisoners(maxSize))

  end subroutine init

  !!
  !! Returns .true. if dungeon is empty
  !!
  elemental function isEmpty(self) result(isIt)
    class(particleDungeon), intent(in) :: self
    logical(defBool)                   :: isIt

    isIt = self % pop == 0

  end function isEmpty

  !!
  !!
  !!
  subroutine killPrisoner_shortInt(self, idx)
    class(particleDungeon), intent(inout) :: self
    integer(shortInt), intent(in)         :: idx

    if (associated(self % prisoners(idx) % ptr)) then
      call self % prisoners(idx) % ptr % kill()
      deallocate(self % prisoners(idx) % ptr)

    end if

  end subroutine killPrisoner_shortInt

  !!
  !!
  !!
  subroutine killPrisoner_shortIntArray(self, idxs)
    class(particleDungeon), intent(inout)       :: self
    integer(shortInt), dimension(:), intent(in) :: idxs
    integer(shortInt)                           :: i

    do i = 1, size(idxs)
      call self % killPrisoner_shortInt(idxs(i))

    end do

  end subroutine killPrisoner_shortIntArray

  !!
  !! Deallocate memory and return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(particleDungeon), intent(inout) :: self
    integer(shortInt)                     :: i

    ! Reset settings
    self % pop = 0
    self % k_eff = ONE

    ! Deallocate memeory.
    if (allocated(self % prisoners)) then
      do i = 1, size(self % prisoners)
        if (associated(self % prisoners(i) % ptr)) then
          call self % prisoners(i) % ptr % kill()
          deallocate(self % prisoners(i) % ptr)

        end if

      end do
      deallocate(self % prisoners)

    end if

  end subroutine kill

  !!
  !! Normalise total number of particles in the dungeon to match the provided number.
  !! Randomly duplicate or remove particles to match the number.
  !! Does not take weight of a particle into account!
  !!
  subroutine normSize(self, N, rand)
    class(particleDungeon), intent(inout)         :: self
    integer(shortInt), intent(in)                 :: N
    type(RNG), intent(inout)                     :: rand
    class(physicalParticleState), pointer         :: ptr
    integer(shortInt)                             :: excessP, i, idx, j, maxBroodId, nCopies, nDuplicates, nPrisoners
    integer(shortInt), dimension(:), allocatable  :: duplicates
    character(*), parameter                       :: here = 'normSize (particleDungeon_class.f90)'

    ! Protect against invalid N
    nPrisoners = size(self % prisoners)
    if (nPrisoners < N) then
      call fatalError(here, 'Requested size: '//numToChar(N)//&
                            ' is greater than maximum size: '//numToChar(nPrisoners)//'.')

    else if (N < 1) then
      call fatalError(here, 'Requested size: '//numToChar(N)//' is negative.')

    end if

    ! Return immediately if there are no prisoners in the dungeon (degenerate case).
    if (self % pop == 0) return

    ! Determine the maximum brood ID and sort the dungeon
    maxBroodId = 0
    do i = 1, self % pop
      maxBroodId = max(maxBroodId, self % prisoners(i) % ptr % getBroodId())

    end do
    call self % sortByBroodID(maxbroodId)

    ! Calculate excess particles to be removed
    excessP = self % pop - N

    if (0 < excessP) then ! Reduce population with reservoir sampling
      do i = N + 1, self % pop
        ! Select new index. Copy data if it is in the safe zone (<= N).
        call rand % generate(idx, i, 1)
        if (idx <= N) then
          ptr => self % prisoners(idx) % ptr
          self % prisoners(idx) % ptr => self % prisoners(i) % ptr
          self % prisoners(i) % ptr => ptr

        end if

      end do

      call self % killPrisoner([(i, i = N + 1, self % pop)])
      self % pop = N

    else if (excessP < 0) then ! Clone randomly selected particles
      ! For a massive undersampling duplicate (or n-plicate) particles
      excessP = -excessP
      nCopies = excessP / self % pop
      nDuplicates = modulo(excessP, self % pop)
      nPrisoners = self % pop

      ! Copy all particles maximum possible number of times
      do i = 1, nCopies
        do j = 1, nPrisoners
          call self % replace(j, nPrisoners * i + j)

        end do

      end do

      ! Choose the remainder particles to duplicate without replacement
      duplicates = [(i, i = 1, nDuplicates)]
      do i = nDuplicates + 1, self % pop
        call rand % generate(idx, i, 1)
        if (idx <= nDuplicates) duplicates(idx) = i

      end do
      self % pop = self % pop * (nCopies + 1)

      ! Copy the duplicated particles at the end
      do i = 1, nDuplicates
        call self % replace(duplicates(i), self % pop + i)

      end do
      self % pop = N

    end if

  end subroutine normSize

  !!
  !! Normalise total weight of the particles in the dungeon to match provided value
  !!
  subroutine normWeight(self, value)
    class(particleDungeon), intent(inout) :: self
    real(defReal), intent(in)             :: value
    integer(shortInt)                     :: i
    real(defReal), dimension(self % pop)  :: weights

    do i = 1, self % pop
      weights(i) = self % prisoners(i) % ptr % getWeight()

    end do
    weights = weights * value / sum(weights)

    do i = 1, self % pop
      call self % prisoners(i) % ptr % setWeight(weights(i))

    end do

  end subroutine normWeight

  !!
  !! Returns number of neutrons in the dungeon
  !!
  function popSize(self) result(pop)
    class(particleDungeon), intent(in) :: self
    integer(shortInt)                  :: pop

    pop = self % pop

  end function popSize

  !!
  !! Returns total population weight
  !!
  function popWeight(self) result(wgt)
    class(particleDungeon), intent(in) :: self
    integer(shortInt)                  :: i
    real(defReal)                      :: wgt

    wgt = ZERO
    do i = 1, self % pop
      wgt = wgt + self % prisoners(i) % ptr % getWeight()

    end do

  end function popWeight
  
  !!
  !! Prints the position of fission sites to a file
  !! Used initially for looking at clustering
  !!
  subroutine printToFile(self, name)
    class(particleDungeon), intent(in) :: self
    character(*), intent(in)           :: name
    character(256)                     :: filename
    integer(shortInt)                  :: i
    integer(shortInt), parameter       :: unit = 10

    filename = trim(name)//'.txt'
    open(unit = unit, file = filename, status = 'new')

    ! Print out each particle co-ordinate
    do i = 1, self % pop
      write(unit, *) self % prisoners(i) % ptr % getGlobalPosition()
      write(unit, *) self % prisoners(i) % ptr % getGlobalDirection()
      select type(ptr => self % prisoners(i) % ptr)
        class is(CEParticleState)
          write(unit, *) ptr % getEnergy()

        class is(MGParticleState)
          write(unit, *) ptr % getEnergyGroup()

        class default
          ! Do nothing.

      end select
      write(unit, *) self % prisoners(i) % ptr % getBroodId()

    end do

    ! Close the file
    close(unit)

  end subroutine printToFile

  !!
  !! Pop the particle from the top of the dungeon.
  !! Makes particle alive at exit
  !!
  subroutine release(self, p)
    class(particleDungeon), intent(inout)             :: self
    class(physicalParticle), allocatable, intent(out) :: p
    integer(shortInt)                                 :: pop

    !$omp atomic capture
    ! Decrease population
    pop = self % pop
    self % pop = self % pop - 1
    !$omp end atomic

    ! Generate new particle from state then free space in dungeon.
    p = new_physicalParticle(self % prisoners(pop) % ptr)
    call self % killPrisoner(pop)

  end subroutine release

  !!
  !! Pop the particle from the top of the dungeon with a critical operation.
  !! Makes particle alive at exit
  !!
  subroutine releaseCritical(self, p)
    class(particleDungeon), intent(inout)             :: self
    class(physicalParticle), allocatable, intent(out) :: p
    integer(shortInt)                                 :: pop

    !$omp critical (dungeon)
    ! Decrease population
    pop = self % pop
    self % pop = self % pop - 1

    ! Generate new particle from state.
    p = new_physicalParticle(self % prisoners(pop) % ptr)
    !$omp end critical (dungeon)

    ! Free space in dungeon.
    call self % killPrisoner(pop)

  end subroutine releaseCritical

  !!
  !! Replace data of particle prisoner at the index idx with particle
  !!
  subroutine replace_particle(self, p, idx)
    class(particleDungeon), intent(inout)    :: self
    class(physicalParticle), intent(in)      :: p
    integer(shortInt), intent(in)            :: idx

    ! Load new particle
    call self % allocatePrisoner(p, idx, .true.)

  end subroutine replace_particle

  !!
  !! Replace data of particle prisoner at the index idx with phaseCoords
  !!
  subroutine replace_particleState(self, state, idx)
    class(particleDungeon), intent(inout)   :: self
    class(transportObjectState), intent(in) :: state
    integer(shortInt), intent(in)           :: idx

    ! Load new particle
    call self % allocatePrisoner(state, idx, .true.)

  end subroutine replace_particleState

  !!
  !!
  !!
  subroutine replace_prisoner(self, sourceIdx, targetIdx)
    class(particleDungeon), intent(inout) :: self
    integer(shortInt), intent(in)         :: sourceIdx, targetIdx
    character(*), parameter               :: HERE = 'replace_prisoner (particleDungeon_class.f90)'

    ! Check that source is associated first.
    if (.not. associated(self % prisoners(sourceIdx) % ptr)) &
    call fatalError(HERE, 'Attempting to replace prisoner by unassociated pointer.')

    ! Kill prisoner at targetIdx then copy.
    call self % killPrisoner(targetIdx)
    allocate(self % prisoners(targetIdx) % ptr, source = self % prisoners(sourceIdx) % ptr)

  end subroutine replace_prisoner

  !!
  !!
  !!
  elemental subroutine setKEff(self, k_eff)
    class(particleDungeon), intent(inout) :: self
    real(defReal), intent(in)             :: k_eff

    self % k_eff = k_eff

  end subroutine setKEff

  !!
  !! Set size of the dungeon to n
  !!
  !! Sets population to arbitrary size n
  !! All stored particles revert to default initialisation state
  !!
  !! Args:
  !!   n [in] -> Requested size of the population
  !!
  !! Errors:
  !!   fatalError if n is invalid (not +ve)
  !!
  subroutine setSize(self, n)
    class(particleDungeon), intent(inout) :: self
    integer(shortInt), intent(in)         :: n
    integer(shortInt)                     :: i, nPrisoners
    character(*), parameter               :: here = 'setSize (particleDungeon_class.f90)'

    if (n < 1) call fatalError(here, 'Requested population: '//numToChar(n)//' is negative.')

    ! Set population
    self % pop = n

    ! Make sure enough space is available
    if (allocated(self % prisoners)) then
      nPrisoners = size(self % prisoners)
      if (nPrisoners < n) then
        call self % killPrisoner([(i, i = 1, nPrisoners)])
        deallocate(self % prisoners)
        allocate(self % prisoners(n))

      end if

    else
      allocate(self % prisoners(n))

    end if

  end subroutine setSize

  !!
  !! Reorder the dungeon so the brood ID is in the ascending order
  !!
  !! Args:
  !!   maxBroodId [in] -> Maximum brood ID
  !!
  subroutine sortByBroodID(self, maxBroodId)
    class(particleDungeon), intent(inout)         :: self
    integer(shortInt), intent(in)                 :: maxBroodId
    integer(shortInt), dimension(maxBroodId)      :: counts
    integer(shortInt)                             :: broodId, i, loc, c
    integer(shortInt), dimension(self % pop)      :: perm
    class(physicalParticleState), pointer         :: ptr
    character(*), parameter :: Here = 'sortBybroodID (particleDungeon_class.f90)'

    ! Count number of particles with each brood ID
    counts = 0
    do i = 1, self % pop
      broodId = self % prisoners(i) % ptr % getBroodId()
      if (broodId < 1 .or. maxBroodId < broodId) call fatalError(Here, 'Out of range brood id: '//numToChar(broodId)//'.')
      counts(broodId) = counts(broodId) + 1

    end do

    ! Convert to starting index
    loc = 1
    do i = 1, maxBroodId
      c = counts(i)
      counts(i) = loc
      loc = loc + c

    end do

    ! Create the permutation array
    do i = 1, self % pop
      broodId = self % prisoners(i) % ptr % getBroodId()
      loc = counts(broodId)
      counts(broodId) = counts(broodId) + 1
      perm(loc) = i

    end do

    ! Permute particles
    do i = 1, self % pop
      loc = perm(i)

      ! If the element was already swapped follow it to its location
      do while (loc < i)
        loc = perm(loc)

      end do

      ! Swap elements
      if (loc /= i) then
        ptr => self % prisoners(i) % ptr
        self % prisoners(i) % ptr => self % prisoners(loc) % ptr
        self % prisoners(loc) % ptr => ptr

      end if

    end do

  end subroutine sortByBroodID

end module particleDungeon_class