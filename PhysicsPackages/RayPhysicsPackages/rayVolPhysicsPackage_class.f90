module rayVolPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,    only : fatalError, numToChar, rotateVector, printFishLineR
  use hashFunctions_func,   only : FNV_1
  use dictionary_class,     only : dictionary
  use rng_class,            only : RNG
  use physicsPackage_inter, only : init_super => init, initPhysicsPackagePayload, kill_super => kill, physicsPackage
  use outputFile_class,     only : outputFile

  ! Timers
  use timer_mod,            only : registerTimer, timerStart, timerStop, timerTime, timerReset, secToChar

  ! Geometry
  use coordList_class,      only : coordList
  use geometry_inter,       only : geometry, distCache
  use geometryReg_mod,      only : gr_geomPtr => geomPtr, gr_geomIdx => geomIdx, gr_kill => kill
  use geometryFactory_func, only : new_geometry

  ! Nuclear Data
  use materialMenu_mod,     only : mm_nMat => nMat, mm_matName => matName
  use nuclearDataReg_mod,   only : ndReg_init => init, ndReg_getMatNames => getMatNames

  implicit none
  private

  ! Parameters
  integer(shortInt), parameter :: SCORE = 1, CSUM = 2, CSUM2 = 3

  !!
  !! Physics package to perform ray-tracking based volume calculation
  !!
  !! Calculates relative volume of diffrent materials in the problem by performing
  !! random ray tracing in the geometry. The volume is normalised that the total domain
  !! volume is 1.0.
  !!
  !! Rays travel by a random, exponentially distributed distance with a user-defined mean free
  !! path. After each segment, the ray scatters isotropically. At each scattering event
  !! ray mat be terminated with provided probability `abs_prob`. Ray will also be killed
  !! if it has reached the OUTSIDE material.
  !!
  !! Simulation can be performed in the ROBUST mode, which is intended to be used for debugging.
  !! If it is enabled, for each ray segment the material idx at the middle of the segment is
  !! obtained and checked against what is stored in coords. If they do not match it means that
  !! for some reason, the geometry stopped to track material composition correctly for the ray.
  !!
  !! This Physics Package exists to serve as a geometry debugging and benchmarking tool.
  !! The calculation it preforms is unlikley to serve any practical purpose.
  !!
  !! Sample Input Dictionary:
  !!   PP {
  !!     type rayVolPhysicsPackage;
  !!     mfp 0.3;       // Mean length of ray segments
  !!     abs_prob 0.1;  // Ray absorbtion probability after each segment
  !!     pop 2000;      // Number of rays per cycle
  !!     cycles 100;    // Number of cycles
  !!     robust 1;      // 1 for true; 0 for false; Enable robust mode
  !!     cache  1;      // 1 for treu; 0 for false; Enable distance caching
  !!     #seed 86868;#  // Optional RNG seed
  !!     geometry {<Geometry Definition>}
  !!     nuclearData {<Nuclear data definition. Requires material names only>}
  !!   }
  !!
  !! Private Members
  !!   geom      -> Pointer to the geometry
  !!   geomIdx   -> Index of the geometry in geometry Registry
  !!   rand      -> Random number generator
  !!   timerMain -> Index of the timer defined to measure calculation time
  !!   mfp       -> Mean length of the ray segment
  !!   abs_prob  -> Ray absorption probability after every segment
  !!   N_cycles  -> Number of cycles
  !!   robust    -> Flag to enable/disable robust mode
  !!   cache     -> Flag to enable/disable distance caching
  !!   res       -> Array to accumulate total track length in each material. Contains
  !!     score (SCORE), cumulative sum over cycles (CSUM) and cumulative sume of squares over
  !!     cycles (CSUM2)
  !!   totDist   -> Accumulator for total track distance in a single cycle
  !!   ray_speed -> Cumulative sume and cumulative sum of squares over cycles for the mean speed
  !!     of ray. Result in meters over CPU second.
  !!
  !! Interface:
  !!   physicsPackage interface
  !!
  type, public, extends(physicsPackage)         :: rayVolPhysicsPackage
    private
    ! Components
    type(RNG)                                   :: rand

    ! Settings
    real(defReal)                               :: abs_prob = ZERO, mfp = ZERO
    logical(defBool)                            :: cache = .false., robust = .false.

    ! Results space
    real(defReal), dimension(:, :), allocatable :: res
    real(defReal)                               :: totDist = ZERO
    real(defReal), dimension(2)                 :: ray_speed = ZERO

  contains
    ! Superclass procedures
    procedure :: collectSpecificResults
    procedure :: init
    procedure :: run
    procedure :: kill

    ! Private procedures
    procedure, private :: cycles
    procedure, private :: trackRay
    procedure, private :: printResults
    procedure, private :: printSettings
  end type rayVolPhysicsPackage

contains
  !!
  !!
  !!
  subroutine collectSpecificResults(self, out)
    class(rayVolPhysicsPackage), intent(in) :: self
    type(outputFile), intent(inout)         :: out

    ! Do nothing for now.

  end subroutine collectSpecificResults

  !!
  !! Initialise Physics Package from dictionary
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine init(self,payload)
    class(rayVolPhysicsPackage), intent(inout)   :: self
    class(initPhysicsPackagePayload), intent(in) :: payload
    character(*), parameter                      :: Here = 'init (rayVolPhysicsPackage_class.f90)'

    ! Initialise superclass.
    call init_super(self, payload)

    ! Load settings
    call payload % dict % get(self % mfp, 'mfp')
    call payload % dict % get(self % abs_prob, 'abs_prob')
    call payload % dict % get(self % robust, 'robust')
    call payload % dict % get(self % cache, 'cache')

    ! Check settings
    if (self % mfp < ZERO) then
      call fatalError(Here, 'Was given negative mean free path (mfp): '//numToChar(self % mfp)//'.')

    else if (self % abs_prob <= ZERO .or. ONE < self % abs_prob) then
      call fatalError(Here, 'Absorption probability is outside range [0, 1): '//numToChar(self % abs_prob)//'.')

    end if

    ! Initialise RNG
    call self % rand % init(self % getInitialSeed())

    ! Allocate results space
    allocate(self % res(mm_nMat(), 3))
    self % res = ZERO

  end subroutine init

  !!
  !! Run calculation
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine run(self)
    class(rayVolPhysicsPackage), intent(inout) :: self

    call self % printSettings()
    call self % cycles(self % rand)
    call self % printResults()

  end subroutine run

  !!
  !! Perform cycles of the stochatic volume calculation with ray tracing
  !!
  !! Randomly places the starting point based on uniform distribution.
  !!
  !! Args:
  !!   rand [inout] -> Initialised random number generator
  !!
  !! NOTE:
  !!   RNG needs to be given as an argument `class(RNG)` to prevent inlining. Compiler (gcc 8.3)
  !!   produced erroneous code withou it. Same random number would be produced for diffrent calls
  !!   of `get` function.
  !!
  subroutine cycles(self, rand)
    class(rayVolPhysicsPackage), intent(inout) :: self
    class(RNG), intent(inout)                  :: rand
    type(coordList)                            :: coords
    real(defReal), dimension(3)                :: randomNumbers, bottom, top
    real(defReal), dimension(3)                :: r, u
    real(defReal)                              :: mu, phi
    integer(shortInt)                          :: gen, ray, matIdx, uniqueId, i, nCycles, nParticles, timerMain
    type(RNG), save                            :: pRNG
    real(defReal)                              :: elapsed_T, end_T, T_toEnd, av_speed, cycle_T
    character(*), parameter                    :: Here = 'cycles (rayVolPhysicsPackage_class.f90)'
    !$omp threadprivate(pRNG)

    !$omp parallel
    pRNG = rand
    !$omp end parallel

    ! Reset and start timer
    timerMain = self % getTimerMain()
    call timerReset(timerMain)
    call timerStart(timerMain)

    ! Get lower an upper corner of bounding box
    associate (aabb => self % getGeometryBounds())
      bottom = aabb(1:3)
      top    = aabb(4:6)
    end associate

    ! Perform clculation
    nCycles = self % getCyclesNumber()
    nParticles = self % getParticlesNumber()
    do gen = 1, nCycles
      !$omp parallel do private(r, u, mu, phi, i, randomNumbers, matIdx, uniqueID, coords)
      do ray = 1, nParticles

        ! Set seed
        call pRNG % stride((gen - 1) * nParticles + ray)

        ! Find starting point that is inside the geometry
        i = 0
        call rand % generateMu(mu)
        call rand % generatePhi(phi)
        u = rotateVector([ONE, ZERO, ZERO], mu, phi)

        rejection : do
          call rand % generate(randomNumbers)
          r = bottom + (top - bottom) * randomNumbers

          ! Exit if point is inside the geometry
          call self % whatIsAt(r, u, uniqueId, matIdx)
          if (matIdx /= OUTSIDE_MAT) exit rejection

          i = i + 1
          if (i > 1000) then
            call fatalError(Here, 'Infinate loop when searching ray start in the geometry.')
          end if
        end do rejection

        ! Place in the geometry & process the ray
        call coords % init(r, u)
        call self % placeCoord(coords)
        call self % trackRay(coords, pRNG)

      end do
      !$omp end parallel do

      ! Calculate times
      call timerStop(timerMain)
      cycle_T = timerTime(timerMain) - elapsed_T
      elapsed_T = timerTime(timerMain)

      ! Predict time to end
      end_T = real(nCycles, defReal) * elapsed_T / gen
      T_toEnd = max(ZERO, end_T - elapsed_T)

      ! Calculate average tracking speed
      av_speed = self % totDist / cycle_T * 1.0E-3_defReal

      ! Display progress
      call printFishLineR(gen)
      print *
      print *, 'Cycle: ', numToChar(gen), ' of ', numToChar(nCycles)
      print *, 'Pop: ', numToChar(nParticles)
      print '(A, ES12.5)', ' Av. Ray speed: [m/s]: ', av_speed
      print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
      print *, 'End time:     ', trim(secToChar(end_T))
      print *, 'Time to end:  ', trim(secToChar(T_toEnd))

      ! Process scores
      self % res(:, SCORE) = self % res(:, SCORE) / self % totDist
      self % res(:, CSUM)  = self % res(:, CSUM)  + self % res(:, SCORE)
      self % res(:, CSUM2) = self % res(:, CSUM2) + self % res(:, SCORE)**2
      self % res(:, SCORE) = ZERO
      self % totDist = ZERO

      ! Average ray speed
      self % ray_speed(1) = self % ray_speed(1) + av_speed
      self % ray_speed(2) = self % ray_speed(2) + av_speed * av_speed

    end do

  end subroutine cycles

  !!
  !! Track single ray through the geometry
  !!
  !! Seperate this functionality in separate function to improve
  !! clarity of `cycles` procedure.
  !!
  !! Moves ray through the geometry untill it leaks or is absorbed.
  !! Accumulates track-length information along the way.
  !!
  !! Args:
  !!   coords [inout] -> Coordinates of ray placed in the geometry.
  !!   rand [inout] -> Initialised random number generator
  !!
  subroutine trackRay(self, coords, rand)
    class(rayVolPhysicsPackage), intent(inout) :: self
    class(RNG), intent(inout)                  :: rand
    type(coordList), intent(inout)             :: coords
    real(defReal)                              :: distance, mu, phi, maxDist, randomNumber
    real(defReal), dimension(3)                :: r, r_pre, u_pre
    integer(shortInt)                          :: event, matIdx, uniqueId, mat_mid, unique_mid
    type(distCache)                            :: cache_space
    character(*), parameter :: Here = 'trackRay (rayVolPhysicsPackage_class.f90)'

    ! Keep compiler happy
    r_pre = ZERO
    r = ZERO
    u_pre = ZERO

    hist : do
      ! Sample distance
      call rand % generateDistance(self % mfp, distance)

      event = LOST_EV
      do while (event /= COLL_EV)
        ! Save pre-movement state
        matIdx = coords % getMatIdx()
        uniqueId = coords % getUniqueId()
        maxDist = distance
        if (self % robust) then
          r_pre = coords % getPosition(1)
          u_pre = coords % getDirection(1)

        end if

        ! Move in geometry
        if (self % cache) then
          call self % move(coords, distance, event, cache_space)

        else
          call self % move(coords, distance, event)

        end if

        ! If robust verify matIdx in the mid point
        if (self % robust) then
          r = r_pre + u_pre * HALF * distance
          call self % whatIsAt(r, u_pre, unique_mid, mat_mid)

          if (matIdx /= mat_mid) then
            print *, "EVENT: ", event
            print *, "PRE MOVE: ", r_pre
            print *, "DIR", coords % getDirection(1)
            print *, "POST MOVE", coords % getPosition(1)
            print *, "MAT CHECKED AT:", r
            print *, "WITH DIRECTION:", u_pre
            print *, "CORRECT MAT:", mat_mid
            print *, "HAS MAT:", matIdx
            call fatalError(Here, 'Ray has lost correct material.')
          end if

        end if

        ! Score result
        !$omp atomic
        self % totDist = self % totDist + distance
        if (matIdx /= VOID_MAT) then
          !$omp atomic
          self % res(matIdx, SCORE) = self % res(matIdx, SCORE) + distance
        end if

        ! Set to remaining distance
        distance = maxDist - distance
      end do

      ! Kill the ray
      call rand % generate(randomNumber)
      if (randomNumber < self % abs_prob .or. coords % getMatIdx() == OUTSIDE_MAT) exit hist

      ! Scatter the ray
      call rand % generateMu(mu)
      call rand % generatePhi(phi)
      call coords % rotate(mu, phi)

    end do hist

  end subroutine trackRay

  !!
  !! Output calculation results to the console
  !!
  !! Convert cumulative sums to mean and absolute standard deviation and
  !! print them to the console.
  !!
  !! Args:
  !!   None
  !!
  subroutine printResults(self)
    class(rayVolPhysicsPackage), intent(in) :: self
    real(defReal)                           :: mean, SD, var
    real(defReal)                           :: V, V_SD
    integer(shortInt)                       :: i, nCycles

    ! Calculate speed and its SD
    nCycles = self % getCyclesNumber()
    V = self % ray_speed(1) / nCycles
    V_SD = self % ray_speed(2) / nCycles - V * V
    V_SD = ONE / (nCycles - 1) * sqrt(V_SD)

    print *
    print '(A, ES12.5, A, ES12.5)', " Ray speed [m/s]: ", V, " +/- ", V_SD
    print *, "RELATIVE VOLUME FOR MATERIALS: "
    do i = 1, mm_nMat()
      mean = self % res(i, CSUM) / nCycles
      var = self % res(i, CSUM2) / nCycles - mean * mean
      SD = ONE / (nCycles - 1) * sqrt(var)
      print '(A, A, A, ES12.5, A, ES12.5)', " Material: ", mm_matName(i), " Vol", mean, " +/-", SD
    end do

  end subroutine printResults

  !!
  !! Print settings of the Ray-tracking volume calculation
  !!
  !! Args:
  !!   None
  !!
  subroutine printSettings(self)
    class(rayVolPhysicsPackage), intent(in) :: self

    print *, repeat("<>", MAX_COL/2)
    print *, "/\/\ RAY-TRACING RELATIVE VOLUME CALCULATION /\/\"
    print *, "Total Cycles:    ", numToChar(self % getCyclesNumber())
    print *, "Rays per cycle: ", numToChar(self % getParticlesNumber())
    print *, "Initial RNG Seed:   ", numToChar(self % rand % getInitialSeed())
    print *
    print *, repeat("<>", MAX_COL/2)

  end subroutine printSettings

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(rayVolPhysicsPackage), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Clean contents
    !call self % rand % kill()

    self % mfp = ZERO
    self % abs_prob = ZERO
    self % robust = .false.

    if (allocated(self % res)) deallocate(self % res)
    self % totDist = ZERO
    self % ray_speed = ZERO

  end subroutine kill

end module rayVolPhysicsPackage_class
