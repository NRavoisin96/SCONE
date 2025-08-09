module physicsPackage_inter

  use coordList_class,    only : coordList
  use dictionary_class,   only : dictionary
  use fieldFactory_func,  only : new_field
  use geometry_inter,     only : distCache, geometry
  use hashFunctions_func, only : FNV_1
  use numPrecision
  use outputFile_class,   only : outputFile
  use timer_mod,          only : registerTimer
  use visualiser_class,   only : visualiser
  use universalVariables

  implicit none
  private

  ! Public procedures.
  public :: copyPayload, init, kill

  !!
  !! Abstract interface of physics Package
  !! Physics package is controles a calculation flow
  !! Each type of calculation has diffrent physics package
  !!
  type, public,abstract      :: physicsPackage
    private
    integer(shortInt)        :: geomIdx = 0, N_cycles = 0, pop = 0, timerMain = 0
    integer(longInt)         :: seed = 0_longInt
    real(defReal)            :: cpu_time_end = ZERO, cpu_time_start = ZERO
    character(pathLen)       :: outputFile = ''
    character(nameLen)       :: outputFormat = ''
    class(geometry), pointer :: geom => null()
  contains
    procedure                                   :: buildVisualisation
    procedure                                   :: collectResults
    procedure(collectSpecificResults), deferred :: collectSpecificResults
    procedure                                   :: getCyclesNumber
    procedure                                   :: getGeometryBounds
    procedure                                   :: getGeometryIdx
    procedure                                   :: getInitialSeed
    procedure                                   :: getOutputFile
    procedure                                   :: getParticlesNumber
    procedure                                   :: getTimerMain
    procedure                                   :: init
    procedure                                   :: kill
    procedure                                   :: move
    procedure                                   :: placeCoord
    procedure(run), deferred                    :: run
    procedure                                   :: whatIsAt
  end type physicsPackage

  !!
  !!
  !!
  type, public :: initPhysicsPackagePayload
    type(dictionary), pointer :: dict => null()
    integer(shortInt)         :: geometryIdx = 0
    class(geometry), pointer  :: geometry => null()
  end type initPhysicsPackagePayload

  abstract interface
    !!
    !!
    !!
    subroutine collectSpecificResults(self, out)
      import                            :: outputFile, physicsPackage
      class(physicsPackage), intent(in) :: self
      type(outputFile), intent(inout)   :: out
    end subroutine collectSpecificResults

    !!
    !! Run calculation in the physics package
    !!
    subroutine run(self)
      import                               :: physicsPackage
      class(physicsPackage), intent(inout) :: self
    end subroutine run

  end interface

contains
  !!
  !!
  !!
  subroutine buildVisualisation(self, dict)
    class(physicsPackage), intent(inout) :: self
    type(dictionary), intent(in)         :: dict
    type(visualiser)                     :: visualisation

    print *, 'Initialising visualiser'
    call visualisation % init(self % geom, dict)
    print *, 'Constructing visualisation'
    call visualisation % makeViz()
    call visualisation % kill()

  end subroutine buildVisualisation

  !!
  !!
  !!
  subroutine collectResults(self)
    class(physicsPackage), intent(inout) :: self
    type(outputFile)                     :: out
    character(nameLen)                   :: name

    call out % init(self % outputFormat, filename = self % outputFile)

    call cpu_time(self % cpu_time_end)
    name = 'Total_CPU_Time'
    call out % printValue(self % cpu_time_end - self % cpu_time_start, name)

    name = 'seed'
    call out % printValue(self % seed, name)

    name = 'pop'
    call out % printValue(self % pop, name)

    ! Print class specific results.
    call self % collectSpecificResults(out)

  end subroutine collectResults

  !!
  !!
  !!
  subroutine copyPayload(inputPayload, outputPayload)
    class(initPhysicsPackagePayload), intent(in)    :: inputPayload
    class(initPhysicsPackagePayload), intent(inout) :: outputPayload

    outputPayload % dict => inputPayload % dict
    outputPayload % geometryIdx = inputPayload % geometryIdx
    outputPayload % geometry => inputPayload % geometry

  end subroutine copyPayload

  !!
  !!
  !!
  elemental function getCyclesNumber(self) result(nCycles)
    class(physicsPackage), intent(in) :: self
    integer(shortInt)                 :: nCycles

    nCycles = self % N_cycles

  end function getCyclesNumber

  !!
  !!
  !!
  function getGeometryBounds(self) result(bounds)
    class(physicsPackage), intent(in) :: self
    real(defReal), dimension(6)       :: bounds

    bounds = self % geom % bounds()

  end function getGeometryBounds

  !!
  !!
  !!
  elemental function getGeometryIdx(self) result(geometryIdx)
    class(physicsPackage), intent(in) :: self
    integer(shortInt)                 :: geometryIdx

    geometryIdx = self % geomIdx

  end function getGeometryIdx

  !!
  !!
  !!
  elemental function getInitialSeed(self) result(seed)
    class(physicsPackage), intent(in) :: self
    integer(longInt)                  :: seed

    seed = self % seed

  end function getInitialSeed

  !!
  !!
  !!
  elemental function getOutputFile(self) result(outputFile)
    class(physicsPackage), intent(in) :: self
    character(pathLen)                :: outputFile

    outputFile = self % outputFile

  end function getOutputFile

  !!
  !!
  !!
  elemental function getParticlesNumber(self) result(nParticles)
    class(physicsPackage), intent(in) :: self
    integer(shortInt)                 :: nParticles

    nParticles = self % pop

  end function getParticlesNumber

  !!
  !!
  !!
  elemental function getTimerMain(self) result(timerMain)
    class(physicsPackage), intent(in) :: self
    integer(shortInt)                 :: timerMain

    timerMain = self % timerMain

  end function getTimerMain

  !!
  !! Initialise Physics Package from dictionary
  !!
  subroutine init(self, payload)
    class(physicsPackage), intent(inout)          :: self
    class(initPhysicsPackagePayload), intent(in)  :: payload
    character(8)                                  :: date
    character(10)                                 :: time
    character(nameLen)                            :: fieldName
    character(:), allocatable                     :: trimmedFieldName
    character(nameLen), dimension(:), allocatable :: fieldNames
    class(dictionary), pointer                    :: fieldsDict
    integer(shortInt)                             :: i
    type(outputFile)                              :: test_out

    ! Initialise CPU time.
    call cpu_time(self % CPU_time_start)

    ! Register timer
    self % timerMain = registerTimer('transportTime')

    ! Read outputfile path
    call payload % dict % getOrDefault(self % outputFile, 'outputFile', './output')

    ! Get output format and verify
    ! Initialise output file before calculation (so mistake in format will be caught early)
    call payload % dict % getOrDefault(self % outputFormat, 'outputFormat', 'asciiMATLAB')
    call test_out % init(self % outputFormat)

    ! Load number of cycles and population from dictionary.
    call payload % dict % get(self % N_cycles, 'cycles')
    call payload % dict % get(self % pop, 'pop')

    ! Assign geomIdx and geom pointer.
    self % geomIdx = payload % geometryIdx
    self % geom => payload % geometry

    ! Get seed from dictionary or create it from current date and time.
    if (payload % dict % isPresent('seed')) then
      call payload % dict % get(self % seed, 'seed')

    else
      call date_and_time(date, time)
      call FNV_1(date // time, self % seed)

    end if

    ! Build fields.
    if (payload % dict % isPresent('fields')) then
      fieldsDict => payload % dict % getDictPtr('fields')
      call fieldsDict % keys(fieldNames, 'dict')
      
      do i = 1, size(fieldNames)
        trimmedFieldName = trim(fieldNames(i))
        select case(trimmedFieldName)
          case('T', 'temp', 'Temp', 'temperature', 'Temperature')
            fieldName = nameTemperature

          case default
            fieldName = fieldNames(i)

        end select
        call new_field(fieldsDict % getDictPtr(trimmedFieldName), fieldName)

      end do

    end if

  end subroutine init

  !!
  !!
  !!
  subroutine move(self, coords, d, event, cache)
    class(physicsPackage), intent(in)        :: self
    type(coordList), intent(inout)           :: coords
    real(defReal), intent(out)               :: d
    integer(shortInt), intent(out)           :: event
    type(distCache), intent(inout), optional :: cache

    call self % geom % move(coords, d, event, cache)

  end subroutine move

  !!
  !!
  !!
  subroutine kill(self)
    class(physicsPackage), intent(inout) :: self

    ! Local.
    self % geomIdx = 0
    self % N_cycles = 0
    self % pop = 0
    self % timerMain = 0
    self % seed = 0_longInt
    self % cpu_time_end = ZERO
    self % cpu_time_start = ZERO
    self % outputFile = ''
    self % outputFormat = ''
    self % geom => null()

  end subroutine kill

  !!
  !!
  !!
  subroutine placeCoord(self, coords)
    class(physicsPackage), intent(in) :: self
    type(coordList), intent(inout)    :: coords

    call self % geom % placeCoord(coords)

  end subroutine placeCoord

  !!
  !!
  !!
  subroutine whatIsAt(self, r, u, uniqueId, matIdx)
    class(physicsPackage), intent(in)       :: self
    real(defReal), dimension(3), intent(in) :: r, u
    integer(shortInt), intent(inout)        :: uniqueId
    integer(shortInt), intent(out)          :: matIdx

    call self % geom % whatIsAt(matIdx, uniqueId, r, u)

  end subroutine whatIsAt

end module physicsPackage_inter