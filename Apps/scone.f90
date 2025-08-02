program scone

  use commandLineUI,              only : addClOption, clOptionIsPresent, getFromCL, getInputFile, killCommandLineOptions => kill
  use dictionary_class,           only : dictionary
  use dictParser_func,            only : fileToDict
  use genericProcedures,          only : printStart
  use geometryFactory_func,       only : new_geometry
  use geometryReg_mod,            only : geomIdx, geomPtr, killGeometryRegistry => kill
  use nuclearDataReg_mod,         only : initNuclearDataRegistry => init, killNuclearDataRegistry => kill
  use numPrecision
  use openmp_func,                only : ompSetNumThreads
  use physicsPackage_inter,       only : initPhysicsPackagePayload, physicsPackage
  use physicsPackageFactory_func, only : new_physicsPackage
  use timer_mod,                  only : killTimer, registerTimer, secToChar, timerStart, timerStop, timerTime
  use vizPhysicsPackage_class,    only : vizPhysicsPackage

  implicit none
  type(dictionary), target           :: input
  class(physicsPackage), allocatable :: core
  character(:), allocatable          :: inputPath
  integer(shortInt)                  :: timerIdx
  integer(shortInt)                  :: cores
  type(dictionary), pointer          :: subDict
  character(nameLen)                 :: geometryName
  type(initPhysicsPackagePayload)    :: physicsPackagePayload

  ! Add command line options here
  call addClOption('--plot', 0, ['int'],&
          'Executes geometry plotting specified by a viz dict in the input file')
#ifdef _OPENMP
  call addClOption('--omp', 1, ['int'], &
          'Number of OpenMP threads in a parallel calculation')
#endif

  ! Get path to input file
  call getInputFile(inputPath)

  ! Set Number of threads
  cores = 1
  if (clOptionIsPresent('--omp')) call getFromCL(cores, '--omp', 1)
  call ompSetNumThreads(cores)

  ! Register timer and begin time tracking.
  timerIdx = registerTimer('Main Timer')
  call printStart()
  call timerStart(timerIdx)

  ! Parse file to dictionary.
  call fileToDict(input, inputPath)

  ! Build nuclear data and geometry.
  call initNuclearDataRegistry(input % getDictPtr('nuclearData'))
  
  subDict => input % getDictPtr('geometry')
  geometryName = 'mainGeometry'
  call new_geometry(subDict, geometryName)
  
  ! Assemble payload then initialise physics package.
  physicsPackagePayload % dict => input
  physicsPackagePayload % geometryIdx = geomIdx(geometryName)
  physicsPackagePayload % geometry => geomPtr(physicsPackagePayload % geometryIdx)
  if (clOptionIsPresent('--plot')) then
    allocate(vizPhysicsPackage :: core)
    call core % init(physicsPackagePayload)

  else
    allocate(core, source = new_physicsPackage(physicsPackagePayload))

  endif

  call core % run()

  call timerStop(timerIdx)
  print *, 'Total calculation time: ', trim(secToChar(timerTime(timerIdx)))
  print *, 'Have a good day and enjoy your result analysis!'

  ! Clean up.
  if (allocated(core)) then
    call core % kill()
    deallocate(core)

  end if

  call killCommandLineOptions()
  call killTimer()

  call input % kill()
  if (allocated(inputPath)) deallocate(inputPath)

  call killNuclearDataRegistry()
  call killGeometryRegistry()

end program scone
