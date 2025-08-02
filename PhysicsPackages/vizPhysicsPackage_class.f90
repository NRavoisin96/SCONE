module vizPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError
  use dictionary_class,               only : dictionary

  ! Timers
  use timer_mod,                      only : registerTimer, timerStart, timerStop, &
                                             timerTime, timerReset, secToChar

  ! Physics package interface
  use physicsPackage_inter,           only : initPhysicsPackagePayload, physicsPackage

  ! Geometry
  use geometry_inter,                 only : geometry
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_geomIdx  => geomIdx
  use geometryFactory_func,           only : new_geometry

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat           => nMat
  use nuclearDataReg_mod,             only : ndReg_init        => init ,&
                                             ndReg_getMatNames => getMatNames
  use nuclearDatabase_inter,          only : nuclearDatabase

  ! Visualisation
  use visualiser_class,               only : visualiser

  implicit none
  private

  !!
  !! Physics Package for eigenvalue calculations
  !!
  type, public, extends(physicsPackage) :: vizPhysicsPackage
    private
    ! Building blocks
    class(geometry), pointer :: geom => null()
    integer(shortInt)        :: geomIdx = 0
    type(visualiser)         :: viz

    ! Timer bins
    integer(shortInt) :: timerMain

  contains
    procedure :: init
    procedure :: run
    procedure :: kill

  end type vizPhysicsPackage

contains

  !!
  !! Calls visualiser to generate visualisation
  !!
  subroutine run(self)
    class(vizPhysicsPackage), intent(inout) :: self

    print *, "Constructing visualisation"
    call self % viz % makeViz()
    call self % viz % kill()

  end subroutine

  !!
  !! Initialise from individual components and dictionaries
  !!
  subroutine init(self, payload)
    class(vizPhysicsPackage), intent(inout)     :: self
    type(initPhysicsPackagePayload), intent(in) :: payload
    character(*), parameter                     :: Here = 'init (vizPhysicsPackage_class.f90)'

    ! Register timer
    self % timerMain = registerTimer('transportTime')

    ! Build geometry
    self % geomIdx = payload % geometryIdx
    self % geom => payload % geometry

    ! Call visualisation
    if (payload % dict % isPresent('viz')) then
      print *, "Initialising visualiser"
      call self % viz % init(self % geom, payload % dict % getDictPtr('viz'))

    else
      call fatalError(here,'Must provide viz dict for plotting.')

    endif

  end subroutine init

  !!
  !! Deallocate memory
  !!
  subroutine kill(self)
    class(vizPhysicsPackage), intent(inout) :: self

    ! TODO: This subroutine

  end subroutine kill

end module vizPhysicsPackage_class
