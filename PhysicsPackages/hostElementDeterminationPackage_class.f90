module hostElementDeterminationPackage_class

  use coord_class,          only : coord
  use dictionary_class,     only : dictionary
  use geometry_inter,       only : geometry
  use geometryStd_class,    only : geometryStd
  use geometryFactory_func, only : new_geometry
  use geometryReg_mod,      only : geomIdx, geomPtr
  use hashFunctions_func,   only : FNV_1
  use meshUniverse_class,   only : meshUniverse
  use nuclearDataReg_mod,   only : ndReg_init => init
  use numPrecision
  use physicsPackage_inter, only : physicsPackage
  use rng_class,            only : rng
  use timer_mod,            only : registerTimer
  use universe_inter,       only : universe

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(physicsPackage) :: hostElementDeterminationPackage
    private
    class(geometry), pointer :: geom => null()
    integer(shortInt)        :: geometryIdx = 0, pop = 0, timerMain = 0
    real(defReal)            :: cpu_start_time = ZERO, cpu_end_time = ZERO
    type(RNG), pointer       :: pRNG => null()
  contains
    procedure :: init
    procedure :: run
  end type hostElementDeterminationPackage

contains
  !!
  !!
  !!
  subroutine init(self, dict)
    class(hostElementDeterminationPackage), intent(inout) :: self
    class(dictionary), intent(inout)                      :: dict
    character(8)                                          :: date
    character(10)                                         :: time
    character(nameLen)                                    :: geometryName
    character(:), allocatable                             :: string
    integer(longInt)                                      :: seed
    integer(shortInt)                                     :: seed_temp

    call cpu_time(self % cpu_start_time)

    ! Read calculation settings
    call dict % get(self % pop, 'pop')

    ! Register timer
    self % timerMain = registerTimer('transportTime')

    ! Initialise RNG
    allocate(self % pRNG)

    ! *** It is a bit silly but dictionary cannot store longInt for now
    !     so seeds are limited to 32 bits (can be -ve)
    if (dict % isPresent('seed')) then
      call dict % get(seed_temp, 'seed')

    else
      ! Obtain time string and hash it to obtain random seed
      call date_and_time(date, time)
      string = date//time
      call FNV_1(string, seed_temp)

    end if
    seed = seed_temp
    print *, 'Seed: ', seed
    call self % pRNG % init(seed)

    ! Build Nuclear Data
    call ndReg_init(dict % getDictPtr("nuclearData"))

    ! Build geometry
    geometryName = 'testGeometry'
    call new_geometry(dict % getDictPtr('geometry'), geometryName)
    self % geometryIdx = geomIdx(geometryName)
    self % geom => geomPtr(self % geometryIdx)

  end subroutine init

  !!
  !!
  !!
  subroutine run(self)
    class(hostElementDeterminationPackage), intent(inout) :: self
    type(meshUniverse), pointer                           :: meshUniversePtr
    class(universe), pointer                              :: universePtr
    integer(shortInt)                                     :: i, matIdx, uniqueId
    real(defReal)                                         :: t1, t2
    real(defReal), dimension(3)                           :: r, randomNumbers, u
    real(defReal), dimension(6)                           :: bounds
    type(coord)                                           :: coords
    type(geometryStd), pointer                            :: geometryStdPtr

    ! Downcast.
    select type(temp => self % geom)
      type is(geometryStd)
        geometryStdPtr => temp

    end select

    ! Get pointer to mesh universe then downcast.
    universePtr => geometryStdPtr % geom % unis % getPtr_fast(2)
    select type(temp => universePtr)
      type is(meshUniverse)
        meshUniversePtr => temp

    end select

    u = [ONE, ZERO, ZERO]
    call coords % setDirection(u)
    bounds = self % geom % bounds()
    call cpu_time(t1)
    do i = 1, self % pop
      ! Sample Position
      call self % pRNG % generate(randomNumbers)
      call coords % setPosition((bounds(4:6) - bounds(1:3)) * randomNumbers + bounds(1:3))

      ! Find element occupied by coordinates.
      call meshUniversePtr % mesh % ptr % findHostElement(coords)

    end do

    call cpu_time(t2)
    print*, "-------------------------------------------------------------"
    print*, "/\/\ Host element determination procedure time /\/\"
    print *, 'CPU time: ', t2 - t1
    print*, "-------------------------------------------------------------"

  end subroutine run

end module hostElementDeterminationPackage_class