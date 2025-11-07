module source_inter

  use dictionary_class,           only : dictionary
  use geometry_inter,             only : geometry
  use numPrecision
  use particleDungeon_class,      only : particleDungeon
  use RNG_class,                  only : RNG
  use transportObjectState_class, only : transportObjectState

  implicit none
  private

  !!
  !! Extendable scource class procedures
  !!
  public :: init, kill

  !!
  !! Abstract interface of source for particles
  !!
  !! Source generates particles from specified distributions
  !! for, e.g., fixed source calcs. or to generate initial
  !! distribution for eigenvalue calcs
  !!
  !! Public members:
  !!   geom -> Pointer to the geometry to ensure source is inside and
  !!           for more complicated source distribution
  !!
  !! Interface:
  !!   init              -> initialise the source
  !!   generate          -> generate particles to fill a dungeon
  !!   sampleParticle    -> sample particles from the corresponding distributions
  !!   kill              -> clean up the source
  !!
  type, public,abstract :: source
    private
    class(geometry), pointer            :: geom => null()
  contains
    procedure, non_overridable          :: generate
    procedure                           :: getGeometryPtr
    procedure                           :: init
    procedure                           :: kill
    procedure(sampleParticle), deferred :: sampleParticle
    procedure(sampleState), deferred    :: sampleState
  end type source

  abstract interface
    !!
    !! Sample particle's phase space co-ordinates
    !!
    !! Generates a phase-space state for a single particle
    !!
    !! Args:
    !!   p [inout] -> particle to be over-written
    !!
    !! Result:
    !!   A particle sampled the prescribed source
    !!
    subroutine sampleParticle(self, rand, state)
      import                                                :: RNG, source, transportObjectState
      class(source), intent(inout)                          :: self
      class(RNG), intent(inout)                             :: rand
      class(transportObjectState), allocatable, intent(out) :: state
    end subroutine sampleParticle

    !!
    !!
    !!
    subroutine sampleState(self, state)
      import                                                :: source, transportObjectState
      class(source), intent(in)                             :: self
      class(transportObjectState), allocatable, intent(out) :: state
    end subroutine sampleState

  end interface

contains
  !!
  !! Generate particles to populate a particleDungeon
  !!
  !! Fills a particle dungeon with n particles, sampled
  !! from the corresponding source distributions
  !!
  !! Args:
  !!   dungeon [inout] -> particle dungeon to be populated
  !!   n [in]          -> number of particles to place in dungeon
  !!
  !! Result:
  !!   A dungeon populated with n particles sampled from the source
  !!
  subroutine generate(self, dungeon, n, rand)
    class(source), intent(inout)             :: self
    type(particleDungeon), intent(inout)     :: dungeon
    integer(shortInt), intent(in)            :: n
    class(RNG), intent(inout)                :: rand
    class(transportObjectState), allocatable :: state
    type(RNG)                                :: pRand
    integer(shortInt)                        :: i
    
    ! Set dungeon size to begin.
    call dungeon % setSize(n)

    ! Generate n particles to populate dungeon
    !$omp parallel do private(pRand, state)
    do i = 1, n
      pRand = rand
      call pRand % stride(i)
      call self % sampleParticle(pRand, state)
      call dungeon % replace(state, i)

    end do
    !$omp end parallel do

    ! Advance RNG after particle generation.
    call rand % stride(n)

  end subroutine generate

  !!
  !!
  !!
  function getGeometryPtr(self) result(ptr)
    class(source), intent(in) :: self
    class(geometry), pointer  :: ptr

    ptr => self % geom

  end function getGeometryPtr

  !!
  !! Initialise source from dictionary & geometry
  !!
  !! Args:
  !!   dict [in] -> dict containing point source information
  !!   geom [in] -> pointer to a geometry
  !!
  subroutine init(self, dict, geom)
    class(source), intent(inout)         :: self
    class(dictionary), intent(in)        :: dict
    class(geometry), pointer, intent(in) :: geom

    self % geom => geom

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(source), intent(inout) :: self

    self % geom => null()

  end subroutine kill

end module source_inter