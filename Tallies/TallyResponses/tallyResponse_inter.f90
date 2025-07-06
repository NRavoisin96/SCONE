module tallyResponse_inter

  use numPrecision
  use dictionary_class,      only : dictionary
  use particle_class,        only : particle

  ! Nuclear Data interface
  use nuclearDatabase_inter, only : nuclearDatabase

  implicit none
  private


  !!
  !! Abstract interface for all tallyResponses
  !!
  !! Very simple class, which given a particle returns a real number
  !! Real number is used to weight FLUX sample to score reaction rates etc.
  !! Returns only scalar to move all logic for dealing with multiple responces to tallyClerks.
  !! Thus tallyResponses should be quick and easy to write
  !!
  !! Interface:
  !!   init -> Initialise
  !!   get  -> Get velue of the response
  !!   kill -> Return to uninitialised state
  !!
  type, public,abstract :: tallyResponse
    private

  contains
    procedure(init), deferred :: init
    procedure(get), deferred  :: get
    procedure(kill), deferred :: kill

  end type tallyResponse

  abstract interface

    !!
    !! Initialise Response from dictionary
    !!
    !! Args:
    !!   dict [in] -> DIctionary with the data
    !!
    !! Errors:
    !!   Depend on specific implementation.
    !!   fatalError if there is a mistake in definition
    !!
    subroutine init(self, dict)
      import :: tallyResponse, &
                dictionary
      class(tallyResponse), intent(inout) :: self
      class(dictionary), intent(in)       :: dict

    end subroutine init

    !!
    !! Get value of response
    !!
    !! Args:
    !!   p [in]         -> Particle to provide state
    !!   xsData [inout] -> Nuclear Database used by the particle
    !!
    !! Result:
    !!   Value of the response for particle p
    !!
    !! Errors:
    !!   Depend on specific implementation
    !!
    subroutine get(self, p, val, xsData)
      import :: defReal, nuclearDatabase, particle, tallyResponse
      class(tallyResponse), intent(in)                :: self
      class(particle), intent(in)                     :: p
      real(defReal), intent(out)                      :: val
      class(nuclearDatabase), intent(inout), optional :: xsData
    end subroutine get

    !!
    !! Return to uninitialised state
    !!
    !! Args:
    !!   None
    !!
    !! Errors:
    !!   None
    !!
    elemental subroutine kill(self)
      import :: tallyResponse
      class(tallyResponse), intent(inout) :: self
    end subroutine kill

  end interface

end module tallyResponse_inter
