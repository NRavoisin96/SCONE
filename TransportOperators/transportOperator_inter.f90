module transportOperator_inter

  use dictionary_class,       only : dictionary
  use errors_mod,             only : fatalError
  use geometry_inter,         only : geometry, distCache
  use geometryReg_mod,        only : gr_geomPtr => geomPtr
  use nuclearDatabase_inter,  only : nuclearDatabase
  use nuclearDataReg_mod,     only : ndReg_get => get
  use numPrecision
  use tallyAdmin_class,       only : tallyAdmin
  use transportObject_inter,  only : transportObject
  use universalVariables

  implicit none
  private

  ! Public procedures.
  public :: init, kill

  !!
  !! This is an abstract interface for all types of transport processing
  !!   -> This interface only deals with scalar processing of particle transport
  !!   -> Assumes that particle moves without any external forces (assumes that particle
  !!      moves along straight lines between colisions)
  !!
  !! Public interface:
  !!   transport(p, tally, thisCycle, nextCycle) -> given particle, tally and particle dungeons
  !!     for particles in this and next cycle performs movement of a particle in the geometry.
  !!     Sends transition report to the tally. Sends history report as well if particle dies.
  !!   init(dict, geom) -> initialises transport operator from a dictionary and pointer to a
  !!                       geometry
  !!
  !! Customisable procedures or transport actions
  !!   transit(p, tally, thisCycle, nextCycle) -> implements movement from collision to collision
  !!
  type, abstract, public :: transportOperator
    private
    class(nuclearDatabase), pointer :: xsData => null() ! Nuclear Data block pointer -> public so it can be used by subclasses (protected member)
    class(geometry), pointer        :: geom => null()   ! Geometry pointer -> public so it can be used by subclasses (protected member)
  contains
    ! Public interface
    procedure, non_overridable   :: transport
    ! Extentable initialisation and deconstruction procedure
    procedure                    :: getTrackingXS
    procedure                    :: getTrackMatXS
    procedure                    :: init
    procedure                    :: kill
    procedure                    :: move
    procedure                    :: teleport
    ! Customisable deferred procedures
    procedure(transit), deferred :: transit
  end type transportOperator

  abstract interface
    !!
    !! Move particle from collision to collision
    !!  Kill particle if needed
    !!
    subroutine transit(self, object, tally)
      import                                  :: tallyAdmin, transportObject, transportOperator
      class(transportOperator), intent(inout) :: self
      class(transportObject), intent(inout)   :: object
      type(tallyAdmin), intent(inout)         :: tally
    end subroutine transit

  end interface

contains
  !!
  !!
  !!
  function getTrackingXS(self, object, materialIdx, what) result(trackingXS)
    class(transportOperator), intent(in) :: self
    class(transportObject), intent(in)   :: object
    integer(shortInt), intent(in)        :: materialIdx, what
    real(defReal)                        :: trackingXS

    trackingXS = self % xsData % getTrackingXS(object, materialIdx, what)

  end function getTrackingXS

  !!
  !!
  !!
  function getTrackMatXS(self, object, materialIdx) result(trackMatXS)
    class(transportOperator), intent(in) :: self
    class(transportObject), intent(in)   :: object
    integer(shortInt), intent(in)        :: materialIdx
    real(defReal)                        :: trackMatXS

    trackMatXS = self % xsData % getTrackMatXS(object, materialIdx)

  end function getTrackMatXS

  !!
  !! Initialise transport operator from dictionary and geometry
  !!
  subroutine init(self, dict)
    class(transportOperator), intent(inout)  :: self
    class(dictionary), intent(in)            :: dict

    ! Do nothing

  end subroutine init

  !!
  !! Free memory. Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(transportOperator), intent(inout) :: self

    self % geom => null()
    self % xsData => null()

  end subroutine kill

  !!
  !!
  !!
  subroutine move(self, object, distance, event, cache)
    class(transportOperator), intent(in)     :: self
    class(transportObject), intent(in)       :: object
    real(defReal), intent(inout)             :: distance
    integer(shortInt), intent(out)           :: event
    type(distCache), intent(inout), optional :: cache

    call self % geom % move(object % getCoordsPtr(), distance, event, cache)

  end subroutine move

  !!
  !!
  !!
  subroutine teleport(self, object, distance)
    class(transportOperator), intent(in)  :: self
    class(transportObject), intent(inout) :: object
    real(defReal), intent(in)             :: distance

    call self % geom % teleport(object % getCoordsPtr(), distance)

  end subroutine teleport

  !!
  !! Master non-overridable subroutine to perform transport
  !!  Performs everything common to all types of transport
  !!
  subroutine transport(self, object, tally)
    class(transportOperator), intent(inout) :: self
    class(transportObject), intent(inout)   :: object
    type(tallyAdmin), intent(inout)         :: tally
    character(*), parameter                 :: Here = 'transport (transportOperator_inter.f90)'

    ! Get nuclear data pointer form the particle
    self % xsData => ndReg_get(object % getType())

    ! Save geometry pointer
    self % geom => gr_geomPtr(object % getGeometryIdx())

    ! Save pre-transition state
    call object % savePreTransitionState()

    ! Perform transit
    call self % transit(object, tally)

    ! Send history reports if particle died
    if (object % getIsDead()) call tally % reportHist(object)

  end subroutine transport

end module transportOperator_inter