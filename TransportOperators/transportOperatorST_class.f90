!!
!! Transport operator for surface tracking
!!
module transportOperatorST_class
  
  use dictionary_class,           only : dictionary
  use errors_mod,                 only : fatalError
  use geometry_inter,             only : geometry, distCache
  use nuclearDatabase_inter,      only : nuclearDatabase
  use numPrecision
  use physicalParticle_inter,     only : castPhysicalParticlePtr, physicalParticle
  use tallyAdmin_class,           only : tallyAdmin
  use tallyCodes
  use transportObject_inter,      only : transportObject
  use transportOperator_inter,    only : transportOperator, init_super => init
  use universalVariables

  implicit none
  private

  !!
  !! Transport operator that moves a particle with surface tracking
  !!
  !! Sample Input Dictionary:
  !!   trans { type transportOperatorST; cache 0;}
  !!
  type, public, extends(transportOperator) :: transportOperatorST
    logical(defBool)  :: cache = .true.
  contains
    procedure :: transit => surfaceTracking
    ! Override procedure
    procedure :: init
  end type transportOperatorST

contains

  !!
  !! Performs surface tracking until a collision point is found
  !!
  subroutine surfaceTracking(self, object, tally)
    class(transportOperatorST), intent(inout) :: self
    class(transportObject), intent(inout)     :: object
    type(tallyAdmin), intent(inout)           :: tally
    class(physicalParticle), pointer          :: p
    integer(shortInt)                         :: event, materialIdx
    real(defReal)                             :: inverseSigmaT, distance
    type(distCache)                           :: cache
    character(*), parameter                   :: here = 'surfaceTracking (transportOperatorST_class.f90)'

    ! Downcast transport object to physical particle.
    p => castPhysicalParticlePtr(object, .true.)

    STLoop: do
      ! Obtain the local cross-section
      materialIdx = p % getMaterialIdx()
      if (materialIdx == VOID_MAT) then
        distance = INF

      else
        inverseSigmaT = ONE / self % getTrackingXS(p, materialIdx, MATERIAL_XS)
        call p % generateDistance(inverseSigmaT, distance)

        ! Should never happen! Catches NaN distances
        if (distance /= distance) call fatalError(Here, "Distance is NaN")

      end if

      ! Save state before movement
      call p % savePrePathState()

      ! Move to the next stop.
      if (self % cache) then
        call self % move(p, distance, event, cache)

      else
        call self % move(p, distance, event)

      end if

      ! Send tally report for a path moved
      call tally % reportPath(p, distance)

      ! Kill particle if it has leaked.
      materialIdx = p % getMaterialIdx()
      if (materialIdx == OUTSIDE_FILL) then
        call p % setFate(LEAK_FATE)
        call p % setIsDead(.true.)
        
      end if

      ! Give error if the particle somehow ended in an undefined material
      if (materialIdx == UNDEF_MAT) then
        print *, p % getGlobalPosition()
        call fatalError(here, 'Particle is in undefined material')

      end if

      ! Return if particle stoped at collision (not cell boundary)
      if (event == COLL_EV .or. p % getIsDead()) exit STLoop

    end do STLoop

    call tally % reportTrans(p)

  end subroutine surfaceTracking

  !!
  !! Initialise ST operator from a dictionary
  !!
  !! See transportOperator_inter for details
  !!
  subroutine init(self, dict)
    class(transportOperatorST), intent(inout) :: self
    class(dictionary), intent(in)             :: dict

    ! Initialise superclass
    call init_super(self, dict)
    if (dict % isPresent('cache')) call dict % get(self % cache, 'cache')

  end subroutine init

end module transportOperatorST_class