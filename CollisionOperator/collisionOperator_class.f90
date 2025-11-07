module collisionOperator_class

  use collisionProcessor_inter,       only : collisionProcessor
  use collisionProcessorFactory_func, only : new_collisionProcessor
  use dictionary_class,               only : dictionary
  use errors_mod,                     only : fatalError
  use genericProcedures,              only : numToChar
  use intMap_class,                   only : intMap
  use numPrecision
  use particleDungeon_class,          only : particleDungeon
  use physicalParticle_inter,         only : castPhysicalParticlePtr, physicalParticle
  use tallyAdmin_class,               only : tallyAdmin
  use transportObject_inter,          only : transportObject
  use universalVariables,             only : NOT_PRESENT

  implicit none
  private

  !!
  !! Local helper type to store polymorphic collisionProcessors in an array
  !!
  type, private :: collProc
    class(collisionProcessor), allocatable :: proc
  end type collProc

  !!
  !! Scalar collision operator
  !!  -> Maps particles of diffrent types to approperiate implementation
  !!     of collision physics (collisionProcessor).
  !!  -> Uses lookup table for to quickly map diffrent combination of physical(neutron, photon ...)
  !!     and processing type (CE, MG) to approperiate physics.
  !!  -> Can store up to 3 diffrent physics types
  !!  -> Gives fatal error if particle type is not recognised or not-supported
  !!
  !! Sample dictionary input( provisional will change):
  !!  collOpName {
  !!    #neutronCE {<collisonProcessorDefinition>} #
  !!    #neutronMG {<collisonProcessorDefinition>} #
  !!  }
  !!
  type, public :: collisionOperator
    private
    type(collProc), dimension(:), allocatable :: collisionProcessors
    type(intMap)                              :: particleTypeToCollisionProcessorMap
  contains
    ! Build procedures
    procedure :: init
    procedure :: kill
    ! Use procedures
    procedure :: collide
  end type collisionOperator

contains
  !!
  !! Initialise collision operator
  !!
  subroutine init(self, dict)
    class(collisionOperator), intent(inout)       :: self
    class(dictionary), intent(in)                 :: dict
    character(nameLen), dimension(:), allocatable :: processorNames
    character(nameLen)                            :: processorName
    integer(shortInt)                             :: i, nProcessors, particleType
    character(*), parameter                       :: here = 'init (collisionOperator_class.f90)'

    ! First retrieve all the processors' definitions from the dictionary and compute the number of processors
    ! to be generated.
    call dict % keys(processorNames, 'dict')
    nProcessors = size(processorNames)

    ! Allocate memory and initialise intMap.
    allocate(self % collisionProcessors(nProcessors))
    call self % particleTypeToCollisionProcessorMap % init(nProcessors)

    ! Loop over all collision processors and build them.
    do i = 1, nProcessors
      call new_collisionProcessor(dict % getDictPtr(processorNames(i)), self % collisionProcessors(i) % proc, &
                                  particleType)

      ! Add the particle type to the map. Return an error if a processor for the current particle type already exists.
      if (self % particleTypeToCollisionProcessorMap % getOrDefault(particleType, NOT_PRESENT) /= NOT_PRESENT) &
      call fatalError(here, 'A collision processor already exists for this type of particle.')
      call self % particleTypeToCollisionProcessorMap % add(particleType, i)

    end do

  end subroutine init

  !!
  !! Clear collision operator. Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(collisionOperator), intent(inout) :: self
    integer(shortInt)                       :: i

    ! Deallocate collision processors.
    if (allocated(self % collisionProcessors)) then
      do i = 1, size(self % collisionProcessors)
        if (allocated(self % collisionProcessors(i) % proc)) deallocate(self % collisionProcessors(i) % proc)

      end do
      deallocate(self % collisionProcessors)

    end if

    ! Kill map.
    call self % particleTypeToCollisionProcessorMap % kill()

  end subroutine kill

  !!
  !! Determine type of the particle and call approperiate collisionProcessor
  !!
  subroutine collide(self, object, tally, thisCycle, nextCycle)
    class(collisionOperator), intent(inout) :: self
    class(transportObject), intent(inout)   :: object
    type(tallyAdmin), intent(inout)         :: tally
    class(particleDungeon), intent(inout)   :: thisCycle, nextCycle
    class(physicalParticle), pointer        :: p
    integer(shortInt)                       :: idx
    character(*), parameter                 :: here = 'collide (collisionOperator_class.f90)'

    ! Downcast transportObject to physicalParticle.
    p => castPhysicalParticlePtr(object, .true.)

    ! Get index from map.
    idx = self % particleTypeToCollisionProcessorMap % getOrDefault(p % getType(), NOT_PRESENT)
    if (idx == NOT_PRESENT) call fatalError(here, 'Physics are not defined for this type of particle.')

    ! Call physics.
    call self % collisionProcessors(idx) % proc % collide(p, tally, thisCycle, nextCycle)

  end subroutine collide

end module collisionOperator_class