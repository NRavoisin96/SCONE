module physicalParticleFactory_func

  use CENeutron_class,             only : CENeutron
  use CENeutronState_class,        only : CENeutronState
  use CEParticleState_class,       only : CEParticleState
  use errors_mod,                  only : fatalError
  use MGNeutron_class,             only : MGNeutron
  use MGNeutronState_class,        only : MGNeutronState
  use MGParticleState_class,       only : MGParticleState
  use numPrecision
  use physicalParticle_inter,      only : physicalParticle
  use physicalParticleState_class, only : physicalParticleState
  use testCEParticle_class,        only : testCEParticle
  use testMGParticle_class,        only : testMGParticle
  use testPhysicalParticle_class,  only : testPhysicalParticle

  implicit none
  private

  ! Public procedures.
  public :: new_physicalParticle

  character(nameLen), dimension(*), parameter :: AVAILABLE_PHYSICAL_PARTICLE_STATES = ['CENeutronState       ', &
                                                                                       'CEParticleState      ', &
                                                                                       'CEPhotonState        ', &
                                                                                       'MGNeutronState       ', &
                                                                                       'MGParticleState      ', &
                                                                                       'MGPhotonState        ', &
                                                                                       'physicalParticleState']

contains
  !!
  !!
  !!
  function new_physicalParticle(state) result(new)
    class(physicalParticleState), intent(in) :: state
    class(physicalParticle), allocatable     :: new
    character(*), parameter                  :: HERE = 'new_physicalParticle (physicalParticleFactory_func.f90)'

    ! Allocate appropriate physical particle depending on the input state.
    select type(state)
      type is(CENeutronState)
        allocate(CENeutron :: new)

      type is(CEParticleState)
        allocate(testCEParticle :: new)

      type is(MGNeutronState)
        allocate(MGNeutron :: new)

      type is(MGParticleState)
        allocate(testMGParticle :: new)

      type is(physicalParticleState)
        allocate(testPhysicalParticle :: new)

      class default
        print *, AVAILABLE_PHYSICAL_PARTICLE_STATES
        call fatalError(HERE, 'Unrecognised physical particle state type.')

    end select

    ! Initialise particle from state.
    call new % init(state)

  end function new_physicalParticle

end module physicalParticleFactory_func