module configSource_inter

  use numPrecision
  use RNG_class,                  only : RNG
  use source_inter,               only : source,  kill_super => kill
  use transportObjectState_class, only : transportObjectState
  
  implicit none
  private

  !!
  !! Extendable scource class procedures
  !!
  public :: kill

  !!
  !! Configurable source
  !!
  !! Generates a single particle sample by calling a number of subroutines
  !! related to each component
  !!
  !! A single sample of a particle state is created by calling `sample****` procedures in order
  !! given in `sampleParticle` function. Note that `sampleEnergyAngle` is called AFTER
  !! `sampleEnergy`.
  !!
  !! Interface:
  !!   source_inter Interface
  !!   sampleType        -> sets the particle type
  !!   samplePosition    -> samples the particle's position in the geometry
  !!   sampleEnergy      -> samples the particle's energy
  !!   sampleEnergyAngle -> samples the particle's energy and angle from corresponding distr.
  !!
  type, public, abstract, extends(source) :: configSource

  contains
    procedure                              :: kill
    procedure(sampleEnergy), deferred      :: sampleEnergy
    procedure(sampleEnergyAngle), deferred :: sampleEnergyAngle
    procedure(samplePosition), deferred    :: samplePosition
    procedure                              :: sampleState
    procedure(sampleType), deferred        :: sampleType
  end type configSource

  abstract interface
    !!
    !! Sample particle Energy/Group
    !!
    !! Sets energy of a particle to a CE value or a MG index
    !! Also sets 'isMG' flag to .true. or .false.
    !!
    !! Inputs:
    !!   p [inout] -> particleState to be given a position
    !!   rand [in] -> random number generator
    !!
    subroutine sampleEnergy(self, state, rand)
      import                                     :: configSource, RNG, transportObjectState
      class(configSource), intent(inout)         :: self
      class(transportObjectState), intent(inout) :: state
      type(RNG), intent(inout)                  :: rand
    end subroutine sampleEnergy

    !!
    !! Sample particle Energy/Group and angle Angle
    !!
    !! Sets diraction of a particle together with its energy.
    !! Sampling of energy is optional if Angle & Energy are uncorrelated
    !! Is called after `sampleEnergy`, to overwrite value provided by that subroutine
    !!
    !! Inputs:
    !!   p [inout] -> particleState to be given a position
    !!   rand [in] -> random number generator
    !!
    subroutine sampleEnergyAngle(self, state, rand)
      import                                     :: configSource, RNG, transportObjectState
      class(configSource), intent(inout)         :: self
      class(transportObjectState), intent(inout) :: state
      type(RNG), intent(inout)                  :: rand
    end subroutine sampleEnergyAngle

    !!
    !! Sample particle position
    !!
    !! Sets position of the particle p
    !!
    !! Inputs:
    !!   p [inout] -> particleState to be given a position
    !!   rand [in] -> random number generator
    !!
    subroutine samplePosition(self, state, rand)
      import                                     :: configSource, RNG, transportObjectState
      class(configSource), intent(inout)         :: self
      class(transportObjectState), intent(inout) :: state
      type(RNG), intent(inout)                  :: rand
    end subroutine samplePosition

    !!
    !! Sample Type of a particle
    !!
    !! Sets 'Type' in the particleState p (e.g. P_NEUTRON)
    !!
    !! Inputs:
    !!   p [inout] -> particleState to be given a type
    !!   rand [in] -> random number generator
    !!
    subroutine sampleType(self, state, rand)
      import                                     :: configSource, RNG, transportObjectState
      class(configSource), intent(inout)         :: self
      class(transportObjectState), intent(inout) :: state
      type(RNG), intent(inout)                  :: rand
    end subroutine sampleType

  end interface

contains
  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(configSource), intent(inout) :: self

    call kill_super(self)

  end subroutine kill

  !!
  !! Sample particle's phase space co-ordinates
  !!
  !! See source_inter for details
  !!
  subroutine sampleState(self, rand, state)
    class(configSource), intent(inout)                    :: self
    type(RNG), intent(inout)                             :: rand
    class(transportObjectState), allocatable, intent(out) :: state

    call self % sampleType(state, rand)
    call self % samplePosition(state, rand)
    call self % sampleEnergyAngle(state, rand)
    call self % sampleEnergy(state, rand)
    call state % setTime(ZERO)
    call state % setWeight(ONE)

  end subroutine sampleState

end module configSource_inter
