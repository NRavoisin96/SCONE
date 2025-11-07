module MGParticle_inter

  use errors_mod,             only : fatalError
  use MGParticleState_class,  only : castMGParticleStatePtr, MGParticleState
  use numPrecision
  use physicalParticle_inter, only : physicalParticle

  implicit none
  private

  ! Public procedures.
  public :: castMGParticlePtr

  !!
  !!
  !!
  type, public, abstract, extends(physicalParticle) :: MGParticle
    private
  contains
    procedure :: getEnergyGroup
    procedure :: getSpeed
    procedure :: setEnergyGroup
  end type MGParticle

contains
  !!
  !!
  !!
  function castMGParticlePtr(source, fatal) result(ptr)
    class(physicalParticle), intent(in)    :: source
    logical(defBool), intent(in), optional :: fatal
    class(MGParticle), pointer             :: ptr
    logical(defBool)                       :: throwError
    character(*), parameter                :: here = 'castMGParticlePtr (MGParticle_inter.f90)'

    select type(temp => source)
      class is(MGParticle)
        ptr => temp

      class default
        ptr => null()

    end select

    ! Throw error if requested.
    throwError = .false.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) call fatalError(here, "Physical particle is not of class 'MGParticle'.")

  end function castMGParticlePtr

  !!
  !!
  !!
  function getEnergyGroup(self) result(energyGroup)
    class(MGParticle), intent(in)   :: self
    integer(shortInt)               :: energyGroup
    class(MGParticleState), pointer :: MGParticleStatePtr

    MGParticleStatePtr => castMGParticleStatePtr(self % getCurrentStatePtr(), .true.)
    energyGroup = MGParticleStatePtr % getEnergyGroup()

  end function getEnergyGroup

  !!
  !!
  !!
  function getSpeed(self) result(speed)
    class(MGParticle), intent(in) :: self
    real(defReal)                 :: speed
    character(*), parameter       :: HERE = 'getSpeed (MGParticle_inter.f90)'

    call fatalError(HERE, 'Unsupported procedure.')
    speed = ZERO

  end function getSpeed

  !!
  !!
  !!
  subroutine setEnergyGroup(self, energyGroup)
    class(MGParticle), intent(inout) :: self
    integer(shortInt), intent(in)    :: energyGroup
    class(MGParticleState), pointer  :: MGParticleStatePtr

    MGParticleStatePtr => castMGParticleStatePtr(self % getCurrentStatePtr(), .true.)
    call MGParticleStatePtr % setEnergyGroup(energyGroup)

  end subroutine setEnergyGroup

end module MGParticle_inter