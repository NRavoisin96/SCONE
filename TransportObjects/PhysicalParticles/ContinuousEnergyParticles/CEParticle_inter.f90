module CEParticle_inter

  use CECollisionData_class,      only : castCECollisionDataPtr, CECollisionData
  use CEParticleState_class,      only : buildCEParticleStatePayload, castCEParticleStatePtr, CEParticleState
  use collisionData_class,        only : collisionData
  use errors_mod,                 only : fatalError
  use numPrecision
  use physicalParticle_inter,     only : physicalParticle, prepareCollisionData_super => prepareCollisionData
  use transportObjectState_class, only : buildTransportObjectStatePayload
  use universalVariables,         only : lightSpeed

  implicit none
  private

  ! Public procedures.
  public :: castCEParticlePtr

  !!
  !!
  !!
  type, public, abstract, extends(physicalParticle) :: CEParticle
    private
  contains
    procedure :: allocatePayload
    procedure :: getEnergy
    procedure :: getSpeed
    procedure :: prepareCollisionData
    procedure :: setEnergy
  end type CEParticle

contains
  !!
  !!
  !!
  subroutine allocatePayload(self, payload)
    class(CEParticle), intent(in)                                     :: self
    class(buildTransportObjectStatePayload), allocatable, intent(out) :: payload

    allocate(buildCEParticleStatePayload :: payload)

  end subroutine allocatePayload

  !!
  !!
  !!
  function castCEParticlePtr(source, fatal) result(ptr)
    class(physicalParticle), intent(in)    :: source
    logical(defBool), intent(in), optional :: fatal
    class(CEParticle), pointer             :: ptr
    logical(defBool)                       :: throwError
    character(*), parameter                :: here = 'castCEParticlePtr (CEParticle_inter.f90)'

    select type(temp => source)
      class is(CEParticle)
        ptr => temp

      class default
        ptr => null()

    end select

    ! Throw error if requested.
    throwError = .false.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) call fatalError(here, "Physical particle is not of class 'CEParticle'.")

  end function castCEParticlePtr

  !!
  !!
  !!
  function getEnergy(self) result(energy)
    class(CEParticle), intent(in)   :: self
    class(CEParticleState), pointer :: CEParticleStatePtr
    real(defReal)                   :: energy

    CEParticleStatePtr => castCEParticleStatePtr(self % getCurrentStatePtr(), .true.)
    energy = CEParticleStatePtr % getEnergy()

  end function getEnergy

  !!
  !!
  !!
  function getSpeed(self) result(speed)
    class(CEParticle), intent(in) :: self
    real(defReal)                 :: speed
    type(CEParticleState), pointer:: CEParticleStatePtr

    CEParticleStatePtr => castCEParticleStatePtr(self % getCurrentStatePtr(), .true.)
    speed = sqrt(TWO * CEParticleStatePtr % getEnergy() / CEParticleStatePtr % getMass()) * lightSpeed

  end function getSpeed

  !!
  !!
  !!
  subroutine prepareCollisionData(self, collDat)
    class(CEParticle), intent(in)       :: self
    class(collisionData), intent(inout) :: collDat
    type(CECollisionData), pointer      :: CECollisionDataPtr

    call prepareCollisionData_super(self, collDat)
    CECollisionDataPtr => castCECollisionDataPtr(collDat)
    CECollisionDataPtr % initialEnergy = self % getEnergy()

  end subroutine prepareCollisionData

  !!
  !!
  !!
  subroutine setEnergy(self, energy)
    class(CEParticle), intent(inout) :: self
    real(defReal), intent(in)        :: energy
    class(CEParticleState), pointer  :: CEParticleStatePtr

    CEParticleStatePtr => castCEParticleStatePtr(self % getCurrentStatePtr(), .true.)
    call CEParticleStatePtr % setEnergy(energy)

  end subroutine setEnergy

end module CEParticle_inter