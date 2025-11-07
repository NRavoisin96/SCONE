!!
!! Transport operator for delta tracking
!!
module transportOperatorDT_class
  
  use dictionary_class,        only : dictionary
  use errors_mod,              only : fatalError
  use numPrecision
  use physicalParticle_inter,  only : castPhysicalParticlePtr, physicalParticle
  use RNG_class,               only : RNG
  use tallyAdmin_class,        only : tallyAdmin
  use tallyCodes
  use transportObject_inter,   only : transportObject
  use transportOperator_inter, only : transportOperator, init_super => init
  use universalVariables

  implicit none
  private

  !!
  !! Transport operator that moves a particle with delta tracking
  !!
  type, public, extends(transportOperator) :: transportOperatorDT
    private
  contains
    procedure :: transit => deltaTracking
  end type transportOperatorDT

contains

  !!
  !! Performs delta tracking until a real collision point is found
  !!
  subroutine deltaTracking(self, object, tally)
    class(transportOperatorDT), intent(inout) :: self
    class(transportObject), intent(inout)     :: object
    type(tallyAdmin), intent(inout)           :: tally
    class(physicalParticle), pointer          :: p
    integer(shortInt)                         :: materialIdx
    real(defReal)                             :: majorant_inv, sigmaT, distance, randomNumber
    type(RNG), pointer                        :: RNGPtr
    character(*), parameter                   :: Here = 'deltaTracking (transportOperatorDT_class.f90)'

    ! Downcast transport object into a physical particle.
    p => castPhysicalParticlePtr(object, .true.)

    ! Get pointer to RNG.
    RNGPtr => p % getRNGPtr()

    ! Get majorant XS inverse: 1/Sigma_majorant
    majorant_inv = ONE / self % getTrackingXS(p, p % getMaterialIdx(), MAJORANT_XS)

   ! Should never happen! Prevents Inf distances
    if (abs(majorant_inv) > huge(majorant_inv)) call fatalError(Here, "Majorant cross section is ZERO.")

    DTLoop: do
      call RNGPtr % generateDistance(majorant_inv, distance)

      ! Move partice in the geometry and get materialIdx following teleportation.
      call self % teleport(p, distance)
      materialIdx = p % getMaterialIdx()

      select case(materialIdx)
        case(OUTSIDE_FILL)
          ! If particle has leaked, exit
          call p % setFate(LEAK_FATE)
          call p % setIsDead(.true.)
          return

        case(VOID_MAT)
          ! Check for void.
          call tally % reportInColl(p, .true.)

        case(UNDEF_MAT)
          ! Give error if the particle somehow ended in an undefined material.
          print *, p % getGlobalPosition()
          call fatalError(Here, 'Particle is in undefined material.')

        case default
          ! Obtain the local cross-section
          sigmaT = self % getTrackMatXS(p, materialIdx)

          ! Roll RNG to determine if the collision is real or virtual
          ! Exit the loop if the collision is real, report collision if virtual
          call RNGPtr % generate(randomNumber)
          if (randomNumber < sigmaT * majorant_inv) then
            exit DTLoop

          else
            call tally % reportInColl(p, .true.)

          end if
          
      end select

    end do DTLoop

    call tally % reportTrans(object)

  end subroutine deltaTracking

end module transportOperatorDT_class