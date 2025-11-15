module CECollisionData_class

  use collisionData_class, only : collisionData
  use errors_mod,          only : fatalError
  use numPrecision

  implicit none
  private

  ! Public procedures.
  public :: castCECollisionDataPtr

  !!
  !!
  !!
  type, public, extends(collisionData) :: CECollisionData
    real(defReal) :: finalEnergy = ZERO, initialEnergy = ZERO, reactionEnergy = ZERO
  end type CECollisionData

contains
  !!
  !!
  !!
  function castCECollisionDataPtr(source, fatal) result(ptr)
    class(collisionData), intent(in)       :: source
    logical(defBool), intent(in), optional :: fatal
    logical(defBool)                       :: throwError
    type(CECollisionData), pointer         :: ptr
    character(*), parameter                :: HERE = 'castCECollisionDataPtr (neutronCECollisionProcessor_inter.f90)'

    select type(temp => source)
      type is(CECollisionData)
        ptr => temp

      class default
        ptr => null()

    end select

    ! Throw error unless specified otherwise.
    throwError = .true.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) &
    call fatalError(HERE, "Collision data is not of type 'CECollisionData'.")

  end function castCECollisionDataPtr

end module CECollisionData_class