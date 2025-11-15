module MGCollisionData_class

  use collisionData_class, only : collisionData
  use errors_mod,          only : fatalError
  use numPrecision

  implicit none
  private

  ! Public procedures.
  public :: castMGCollisionDataPtr

  !!
  !!
  !!
  type, public, extends(collisionData) :: MGCollisionData
    integer(shortInt) :: G_in = 0, G_out = 0
  end type MGCollisionData

contains
  !!
  !!
  !!
  function castMGCollisionDataPtr(source, fatal) result(ptr)
    class(collisionData), intent(in)       :: source
    logical(defBool), intent(in), optional :: fatal
    logical(defBool)                       :: throwError
    type(MGCollisionData), pointer         :: ptr
    character(*), parameter                :: HERE = 'castMGCollisionDataPtr (neutronMGCollisionProcessor_inter.f90)'

    select type(temp => source)
      type is(MGCollisionData)
        ptr => temp

      class default
        ptr => null()

    end select

    ! Throw error unless specified otherwise.
    throwError = .true.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) &
    call fatalError(HERE, "Collision data is not of type 'MGCollisionData'.")

  end function castMGCollisionDataPtr  

end module MGCollisionData_class