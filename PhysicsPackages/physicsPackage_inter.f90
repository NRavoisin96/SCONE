module physicsPackage_inter

  use dictionary_class, only : dictionary
  use geometry_inter,   only : geometry
  use numPrecision

  implicit none
  private

  !!
  !! Abstract interface of physics Package
  !! Physics package is controles a calculation flow
  !! Each type of calculation has diffrent physics package
  !!
  type, public,abstract :: physicsPackage
    private
  contains
    procedure(init), deferred :: init
    procedure(kill), deferred :: kill
    procedure(run),deferred   :: run
  end type physicsPackage

  !!
  !!
  !!
  type, public :: initPhysicsPackagePayload
    type(dictionary), pointer :: dict => null()
    integer(shortInt)         :: geometryIdx = 0
    class(geometry), pointer  :: geometry => null()
  end type initPhysicsPackagePayload

  abstract interface
    !!
    !! Initialise Physics Package from dictionary
    !!
    subroutine init(self, payload)
      import                                      :: initPhysicsPackagePayload, physicsPackage
      class(physicsPackage), intent(inout)        :: self
      type(initPhysicsPackagePayload), intent(in) :: payload
    end subroutine init

    !!
    !! Deallocate memory used by physicsPackage
    !!
    subroutine kill(self)
      import :: physicsPackage
      class(physicsPackage), intent(inout) :: self
    end subroutine kill


    !!
    !! Run calculation in the physics package
    !!
    subroutine run(self)
      import :: physicsPackage
      class(physicsPackage), intent(inout) :: self
    end subroutine run
  end interface

end module physicsPackage_inter