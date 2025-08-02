!!
!! Module to generate physicsPackages
!!
module physicsPackageFactory_func

  use numPrecision
  use genericProcedures,               only : fatalError
  use dictionary_class,                only : dictionary

  ! Physics Package interface
  use physicsPackage_inter,            only : initPhysicsPackagePayload, physicsPackage

  ! Implementations
  use eigenPhysicsPackage_class,       only : eigenPhysicsPackage
  use fixedSourcePhysicsPackage_class, only : fixedSourcePhysicsPackage
  use vizPhysicsPackage_class,         only : vizPhysicsPackage
  use rayVolPhysicsPackage_class,      only : rayVolPhysicsPackage
!  use dynamPhysicsPackage_class, only : dynamPhysicsPackage

  implicit none
  private

  ! List that contains all accaptable types of Physics Packages
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen), dimension(*), parameter :: AVAILABLE_physicsPackages = ['eigenPhysicsPackage      ',&
                                                                              'fixedSourcePhysicsPackage',&
                                                                              'vizPhysicsPackage        ',&
                                                                              'rayVolPhysicsPackage     ']

  !!
  !! Public interface
  !!
  public :: new_physicsPackage
  public :: new_physicsPackage_ptr

contains

  !!
  !! Return a new allocatable instance of a physicsPackage
  !!
  function new_physicsPackage(payload) result(new)
    type(initPhysicsPackagePayload), intent(in) :: payload
    class(physicsPackage), allocatable          :: new
    character(nameLen)                          :: type
    character(*), parameter                     :: Here = 'new_physicsPackage (physicsPackageFactory_func.f90)'

    ! Obtain string that specifies type to be built
    call payload % dict % get(type, 'type')

    ! Allocate approperiate subclass of physicsPackage
    select case(type)
      case('eigenPhysicsPackage')
        allocate(eigenPhysicsPackage :: new)

      case('fixedSourcePhysicsPackage')
        allocate(fixedSourcePhysicsPackage :: new)

      case('vizPhysicsPackage')
        allocate(vizPhysicsPackage :: new)

      case('rayVolPhysicsPackage')
        allocate(rayVolPhysicsPackage :: new)

      case default
        print *, AVAILABLE_physicsPackages
        call fatalError(Here, 'Unrecognised type of Physics Package: ' //trim(type)//'.')

    end select

    ! Initialise new physics package
    call new % init(payload)

  end function new_physicsPackage

  !!
  !! Returns allocated pointer to Physics Package build with the dictionary dict
  !!
  function new_physicsPackage_ptr(payload) result(new)
    type(initPhysicsPackagePayload), intent(in) :: payload
    class(physicsPackage), pointer              :: new
    character(*), parameter                     :: Here = 'new_physicsPackage_ptr (physicsPackageFactory_func.f90)'

    ! Allocate pointer and copy data from local allocatable
    allocate(new, source = new_physicsPackage(payload))

  end function new_physicsPackage_ptr

end module physicsPackageFactory_func
