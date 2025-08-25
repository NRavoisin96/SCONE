module tallyClerkFactory_func

  use centreOfMassClerk_class,         only : centreOfMassClerk
  use collisionClerk_class,            only : collisionClerk
  use collisionProbabilityClerk_class, only : collisionProbabilityClerk
  use dancoffBellClerk_class,          only : dancoffBellClerk
  use dictionary_class,                only : dictionary
  use errors_mod,                      only : fatalError
  use keffAnalogClerk_class,           only : keffAnalogClerk
  use keffImplicitClerk_class,         only : keffImplicitClerk
  use mgXsClerk_class,                 only : mgXsClerk
  use numPrecision
  use shannonEntropyClerk_class,       only : shannonEntropyClerk
  use simpleFMClerk_class,             only : simpleFMClerk
  use tallyClerk_inter,                only : tallyClerk
  use temperatureClerk_class,          only : temperatureClerk
  use trackClerk_class,                only : trackClerk

  implicit none
  private

  public :: new_tallyClerk

  ! List that contains all accaptable types of tallyClerks
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen), dimension(*), parameter :: AVAILABLE_tallyClerks = ['centreOfMassClerk        ', &
                                                                          'collisionClerk           ', &
                                                                          'collisionProbabilityClerk', &
                                                                          'dancoffBellClerk         ', &
                                                                          'keffAnalogClerk          ', &
                                                                          'keffImplicitClerk        ', &
                                                                          'mgXsClerk                ', &
                                                                          'shannonEntropyClerk      ', &
                                                                          'simpleFMClerk            ', &
                                                                          'temperatureClerk         ', &
                                                                          'trackClerk               ']

contains

  !!
  !! Allocate new allocatable tallyClerk to a specific type
  !! If new is allocated it deallocates it
  !!
  subroutine new_tallyClerk(new, dict, name)
    class(tallyClerk), allocatable, intent(inout) :: new
    class(dictionary), intent(in)                :: dict
    character(nameLen), intent(in)                :: name
    character(nameLen)            :: type
    character(100), parameter      :: Here = 'new_tallyClerk (tallyClerkFactory_func.f90)'

    ! Deallocate new if allocated
    if (allocated(new)) deallocate(new)

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of tallyClerk
    select case(type)
      case('centreOfMassClerk')
        allocate(centreOfMassClerk :: new)

      case('collisionClerk')
        allocate(collisionClerk :: new)

      case('collisionProbabilityClerk')
        allocate(collisionProbabilityClerk :: new)

      case('dancoffBellClerk')
        allocate(dancoffBellClerk :: new)

      case('keffAnalogClerk')
        allocate(keffAnalogClerk :: new)

      case('keffImplicitClerk')
        allocate(keffImplicitClerk :: new)

      case('mgXsClerk')
        allocate(mgXsClerk :: new)

      case('shannonEntropyClerk')
        allocate(shannonEntropyClerk :: new)

      case('simpleFMClerk')
        allocate(simpleFMClerk :: new)

      case('temperatureClerk')
        allocate(temperatureClerk :: new)

      case('trackClerk')
        allocate(trackClerk :: new)

      case default
        print *, AVAILABLE_tallyClerks
        call fatalError(Here, 'Unrecognised type of tallyClerk: '//trim(type)//'.')

    end select

    ! Initialise new clerk
    call new % init(dict, name)

  end subroutine new_tallyClerk

end module tallyClerkFactory_func
