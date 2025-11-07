module collisionProcessorFactory_func

  use collisionProcessor_inter, only : collisionProcessor
  use dictionary_class,         only : dictionary
  use errors_mod,               only : fatalError
  use neutronCEimp_class,       only : neutronCEimp
  use neutronCEstd_class,       only : neutronCEstd
  use neutronMGimp_class,       only : neutronMGimp
  use neutronMGstd_class,       only : neutronMGstd
  use numPrecision
  use universalVariables,       only : P_NEUTRON_CE, P_NEUTRON_MG, P_PHOTON_CE, P_PHOTON_MG

  implicit none
  private

  public :: new_collisionProcessor

  ! List that contains all accaptable types of collisionProcessors
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen), dimension(*), parameter :: AVAILABLE_collisionProcessors = ['neutronCEimp', &
                                                                                  'neutronCEstd', &
                                                                                  'neutronMGimp', &
                                                                                  'neutronMGstd']

contains

  !!
  !! Allocate new allocatable collisionProcessor to a specific type
  !! If new is allocated it deallocates it
  !!
  subroutine new_collisionProcessor(dict, new, particleType)
    class(dictionary), intent(in)                         :: dict
    class(collisionProcessor), allocatable, intent(inout) :: new
    integer(shortInt), intent(out)                        :: particleType
    character(nameLen)                                    :: type
    character(*), parameter                               :: Here = 'new_collisionProcessor (collisionProcessorFactory_func.f90)'

    ! Deallocate new if allocated
    if (allocated(new)) deallocate(new)

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of collisionProcessor
    select case(type)
      case('neutronCEimp')
        allocate(neutronCEimp :: new)
        particleType = P_NEUTRON_CE

      case('neutronCEstd')
        allocate(neutronCEstd :: new)
        particleType = P_NEUTRON_CE

      case('neutronMGimp')
        allocate(neutronMGimp :: new)
        particleType = P_NEUTRON_MG

      case('neutronMGstd')
        allocate(neutronMGstd :: new)
        particleType = P_NEUTRON_MG

      case default
        print *, AVAILABLE_collisionProcessors
        call fatalError(Here, 'Unrecognised type of collisionProcessor: '//trim(type)//'.')

    end select

    ! Initialise new processor
    call new % init(dict)

  end subroutine new_collisionProcessor

end module collisionProcessorFactory_func
