module sourceFactory_func

  use dictionary_class,      only : dictionary
  use errors_mod,            only : fatalError
  use CEFissionSource_class, only : CEFissionSource
  use geometry_inter,        only : geometry
  use numPrecision
  use source_inter,          only : source

  implicit none
  private

  public :: new_source

  ! List that contains all accaptable types of sources
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all entries have the same length
  character(nameLen), dimension(*), parameter :: AVAILABLE_sources = ['CEFissionSource']

contains

  !!
  !! Allocate new allocatable source to a specific type
  !! If new is allocated it deallocates it
  !!
  subroutine new_source(new, dict, geom)
    class(source), allocatable, intent(inout) :: new
    class(dictionary), intent(in)             :: dict
    class(geometry), pointer, intent(in)      :: geom
    character(nameLen)                        :: type
    character(*), parameter                   :: Here = 'new_source (sourceFactory_func.f90)'

    ! Deallocate new if allocated
    if (allocated(new)) deallocate(new)

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of source
    select case(type)
      case('CEFissionSource')
        allocate(CEFissionSource :: new)

      case default
        print *, AVAILABLE_sources
        call fatalError(Here, 'Unrecognised type of source: '//trim(type)//'.')

    end select

    ! Initialise new source
    call new % init(dict, geom)

  end subroutine new_source

end module sourceFactory_func
