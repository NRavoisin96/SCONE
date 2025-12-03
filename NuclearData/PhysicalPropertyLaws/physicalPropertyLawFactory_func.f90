module physicalPropertyLawFactory_func

  use constantPropertyLaw_class,   only : constantPropertyLaw
  use dictionary_class,            only : dictionary
  use errors_mod,                  only : fatalError
  use numPrecision
  use physicalPropertyLaw_inter,   only : physicalPropertyLaw
  use polynomialPropertyLaw_class, only : polynomialPropertyLaw

  implicit none
  private

  ! Parameters.
  character(nameLen), dimension(*), parameter :: AVAILABLE_PHYSICAL_PROPERTY_LAWS = ['constantPropertyLaw  ', &
                                                                                     'polynomialPropertyLaw']

  ! Public procedures.
  public :: new_physicalPropertyLaw

contains
  !!
  !!
  !!
  subroutine new_physicalPropertyLaw(dict, new)
    type(dictionary), intent(in)                         :: dict
    class(physicalPropertyLaw), allocatable, intent(out) :: new
    character(nameLen)                                   :: type
    character(*), parameter                              :: HERE = 'new_physicalPropertyLaw (physicalPropertyLaw_func.f90)'

    ! Retrieve type from dictionary, allocate then initialise.
    call dict % get(type, 'type')
    
    select case(type)
      case('constantPropertyLaw')
        allocate(constantPropertyLaw :: new)

      case('polynomialPropertyLaw')
        allocate(polynomialPropertyLaw :: new)

      case default
        print '(A)', AVAILABLE_PHYSICAL_PROPERTY_LAWS
        call fatalError(HERE, 'Unrecognised physical property law type: '//trim(type)//'.')

    end select

    call new % init(dict)

  end subroutine new_physicalPropertyLaw

end module physicalPropertyLawFactory_func