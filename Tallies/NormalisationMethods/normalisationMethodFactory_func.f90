module normalisationMethodFactory_func

  use constantNormalisationMethod_class,       only : constantNormalisationMethod
  use dictionary_class,                        only : dictionary
  use genericProcedures,                       only : fatalError
  use normalisationMethod_inter,               only : normalisationMethod
  use numPrecision
  use volumeWeightedNormalisationMethod_class, only : volumeWeightedNormalisationMethod

  implicit none
  private

  ! List that contains all accaptable types of normalisation methods
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen), dimension(*), parameter :: AVAILABLE_NORMALISATIONMETHODS = ['constantNormalisationMethod', &
                                                                                   'volumeNormalisationMethod  ']

  ! Public interface.
  public :: newNormalisationMethod

contains
  !!
  !!
  !!
  function newNormalisationMethod(dict) result(new)
    class(dictionary), intent(in)           :: dict
    class(normalisationMethod), allocatable :: new
    character(nameLen)                      :: type
    character(*), parameter                 :: here = 'newNormalisationMethod (normalisationMethodFactory_func.f90)'

    call dict % get(type, 'type')
    select case(type)
      case('constantNormalisationMethod')
        allocate(constantNormalisationMethod :: new)

      case('volumeNormalisationMethod')
        allocate(volumeWeightedNormalisationMethod :: new)

      case default
        print '(A)', 'AVAILABLE NORMALISATION METHODS: '
        print '(A)', AVAILABLE_NORMALISATIONMETHODS
        call fatalError(here, 'Unrecognised normalisation method type: '//type//'.')

    end select

    ! Initialise.
    call new % init(dict)

  end function newNormalisationMethod

end module normalisationMethodFactory_func