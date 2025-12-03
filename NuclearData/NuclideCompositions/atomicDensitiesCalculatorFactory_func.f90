module atomicDensitiesCalculatorFactory_func

  use atomicDensitiesCalculator_inter, only : atomicDensitiesCalculator
  use atomicDensityFraction_class,     only : atomicDensityFraction
  use dictionary_class,                only : dictionary
  use errors_mod,                      only : fatalError
  use massFraction_class,              only : massFraction
  use numPrecision
  use rawAtomicDensities_class,        only : rawAtomicDensities

  implicit none
  private

  ! Parameters.
  character(nameLen), dimension(*), parameter :: AVAILABLE_ATOMIC_DENSITIES_CALCULATORS = ['atomicDensityFraction', &
                                                                                           'massFraction         ', &
                                                                                           'rawAtomicDensities   ']

  ! Public procedures.
  public :: new_atomicDensitiesCalculator

contains
  !!
  !!
  !!
  subroutine new_atomicDensitiesCalculator(dict, new)
    type(dictionary), intent(in)                               :: dict
    class(atomicDensitiesCalculator), allocatable, intent(out) :: new
    character(nameLen)                                         :: type
    character(*), parameter :: HERE = 'new_atomicDensitiesCalculator (atomicDensitiesCalculatorFactory_func.f90)'

    ! Retrieve type from dictionary, allocate then initialise.
    call dict % get(type, 'type')
    select case(type)
      case('atomicDensityFraction')
        allocate(atomicDensityFraction :: new)
        
      case('massFraction')
        allocate(massFraction :: new)
        
      case('rawAtomicDensities')
        allocate(rawAtomicDensities :: new)

      case default
        print '(A)', AVAILABLE_ATOMIC_DENSITIES_CALCULATORS
        call fatalError(HERE, 'Unrecognised atomic densities calculator type: '//trim(type)//'.')

    end select

  end subroutine new_atomicDensitiesCalculator

end module atomicDensitiesCalculatorFactory_func