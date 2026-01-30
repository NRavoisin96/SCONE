module accelerationStructureFactory_func

  use accelerationStructure_inter,    only : accelerationStructure
  use dictionary_class,               only : dictionary
  use edgeShelf_class,                only : edgeShelf
  use elementShelf_class,             only : elementShelf
  use errors_mod,                     only : fatalError
  use faceShelf_class,                only : faceShelf
  use numPrecision
  use octreeAcceleration_class,       only : octreeAcceleration
  use patchSearchAcceleration_class,  only : patchSearchAcceleration
  use standardASCGAcceleration_class, only : standardASCGAcceleration
  use vertexShelf_class,              only : vertexShelf

  implicit none
  private

  ! Parameters.
  character(nameLen), dimension(*), parameter :: AVAILABLE_ACCELERATION_STRUCTURES = ['octreeAcceleration      ', &
                                                                                      'patchSearchAcceleration ', &
                                                                                      'standardASCGAcceleration']

  ! Public procedures.
  public :: new_accelerationStructure

contains
  !!
  !!
  !!
  subroutine new_accelerationStructure(dict, vertices, edges, faces, elements, new)
    type(dictionary), intent(in)                           :: dict
    type(vertexShelf), intent(in)                          :: vertices
    type(edgeShelf), intent(inout)                         :: edges
    type(faceShelf), intent(inout)                         :: faces
    type(elementShelf), intent(in)                         :: elements
    class(accelerationStructure), allocatable, intent(out) :: new
    character(nameLen)                                     :: type
    real(defReal)                                          :: t1, t2
    character(*), parameter :: HERE = 'new_accelerationStructure (accelerationStructureFactory_func.f90)'

    call cpu_time(t1)

    ! Get type from dictionary, allocate then initialise.
    call dict % get(type, 'type')
    select case(type)
      case('octreeAcceleration')
        allocate(octreeAcceleration :: new)

      case('patchSearchAcceleration')
        allocate(patchSearchAcceleration :: new)

      case('standardASCGAcceleration')
        allocate(standardASCGAcceleration :: new)

      case default
        print '(A)', AVAILABLE_ACCELERATION_STRUCTURES
        call fatalError(HERE, 'Unrecognised acceleration structure type: '//trim(type)//'.')

    end select
    call new % init(dict, vertices, edges, faces, elements)

    ! end timer for initialisation and print
    call cpu_time(t2)          ! CPU-time
    print*, "-------------------------------------------------------------"
    print*, "/\/\ Initialisation procedure time /\/\"
    print*, "CPU time: ", t2 - t1, " seconds"
    print*, "-------------------------------------------------------------"

  end subroutine new_accelerationStructure

end module accelerationStructureFactory_func