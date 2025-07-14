module accelerationStructureFactory_func

  use accelerationStructure_inter,   only : accelerationStructure
  use dictionary_class,              only : dictionary
  use edgeShelf_class,               only : edgeShelf
  use elementShelf_class,            only : elementShelf
  use faceShelf_class,               only : faceShelf
  use genericProcedures,             only : fatalError
  use noAcceleration_class,          only : noAcceleration
  use numPrecision
  use octreeAcceleration_class,      only : octreeAcceleration
  use patchSearchAcceleration_class, only : patchSearchAcceleration
  use vertexShelf_class,             only : vertexShelf

  implicit none
  private

  ! List which contains acceptable types of acceleration structures.
  ! NOTE: It is necessary to adjust trailing blanks so all entries have the same length.
  character(nameLen), dimension(*), parameter :: AVAILABLE_ACCELERATIONS = ['none       ', &
                                                                            'octree     ', &
                                                                            'patchSearch']

  ! Public interface.
  public :: newAccelerationStructurePtr

contains
  !!
  !!
  !!
  subroutine newAccelerationStructurePtr(dict, edges, elements, faces, vertices, ptr)
    class(dictionary), intent(in)                      :: dict
    type(edgeShelf), intent(in)                        :: edges
    type(elementShelf), intent(in)                     :: elements
    type(faceShelf), intent(in)                        :: faces
    type(vertexShelf), intent(in)                      :: vertices
    class(accelerationStructure), pointer, intent(out) :: ptr
    character(nameLen)                                 :: type
    character(*), parameter :: here = 'newAccelerationStructurePtr (accelerationStructureFactory_func.f90)'

    ! Get type of acceleration structure from dictionary. Default to 'none' if user has not specified one.
    call dict % getOrDefault(type, 'accelerationMethod', 'none')
    select case(type)
        case('none')
            allocate(noAcceleration :: ptr)

        case('octree')
            allocate(octreeAcceleration :: ptr)

        case('patchSearch')
            allocate(patchSearchAcceleration :: ptr)

        case default
            print '(A)', 'AVAILABLE ACCELERATION STRUCTURES: '
            print '(A)', AVAILABLE_ACCELERATIONS
            call fatalError(here, 'Unrecognised type of acceleration structure: '//trim(type)//'.')

    end select

    ! Initialise acceleration structure.
    call ptr % init(dict, edges, elements, faces, vertices)

  end subroutine newAccelerationStructurePtr

end module accelerationStructureFactory_func