module triangulationFactory_func

    use centroidTriangulationMethod_class,  only : centroidTriangulationMethod
    use dictionary_class,                   only : dictionary
    use DompierreTriangulationMethod_class, only : DompierreTriangulationMethod
    use genericProcedures,                  only : fatalError
    use noTriangulation_class,              only : noTriangulation
    use numPrecision
    use triangulationMethod_inter,          only : triangulationMethod

    implicit none
    private

    ! List which contains acceptable types of triangulation methods.
    ! NOTE: It is necessary to adjust trailing blanks so all entries have the same length
    character(nameLen), dimension(*), parameter :: AVAILABLE_TRIANGULATIONS = ['centroidBased', &
                                                                               'Dompierre    ', &
                                                                               'none         ']

    ! Public interface.
    public :: newTriangulationPtr

contains
    !!
    !!
    !!
    subroutine newTriangulationPtr(dict, ptr)
        class(dictionary), intent(in)                    :: dict
        class(triangulationMethod), pointer, intent(out) :: ptr
        character(nameLen)                               :: type
        character(*), parameter                          :: here = 'newTriangulationPtr (triangulationFactory_func.f90)'

        ! Get type of triangulation method from dictionary. Default to 'none' if user has not specified one.
        call dict % getOrDefault(type, 'triangulationMethod', 'none')
        select case(type)
            case('none')
                allocate(noTriangulation :: ptr)

            case('centroidBased')
                allocate(centroidTriangulationMethod :: ptr)

            case('Dompierre')
                allocate(DompierreTriangulationMethod :: ptr)

            case default
                print '(A)', 'AVAILABLE TRIANGULATION METHODS: '
                print '(A)', AVAILABLE_TRIANGULATIONS
                call fatalError(here, 'Unrecognised type of triangulation method: '//trim(type)//'.')

        end select

    end subroutine newTriangulationPtr

end module triangulationFactory_func