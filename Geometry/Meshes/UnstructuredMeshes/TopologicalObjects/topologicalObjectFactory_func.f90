module topologicalObjectFactory_func

  use edge_class,                    only : edge
  use element_class,                 only : buildElementPayload, element
  use extentTopologicalObject_inter, only : buildExtentTopologicalObjectPayload
  use face_class,                    only : buildFacePayload, face
  use genericProcedures,             only : fatalError
  use numPrecision
  use topologicalObject_inter,       only : buildTopologicalObjectPayload, topologicalObjectBox
  use vertex_class,                  only : buildVertexPayload, vertex

  implicit none
  private

  ! List which contains acceptable types of triangulation methods.
  ! NOTE: It is necessary to adjust trailing blanks so all entries have the same length
  character(nameLen), dimension(*), parameter :: AVAILABLE_TOPOLOGICAL_OBJECTS_PAYLOADS = ['buildEdgePayload   ', &
                                                                                           'buildElementPayload', &
                                                                                           'buildFacePayload   ', &
                                                                                           'buildVertexPayload ']

  ! Public interface.
  public :: newTopologicalObjectBox

contains
  !!
  !!
  !!
  subroutine newTopologicalObjectBox(payload, box)
    class(buildTopologicalObjectPayload), intent(inout) :: payload
    type(topologicalObjectBox), intent(out)             :: box
    character(*), parameter                             :: here = 'newTopologicalObjectBox (topologicalObjectFactory_func.f90)'

    ! Allocate pointer inside box to correct class.
    select type(payload)
      type is(buildExtentTopologicalObjectPayload)
        allocate(edge :: box % ptr)

      type is(buildElementPayload)
        allocate(element :: box % ptr)

      type is(buildFacePayload)
        allocate(face :: box % ptr)

      type is(buildVertexPayload)
        allocate(vertex :: box % ptr)

      class default
        print '(A)', 'AVAILABLE TOPOLOGICAL OBJECT PAYLOADS: '
        print '(A)', AVAILABLE_TOPOLOGICAL_OBJECTS_PAYLOADS
        call fatalError(Here, 'Unrecognised topological object payload.')

    end select

    ! Initialise pointer inside box.
    call box % ptr % init(payload)

  end subroutine newTopologicalObjectBox

end module topologicalObjectFactory_func