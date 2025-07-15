module noTriangulation_class

    use topologicalObjectShelf_class, only : topologicalObjectShelf
    use triangulationMethod_inter,    only : triangulationMethod

    implicit none
    private

    !!
    !!
    !!
    type, public, extends(triangulationMethod) :: noTriangulation
        private
    contains
        procedure :: triangulate
    end type noTriangulation

contains
    !!
    !!
    !!
    subroutine triangulate(self, edges, elements, faces, vertices)
        class(noTriangulation), intent(in)          :: self
        type(topologicalObjectShelf), intent(inout) :: edges, elements, faces, vertices

        ! Do nothing.

    end subroutine triangulate

end module noTriangulation_class