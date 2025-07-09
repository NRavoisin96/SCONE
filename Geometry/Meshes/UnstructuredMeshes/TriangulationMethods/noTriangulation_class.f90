module noTriangulation_class

    use edgeShelf_class,           only : edgeShelf
    use elementShelf_class,        only : elementShelf
    use faceShelf_class,           only : faceShelf
    use triangulationMethod_inter, only : triangulationMethod
    use vertexShelf_class,         only : vertexShelf

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
        class(noTriangulation), intent(in) :: self
        type(edgeShelf), intent(inout)     :: edges
        type(elementShelf), intent(inout)  :: elements
        type(faceShelf), intent(inout)     :: faces
        type(vertexShelf), intent(inout)   :: vertices

        ! Do nothing.

    end subroutine triangulate

end module noTriangulation_class