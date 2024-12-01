module edge_class
  
  use numPrecision
  
  implicit none
  private
  
  !!
  !! Edge of a mesh linking two vertices.
  !!
  !! Private members:
  !!   idx            -> Index of the edge.
  !!   startVertexIdx -> Index of the first vertex in the edge.
  !!   endVertexIdx   -> Index of the end vertex in the edge.
  !!   edgeToFaces    -> Array that stores edge-to-faces connectivity information.
  !!   edgeToElements -> Array that stores edge-to-elements connectivity information.
  !!
  type, public :: edge
    private
    integer(shortInt)                            :: idx, startVertexIdx, endVertexIdx
    integer(shortInt), dimension(:), allocatable :: edgeToFaces, edgeToElements
  end type edge
end module edge_class