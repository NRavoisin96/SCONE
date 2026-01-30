module octree_class

  use coord_class,        only : coord
  use dictionary_class,   only : dictionary
  use elementShelf_class, only : elementShelf
  use errors_mod,         only : fatalError
  use faceShelf_class,    only : faceShelf
  use genericProcedures,  only : areEqual, numToChar
  use objectKDTree_class, only : objectKDTree
  use mesh_inter,         only : mesh
  use numPrecision
  use octreeNode_class,   only : octreeNode
  use universalVariables, only : NUDGE
  use vertexShelf_class,  only : vertexShelf

  implicit none
  private

  type, public :: octree
    private
    integer(shortInt) :: nLeaves = 0, maxFacesNumber = 0, maxRefinementLevel = 0
    type(octreeNode)  :: root
  contains
    ! Build procedures.
    procedure :: init
    procedure :: kill
    ! Runtime procedures.
    procedure :: countInside
    procedure :: countOutside
    procedure :: findLeaf
  end type octree

contains
  !!
  !!
  !!
  elemental function countInside(self) result(nInside)
    class(octree), intent(in) :: self
    integer(shortInt)         :: nInside

    nInside = 0
    call self % root % countInside(nInside)

  end function countInside

  !!
  !!
  !!
  elemental function countOutside(self) result(nOutside)
    class(octree), intent(in) :: self
    integer(shortInt)         :: nOutside

    nOutside = 0
    call self % root % countOutside(nOutside)

  end function countOutside

  !!
  !!
  !!
  subroutine findLeaf(self, coords, leaf)
    class(octree), intent(in)              :: self
    type(coord), intent(inout)             :: coords
    type(octreeNode), pointer, intent(out) :: leaf

    ! First find the leaf node containing the coordinates.
    call self % root % findLeaf(coords, leaf, .true.)

  end subroutine findLeaf

  !!
  !!
  !!
  subroutine init(self, dict, vertices, faces, elements)
    class(octree), intent(inout)   :: self
    type(dictionary), intent(in)   :: dict
    type(vertexShelf), intent(in)  :: vertices
    type(faceShelf), intent(in)    :: faces
    type(elementShelf), intent(in) :: elements
    type(objectKDTree)             :: tree
    real(defReal), dimension(6)    :: boundingBoxBounds
    character(*), parameter        :: HERE = 'init (octree_class.f90)'

    call dict % getOrDefault(self % maxRefinementLevel, 'depth', 10)
    if(self % maxRefinementLevel < 1) &
    call fatalError(HERE, 'Depth must be at least 1. Is: '//numToChar(self % maxRefinementLevel)//'.')
    
    call dict % getOrDefault(self % maxFacesNumber, 'nMaxFaces', 4)
    if(self % maxFacesNumber < 1) &
    call fatalError(HERE, 'nMaxFaces must be at least 1. Is: '//numToChar(self % maxFacesNumber)//'.')

    ! Initialise k-d trees from the unstructured mesh faces and elements, then build the Cartesian grid's 
    ! root cell and all its children cells.
    call tree % init(faces % getAllFaceCentroids(), faces % getAllFaceBoundingBoxes())

    boundingBoxBounds = tree % getRootBoundingBoxBounds() + [-NUDGE, -NUDGE, -NUDGE, NUDGE, NUDGE, NUDGE]
    call self % root % init(boundingBoxBounds, tree, vertices, faces, 1, self % maxFacesNumber, self % maxRefinementLevel, &
                            self % nLeaves)

    ! After the refinement process, assign empty cells to their correct element.
    call self % root % assignElement(tree, vertices, faces, elements)

    ! Kill the k-d tree as it is no longer needed.
    call tree % kill()

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(octree), intent(inout) :: self

    ! Local.
    self % nLeaves = 0
    self % maxRefinementLevel = 0
    self % maxFacesNumber = 0
    call self % root % kill()

  end subroutine kill

end module octree_class