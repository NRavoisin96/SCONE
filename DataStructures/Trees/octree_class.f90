module octree_class

  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use coord_class,                  only : coord
  use dictionary_class,             only : dictionary
  use elementShelf_class,           only : elementShelf
  use faceShelf_class,              only : faceShelf
  use genericProcedures,            only : areEqual, fatalError
  use objectKDTree_class,           only : objectKDTree
  use mesh_inter,                   only : mesh
  use numPrecision
  use octreeNode_class,             only : octreeNode
  use universalVariables,           only : NUDGE
  use vertexShelf_class,            only : vertexShelf

  implicit none
  private

  type, public :: octree
    private
    integer(shortInt)  :: nLeaves = 0, maxFacesNumber = 4, maxRefinementLevel = 10
    type(octreeNode)   :: root
  contains
    ! Build procedures.
    procedure          :: init
    procedure          :: kill
    ! Runtime procedures.
    procedure          :: countInside
    procedure          :: countOutside
    procedure          :: findLeaf
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
  subroutine init(self, vertices, faces, elements)
    class(octree), intent(inout)   :: self
    type(vertexShelf), intent(in)  :: vertices
    type(faceShelf), intent(in)    :: faces
    type(elementShelf), intent(in) :: elements
    type(objectKDTree)             :: tree
    type(axisAlignedBoundingBox)   :: boundingBox, rootBoundingBox

    ! Initialise k-d trees from the unstructured mesh faces and elements, then build the Cartesian grid's 
    ! root cell and all its children cells.
    call tree % init(faces % getAllFaceCentroids(), faces % getAllFaceBoundingBoxes())

    rootBoundingBox = tree % getRootBoundingBox()
    call boundingBox % init(rootBoundingBox % getBounds() + [-NUDGE, -NUDGE, -NUDGE, NUDGE, NUDGE, NUDGE])
    call self % root % init(tree, vertices, faces, boundingBox, 1, &
                            self % maxFacesNumber, self % maxRefinementLevel, self % nLeaves)

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
    call self % root % kill()

  end subroutine kill

end module octree_class