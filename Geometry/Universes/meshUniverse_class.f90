module meshUniverse_class
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use numPrecision
  use universalVariables,           only : NUDGE, ZERO, ONE, INF, nameLen
  use genericProcedures,            only : fatalError, numToChar, countCharacters
  use dictionary_class,             only : dictionary
  use box_class,                    only : box
  use coord_class,                  only : coord
  use charMap_class,                only : charMap
  use sphere_class,                 only : sphere
  use surface_inter,                only : surface
  use surfaceShelf_class,           only : surfaceShelf
  use cell_inter,                   only : cell
  use cellShelf_class,              only : cellShelf
  use simpleCell_class,             only : simpleCell
  use mesh_inter,                   only : mesh
  use meshShelf_class,              only : meshShelf
  use universe_inter,               only : universe, kill_super => kill, charToFill
  
  implicit none
  private
  
  !!
  !! Local helper class to group cell data
  !!
  !! Public Members:
  !!   idx -> cellIdx of the cell in cellShelf
  !!   ptr -> Pointer to the cell
  !!
  type, private :: localCell
    integer(shortInt)    :: idx = 0
    class(cell), pointer :: ptr => null()
  end type localCell
  
  !!
  !! Local helper class to group mesh data
  !!
  !! Public Members:
  !!   idx -> meshIdx of the mesh in meshShelf
  !!   ptr -> Pointer to the mesh
  !!
  type, private :: localMesh
    integer(shortInt)    :: idx = 0
    class(mesh), pointer :: ptr => null()
  end type localMesh
  
  !!
  !! Special type of universe containing a mesh geometry. It is defined by a single simpleCell. 
  !! Note: for simplicity, a meshUniverse can only contain a single mesh. If multiple meshes are 
  !! to be used in the same simulation then multiple meshUniverses need to be defined.
  !!
  !! Sample Input Dictionary:
  !!   uni { type meshUniverse;
  !!         id 7;
  !!         # origin (2.0 0.0 0.0);    #
  !!         # rotation (23.0 0.0 0.0); #
  !!         cell 1;
  !!         mesh 2;
  !!         fills (mat1 mat2 mat3 ...);  }
  !!
  !! Notes:
  !!   - Fills are assigned in order as in definition. For example in an OpenFOAM mesh this would
  !!     lead to the following mapping [fill1: cell_zone1, fill2: cell_zone2, ...]. In tracking no
  !!     further checks are made, so be careful.
  !!   - A check is made on initialisation to ensure that the CSG cell does not crop the mesh
  !!     geometry.
  !!
  !! Public Members:
  !!   cell -> Structure that stores cellIdx and a pointer to the cell.
  !!   mesh -> Structure that stores meshIdx and a pointer to the mesh geometry.
  !!
  !! Interface:
  !!   universe interface
  !!
  type, public, extends(universe) :: meshUniverse
    type(localCell)               :: cell
    type(localMesh)               :: mesh
  contains
    ! Superclass procedures
    procedure                     :: init
    procedure                     :: kill
    procedure                     :: findCell
    procedure                     :: distance
    procedure                     :: cross
    procedure                     :: cellOffset
    ! Local procedures
    procedure                     :: checkForCropping
  end type meshUniverse
contains
  
  !!
  !! Initialise Universe
  !!
  !! See universe_inter for details.
  !!
  subroutine init(self, dict, mats, fills, cells, surfs, meshes)
    class(meshUniverse), intent(inout)                        :: self
    class(dictionary), intent(in)                             :: dict
    type(charMap), intent(in)                                 :: mats
    integer(shortInt), dimension(:), allocatable, intent(out) :: fills
    type(cellShelf), intent(inout)                            :: cells
    type(surfaceShelf), intent(inout)                         :: surfs
    type(meshShelf), intent(inout)                            :: meshes
    integer(shortInt)                                         :: cellId, i, meshId, nFills
    character(nameLen), dimension(:), allocatable             :: fillNames
    character(100), parameter                                 :: Here = 'init (meshUniverse_class.f90)'
    
    ! Setup the base class
    ! With: id, origin rotations...
    call self % setupBase(dict)
    
    ! Load meshId, convert meshId to meshIdx and get pointer to the mesh with corresponding idx.
    call dict % get(meshId, 'mesh')
    self % mesh % idx = meshes % getIdx(meshId)
    self % mesh % ptr => meshes % getPtr(self % mesh % idx)

    ! Load cellId, covert cellId to cellIdx and get pointer to the cell with corresponding idx.
    call dict % get(cellId, 'cell')
    self % cell % idx = cells % getIdx(cellId)
    self % cell % ptr => cells % getPtr(self % cell % idx)
    
    ! Check that the CSG cell does not crop the mesh.
    call self % checkForCropping(surfs)
    
    ! Retrieve fills in the dict and get local pointer to the mesh. Call fatalError if the number 
    ! of fills does not match the number of element zones in the mesh.
    call dict % get(fillNames, 'fills')
    nFills = size(fillNames)
    if (self % mesh % ptr % getLocalIdsNumber() /= nFills) call fatalError(Here, &
    'The number of fills does not match the number of element zones in mesh geometry with id: '//numToChar(meshId)//'.')
    
    ! Create fill array. First entry is fill of the CSG cell, remaining entries are the fills for 
    ! the pseudo-cells.
    allocate(fills(nFills + 1))
    fills(1) = cells % getFill(self % cell % idx)
    do i = 1, nFills
      fills(1 + i) = charToFill(fillNames(i), mats, Here)

    end do

  end subroutine init
  
  !!
  !! Find local cell ID given a point
  !!
  !! See universe_inter for details.
  !!
  subroutine findCell(self, coords)
    class(meshUniverse), intent(inout) :: self
    type(coord), intent(inout)         :: coords
    
    ! Set cellIdx to the index of the CSG cell, then find elementIdx and localId within mesh.
    call coords % setCellIdx(self % cell % idx)
    call self % mesh % ptr % findHostElement(coords)

  end subroutine findCell
  
  !!
  !! Returns distance to the next boundary between local cells in the universe
  !!
  !! See universe_inter for details.
  !!
  subroutine distance(self, coords, d, surfIdx)
    class(meshUniverse), intent(inout) :: self
    type(coord), intent(inout)         :: coords
    real(defReal), intent(out)         :: d
    integer(shortInt), intent(out)     :: surfIdx
    logical(defBool)                   :: inside
    
    ! Initialise surfIdx = 0 and compute distance to next mesh crossing. Also check if particle is
    ! inside the mesh.
    surfIdx = 0
    call self % mesh % ptr % distance(d, coords, inside)
    
    ! If particle is outside the mesh then compute distance to the next CSG surface crossing.
    if (.not. inside) call self % cell % ptr % distance(d, surfIdx, coords % getPosition(), coords % getDirection())

  end subroutine distance
  
  !!
  !! Cross between local cells
  !!
  !! See universe_inter for details.
  !!
  !! Note: Introduces extra movement to the particle to push it over boundary
  !!   for more efficent search. Distance is NUGDE.
  !!
  subroutine cross(self, coords, surfIdx)
    class(meshUniverse), intent(inout) :: self
    type(coord), intent(inout)         :: coords
    integer(shortInt), intent(in)      :: surfIdx

    ! Do nothing.

  end subroutine cross
  
  !!
  !! Returns offset for the current cell
  !!
  !! See universe_inter for details.
  !!
  function cellOffset(self, coords) result (offset)
    class(meshUniverse), intent(in) :: self
    type(coord), intent(in)         :: coords
    real(defReal), dimension(3)     :: offset
    
    ! There is no cell offset.
    offset = ZERO

  end function cellOffset
  
  !!
  !! Returns to uninitialised state
  !!
  elemental subroutine kill(self)
    class(meshUniverse), intent(inout) :: self
    
    ! Superclass.
    call kill_super(self)
    
    ! Local.
    self % cell % idx = 0
    self % mesh % idx = 0
    if (associated(self % cell % ptr)) nullify(self % cell % ptr)
    if (associated(self % mesh % ptr)) nullify(self % mesh % ptr)

  end subroutine kill
  
  !! Subroutine 'checkForCropping'
  !!
  !! Basic description:
  !!   Performs checks to ensure that the surface defining the CSG cell does not crop the mesh
  !!   contained in the universe.
  !!
  !! Arguments:
  !!   surfs [in] -> A surfaceShelf.
  !!
  !! Error:
  !!   Calls fatalError if the surface crops the mesh (that is, if the mesh is not fully contained
  !!   in the surface defining the CSG cell).
  !!
  subroutine checkForCropping(self, surfs)
    class(meshUniverse), intent(in)              :: self
    type(surfaceShelf), intent(in)               :: surfs
    class(cell), pointer                         :: cellPtr
    class(surface), pointer                      :: surfPtr
    integer(shortInt), dimension(:), allocatable :: surfIdxs
    type(axisAlignedBoundingBox), pointer        :: boundingBoxPtr
    real(defReal), dimension(3, 2)               :: bounds
    character(*), parameter                      :: Here = 'checkForCropping (meshUniverse_class.f90)'
    
    ! Get local pointer to cell. We need this to select the cell type.
    cellPtr => self % cell % ptr
    select type(cellPtr)
      ! If need to include more cell types in the future do it here.
      type is (simpleCell)
      
      ! Retrieve the index of the surface making the CSG cell.
      surfIdxs = cellPtr % getSurfaces()
      
    end select
    
    ! Call fatalError if the CSG cell uses more than one surface (this is to prevent the mesh from
    ! being put in the non-overlapping region between surfaces).
    if (size(surfIdxs) > 1) call fatalError(Here, 'The CSG cell used in the mesh universe has more than one surface.')
    
    ! Get pointer to the surface of the CSG cell and check that the surface of the CSG cell does not crop it.
    surfPtr => surfs % getPtr(abs(surfIdxs(1)))
    boundingBoxPtr => self % mesh % ptr % getBoundingBoxPtr()
    if (surfPtr % cropsBoundingBox(boundingBoxPtr % getBounds())) then
      bounds = boundingBoxPtr % getBounds()
      print *, 'Minimum x-coordinate: '//numToChar(bounds(1, 1))//'.'
      print *, 'Minimum y-coordinate: '//numToChar(bounds(2, 1))//'.'
      print *, 'Minimum z-coordinate: '//numToChar(bounds(3, 1))//'.'
      print *, 'Maximum x-coordinate: '//numToChar(bounds(1, 2))//'.'
      print *, 'Maximum y-coordinate: '//numToChar(bounds(2, 2))//'.'
      print *, 'Maximum z-coordinate: '//numToChar(bounds(3, 2))//'.'
      call fatalError(Here, 'Surface with id: '//numToChar(surfPtr % getId())//&
                      ' crops the bounding box of the mesh geometry with id: '//numToChar(self % mesh % ptr % getId())//'.')

    end if

    ! Nullify pointers.
    nullify(cellPtr)
    nullify(surfPtr)

  end subroutine checkForCropping
  
end module meshUniverse_class