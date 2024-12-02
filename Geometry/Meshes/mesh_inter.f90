module mesh_inter
  
  use numPrecision
  use universalVariables,  only : NUDGE
  use genericProcedures,   only : fatalError, numToChar, openToRead
  use cellZoneShelf_class, only : cellZoneShelf
  use dictionary_class,    only : dictionary
  use coord_class,         only : coord
  
  implicit none
  private
  
  ! Extendable methods.
  public :: kill
  
  !!
  !! Abstract interface for all meshes.
  !!
  !! A mesh represents a subdivision of the entire space into small 2- or 3-D elements.
  !!
  !! Private Members:
  !!   id                    -> Id of the mesh.
  !!   nElementZones         -> Number of element zones in the mesh.
  !!   boundingBox           -> Axis-aligned bounding box (AABB) of the mesh.
  !!   cellZonesFile         -> .true. if a file is present to explicitly assign element zones.
  !!   cellZones             -> Shelf that stores element zones.
  !!
  !! Interface:
  !!   getBoundingBox        -> Returns the bounding box of the mesh.
  !!   getElementZonesNumber -> Returns the number of element zones in the mesh.
  !!   getId                 -> Returns the id of the mesh.
  !!   setBoundingBox        -> Sets bounding box of the mesh.
  !!   setId                 -> Sets Id of the mesh.
  !!   init                  -> Initialises mesh from input files.
  !!   kill                  -> Returns to uninitialised state.
  !!   findElement           -> Returns local cellId and cellIdx in cellShelf for a given position.
  !!   distance              -> Calculates the distance travelled by a particle within the mesh.
  !!   distanceToNextFace    -> Calculates the distance to the next face in the mesh.
  !!
  type, public, abstract :: mesh
    private
    integer(shortInt), public               :: id = 0, nElementZones = 0
    real(defReal), dimension(6)             :: boundingBox = ZERO
    logical(defBool), public                :: cellZonesFile = .false.
    type(cellZoneShelf), public             :: cellZones
  contains
    ! Build procedures.
    procedure, non_overridable              :: setBoundingBox
    procedure, non_overridable              :: setId
    procedure(init), deferred               :: init
    procedure                               :: kill
    ! Runtime procedures.
    procedure, non_overridable              :: getBoundingBox
    procedure, non_overridable              :: getElementZonesNumber
    procedure, non_overridable              :: getId
    procedure(findElement), deferred        :: findElement
    procedure(distance), deferred           :: distance
    procedure(distanceToNextFace), deferred :: distanceToNextFace
  end type mesh
  
  abstract interface
  
    !! Subroutine 'init'
    !!
    !! Basic description:
    !!   Initialises mesh.
    !!
    !! Arguments:
    !!   folderPath [in] -> Path of the folder where the various mesh files are located.
    !!   name [in]       -> Name of the mesh.
    !!   dict [in]       -> Dictionary with the mesh definition.
    !!
    subroutine init(self, folderPath, dict)
      import                        :: mesh, shortInt, dictionary
      class(mesh), intent(inout)    :: self
      character(*), intent(in)      :: folderPath
      class(dictionary), intent(in) :: dict
      
    end subroutine init
    
    !! Subroutine 'findElement'
    !!
    !! Basic description:
    !!   Returns the index of the tetrahedron occupied by the particle as well as the localId to 
    !!   which the tetrahedron belongs.
    !!
    !! Arguments:
    !!   r [in]           -> Position of the particle.
    !!   u [in]           -> Direction of the particle.
    !!   elementIdx [out] -> Index of the tetrahedron in which the particle is.
    !!   localId [out]    -> Local Id for the given particle.
    !!
    pure subroutine findElement(self, r, u, elementIdx, localId)
      import                                  :: mesh, shortInt, defReal
      class(mesh), intent(in)                 :: self
      real(defReal), dimension(3), intent(in) :: r, u
      integer(shortInt), intent(out)          :: elementIdx, localId

    end subroutine findElement
    
    !! Subroutine 'distance'
    !!
    !! Basic description:
    !!   Returns the distance to the next mesh face intersected by a particle's path.
    !!
    !! Arguments:
    !!   d [out]        -> Distance to the surface intersected by the particle's path.
    !!   coords [inout] -> Coordinates of the particle within the universe (after transformations
    !!                     and with elementIdx already set).
    !!   isInside [out] -> .true. if the particle is inside or entering the mesh. If .false. then 
    !!                     CSG tracking resumes.
    !!
    pure subroutine distance(self, d, coords, isInside)
      import                        :: mesh, coord, defBool, defReal
      class(mesh), intent(in)       :: self
      real(defReal), intent(out)    :: d
      type(coord), intent(inout)    :: coords
      logical(defBool), intent(out) :: isInside

    end subroutine distance

    !! Subroutine 'distanceToNextFace'
    !!
    !! Basic description:
    !!   Returns the distance to the next mesh face intersected by the particle's path. Returns
    !!   INF if the particle does not intersect any face (i.e., if its path is entirely
    !!   contained in the element the particle currently is).
    !!
    !! Arguments:
    !!   d [out]        -> Distance to the next intersected face.
    !!   coords [inout] -> Particle's coordinates.
    !!
    pure subroutine distanceToNextFace(self, d, coords)
      import                     :: mesh, defReal, coord
      class(mesh), intent(in)    :: self
      real(defReal), intent(out) :: d
      type(coord), intent(inout) :: coords

    end subroutine distanceToNextFace

  end interface

contains

  !! Function 'getBoundingBox'
  !!
  !! Basic description:
  !!   Returns the axis-aligned bounding box (AABB) of the mesh.
  !!
  !! Result:
  !!   boundingBox -> AABB of the mesh.
  !!
  pure function getBoundingBox(self) result(boundingBox)
    class(mesh), intent(in)     :: self
    real(defReal), dimension(6) :: boundingBox
    
    boundingBox = self % boundingBox

  end function getBoundingBox

  !! Function 'getElementZonesNumber'
  !!
  !! Basic description:
  !!   Returns the number of element zones in the mesh. This can be (for instance) cell zones in
  !!   OpenFOAM meshes, or other element subdivisions.
  !!
  !! Result:
  !!   nZones: an integer corresponding to the number of element zones in the mesh.
  !!
  elemental function getElementZonesNumber(self) result(nElementZones)
    class(mesh), intent(in) :: self
    integer(shortInt)       :: nElementZones

    nElementZones = self % nElementZones

  end function getElementZonesNumber
  
  !! Subroutine 'setBoundingBox'
  !!
  !! Basic description:
  !!   Sets the axis-aligned bounding box (AABB) of the mesh. Applies NUDGE to prevent the AABB from
  !!   touching any mesh vertex.
  !!
  !! Arguments:
  !!   boundingBox [in] -> An array of six reals whose first three entries correspond to the minimum
  !!                       x-, y-, and z-values of the bounding box and the remaining three
  !!                       correspond to the maximum x-, y- and z-values of the bounding box.
  !!
  pure subroutine setBoundingBox(self, boundingBox)
    class(mesh), intent(inout)              :: self
    real(defReal), dimension(6), intent(in) :: boundingBox
    integer(shortInt)                       :: i
    
    do i = 1, 3
      self % boundingBox(i) = boundingBox(i) - NUDGE
      self % boundingBox(3 + i) = boundingBox(3 + i) + NUDGE

    end do

  end subroutine setBoundingBox
  
  !! Subroutine 'setId'
  !!
  !! Basic description:
  !!   Sets the id of the mesh.
  !!
  !! Arguments:
  !!   id [in] -> Id of the mesh.
  !!
  !! Errors:
  !!   fatalError if id < 1.
  !!
  subroutine setId(self, id)
    class(mesh), intent(inout)    :: self
    integer(shortInt), intent(in) :: id
    character(100), parameter     :: Here = 'setId (mesh_inter.f90)'
    
    ! Catch invalid id and set id.
    if (id < 1) call fatalError(Here, 'Id must be +ve. Is: '//numToChar(id)//'.')
    self % id = id

  end subroutine setId
  
  !! Function 'getId'
  !!
  !! Basic description:
  !!   Returns the id of the mesh.
  !!
  !! Result:
  !!   id -> Id of the mesh.
  !!
  elemental function getId(self) result(id)
    class(mesh), intent(in) :: self
    integer(shortInt)       :: id
    
    id = self % id

  end function getId
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(mesh), intent(inout) :: self
   
    self % id = 0
    self % nElementZones = 0
    self % boundingBox = ZERO
    self % cellZonesFile = .false.
    call self % cellZones % kill()

  end subroutine kill

end module mesh_inter