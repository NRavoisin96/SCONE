module mesh_inter
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use coord_class,                  only : coord
  use dictionary_class,             only : dictionary
  use genericProcedures,            only : fatalError, numToChar, openToRead
  use numPrecision
  use universalVariables,           only : INF, NUDGE
  
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
  !!   id                       -> Id of the mesh.
  !!   nElementZones            -> Number of element zones in the mesh.
  !!   boundingBox              -> Axis-aligned bounding box (AABB) of the mesh.
  !!   cellZonesFile            -> .true. if a file is present to explicitly assign element zones.
  !!   cellZones                -> Shelf that stores element zones.
  !!
  !! Interface:
  !!   getBoundingBox           -> Returns the bounding box of the mesh.
  !!   getElementZonesNumber    -> Returns the number of element zones in the mesh.
  !!   getId                    -> Returns the id of the mesh.
  !!   setId                    -> Sets Id of the mesh.
  !!   init                     -> Initialises mesh from input files.
  !!   kill                     -> Returns to uninitialised state.
  !!   distance                 -> Calculates the distance travelled by a particle within the mesh.
  !!   distanceToBoundary       -> Calculates the distance to the mesh boundary.
  !!   distanceToBoundaryFace   -> Calculates the distance to the intersected mesh boundary face.
  !!   distanceToNextFace       -> Calculates the distance to the next face in the mesh.
  !!   findOccupiedElementIdx   -> Finds the index of the mesh element occupied by a particle and
  !!                               the corresponding local id.
  !!   findElementAndParentIdxs -> Finds the index of the mesh element occupied by a particle and
  !!                               the index of the parent element of the occupied element.
  !!
  type, public, abstract :: mesh
    private
    integer(shortInt)            :: id = 0, nLocalIds = 0
    type(axisAlignedBoundingBox) :: boundingBox
  contains
    ! Build procedures.
    procedure, non_overridable                  :: setBoundingBox
    procedure, non_overridable                  :: setId
    procedure, non_overridable                  :: setLocalIdsNumber
    procedure, non_overridable                  :: setupBase
    procedure(init), deferred                   :: init
    procedure                                   :: kill
    ! Runtime procedures.
    procedure, non_overridable                  :: distance
    procedure, non_overridable                  :: distanceToBoundary
    procedure(distanceToBoundaryFace), deferred :: distanceToBoundaryFace
    procedure(distanceToNextFace), deferred     :: distanceToNextFace
    procedure(findHostElement), deferred        :: findHostElement
    procedure, non_overridable                  :: getBoundingBoxPtr
    procedure, non_overridable                  :: getLocalIdsNumber
    procedure, non_overridable                  :: getId
  end type mesh
  
  abstract interface

    !! Subroutine 'distanceToNextFace'
    !!
    !! Basic description:
    !!   Returns the distance to the next intersected mesh face.
    !!
    !! Arguments:
    !!   d [out]        -> Distance to the next intersected face.
    !!   coords [inout] -> Particle's coordinates.
    !!
    elemental subroutine distanceToNextFace(self, d, coords)
      import                     :: mesh, defReal, coord
      class(mesh), intent(in)    :: self
      real(defReal), intent(out) :: d
      type(coord), intent(inout) :: coords

    end subroutine distanceToNextFace

    !! Subroutine 'distanceToBoundaryFace'
    !!
    !! Basic description:
    !!   Returns the distance to the mesh boundary face intersected by a particle's path. Also
    !!   returns the index of the parent element containing the intersected boundary face.
    !!
    !! Arguments:
    !!   d [out]         -> Distance to the mesh boundary face.
    !!   coords [inout]  -> Particle's coordinates.
    !!   parentIdx [out] -> Index of the parent element containing the boundary face.
    !!
    subroutine distanceToBoundaryFace(self, d, coords)
      import                         :: mesh, defReal, coord
      class(mesh), intent(in)        :: self
      real(defReal), intent(out)     :: d
      type(coord), intent(inout)     :: coords

    end subroutine distanceToBoundaryFace

    !! Subroutine 'findElementAndParentIdxs'
    !!
    !! Basic description:
    !!   Returns the index of the mesh element occupied by a particle. Also returns the index
    !!   of the parent mesh element containing the occupied element.
    !!
    !! Arguments:
    !!   r [in]           -> Particle's location.
    !!   u [in]           -> Particle's direction.
    !!   elementIdx [out] -> Index of the mesh element occupied by the particle.
    !!   parentIdx [out]  -> Index of the parent mesh element containing the occupied element.
    !!
    subroutine findHostElement(self, coords)
      import                                 :: coord, mesh 
      class(mesh), intent(in)                :: self
      type(coord), intent(inout)             :: coords

    end subroutine findHostElement

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

  end interface

contains

  !! Subroutine 'distance'
  !!
  !! Basic description:
  !!   Returns the distance to the next mesh face intersected by a particle's path.
  !!
  !! Arguments:
  !!   d [out]        -> Distance to the surface intersected by the particle's path.
  !!   coords [inout] -> Coordinates of the particle within the universe (after transformations and with elementIdx already set).
  !!   isInside [out] -> .true. if the particle is inside or entering the mesh. If .false. then CSG tracking resumes.
  !!
  subroutine distance(self, d, coords, isInside)
    class(mesh), intent(in)       :: self
    real(defReal), intent(out)    :: d
    type(coord), intent(inout)    :: coords
    logical(defBool), intent(out) :: isInside
    integer(shortInt)             :: elementIdx

    ! Initialise isInside = .true.
    isInside = .true.
    
    ! If particle is already inside a tetrahedron, simply compute the distance to the next mesh face and return.
    if (coords % getElementIdx() > 0) then
      call self % distanceToNextFace(d, coords)

    else
      ! If not, we need to check if the particle enters the mesh. If yes, update localId from index of the parent element and return.
      call self % distanceToBoundary(d, coords)
      if(coords % getElementIdx() == 0) isInside = .false.

    end if

  end subroutine distance

  !! Subroutine 'distanceToBoundary'
  !!
  !! Basic description:
  !!   Returns the distance to the next intersected mesh boundary face.
  !!
  !! Arguments:
  !!   d [out]         -> Distance to the intersected mesh boundary face.
  !!   coords [inout]  -> Particle's coordinates.
  !!   parentIdx [out] -> Index of the parent element containing the intersected mesh boundary face.
  !!
  subroutine distanceToBoundary(self, d, coords)
    class(mesh), intent(in)        :: self
    real(defReal), intent(out)     :: d
    type(coord), intent(inout)     :: coords
    real(defReal), dimension(3)    :: r, rEnd
    real(defReal), dimension(6)    :: bounds

    ! Initialise d = INF and check that the particle's path intersects the mesh's bounding box.
    d = INF
    r = coords % getPosition()
    rEnd = coords % getEndPosition()
    bounds = self % boundingBox % getBounds()
    if (any(r < bounds(1:3) .and. rEnd < bounds(1:3)) .or. any(r > bounds(4:6) .and. rEnd > bounds(4:6))) return

    ! If particle intersects the bounding box, compute the distance to the next intersected boundary face.
    call self % distanceToBoundaryFace(d, coords)

  end subroutine distanceToBoundary

  !! Function 'getBoundingBox'
  !!
  !! Basic description:
  !!   Returns the axis-aligned bounding box (AABB) of the mesh.
  !!
  !! Result:
  !!   boundingBox -> AABB of the mesh.
  !!
  function getBoundingBoxPtr(self) result(boundingBoxPtr)
    class(mesh), target, intent(in)       :: self
    type(axisAlignedBoundingBox), pointer :: boundingBoxPtr
    
    boundingBoxPtr => self % boundingBox

  end function getBoundingBoxPtr

  !! Function 'getElementZonesNumber'
  !!
  !! Basic description:
  !!   Returns the number of element zones in the mesh. This can be (for instance) cell zones in
  !!   OpenFOAM meshes, or other element subdivisions.
  !!
  !! Result:
  !!   nElementZones -> Number of element zones in the mesh.
  !!
  elemental function getLocalIdsNumber(self) result(nLocalIds)
    class(mesh), intent(in) :: self
    integer(shortInt)       :: nLocalIds

    nLocalIds = self % nLocalIds

  end function getLocalIdsNumber
  
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
    class(mesh), intent(inout)               :: self
    type(axisAlignedBoundingBox), intent(in) :: boundingBox
    
    self % boundingBox = boundingBox

  end subroutine setBoundingBox

  !! Subroutine 'setElementZonesNumber'
  !!
  !! Basic description:
  !!   Sets the number of element zones in the mesh.
  !!
  !! Arguments:
  !!   nElementZones [in] -> Number of element zones in the mesh.
  !!
  elemental subroutine setLocalIdsNumber(self, nLocalIds)
    class(mesh), intent(inout)    :: self
    integer(shortInt), intent(in) :: nLocalIds

    self % nLocalIds = nLocalIds

  end subroutine setLocalIdsNumber
  
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

  !! Subroutine 'setupBase'
  !!
  !! Basic description:
  !!   Sets basic mesh components from dictionary.
  !!
  !! Arguments:
  !!   dict [in] -> A dictionary.
  !!
  !! Errors:
  !!   - fatalError if id < 1.
  !!
  subroutine setupBase(self, dict)
    class(mesh), intent(inout)    :: self
    class(dictionary), intent(in) :: dict
    integer(shortInt)             :: id
    character(*), parameter       :: here = 'setupBase (mesh_inter.f90)'

    ! Load id from the dictionary. Call fatal error if id is unvalid.
    call dict % get(id, 'id')
    if (id < 1) call fatalError(Here, 'Mesh Id must be +ve. Is: '//numToChar(id)//'.')
    self % id = id

  end subroutine setupBase
  
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
    self % nLocalIds = 0
    call self % boundingBox % kill()

  end subroutine kill

end module mesh_inter