module mesh_inter
  
  use numPrecision
  use universalVariables,  only : INF, NUDGE
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
  !!   setBoundingBox           -> Sets bounding box of the mesh.
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
  type, public, abstract                          :: mesh
    private
    integer(shortInt), public                     :: id = 0, nElementZones = 0
    real(defReal), dimension(6)                   :: boundingBox = ZERO
    logical(defBool), public                      :: cellZonesFile = .false.
    type(cellZoneShelf), public                   :: cellZones
  contains
    ! Build procedures.
    procedure, non_overridable                    :: setBoundingBox
    procedure, non_overridable                    :: setId
    procedure(init), deferred                     :: init
    procedure                                     :: kill
    ! Runtime procedures.
    procedure, non_overridable                    :: distance
    procedure, non_overridable                    :: distanceToBoundary
    procedure(distanceToBoundaryFace), deferred   :: distanceToBoundaryFace
    procedure(distanceToNextFace), deferred       :: distanceToNextFace
    procedure, non_overridable                    :: findOccupiedElementIdx
    procedure(findElementAndParentIdxs), deferred :: findElementAndParentIdxs
    procedure, non_overridable                    :: getBoundingBox
    procedure, non_overridable                    :: getElementZonesNumber
    procedure, non_overridable                    :: getId
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
    elemental subroutine distanceToBoundaryFace(self, d, coords, parentIdx)
      import                         :: mesh, defReal, coord, shortInt
      class(mesh), intent(in)        :: self
      real(defReal), intent(out)     :: d
      type(coord), intent(inout)     :: coords
      integer(shortInt), intent(out) :: parentIdx

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
    pure subroutine findElementAndParentIdxs(self, r, u, elementIdx, parentIdx)
      import                                  :: mesh, defReal, shortInt  
      class(mesh), intent(in)                 :: self
      real(defReal), dimension(3), intent(in) :: r, u
      integer(shortInt), intent(out)          :: elementIdx, parentIdx

    end subroutine findElementAndParentIdxs

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
  pure subroutine distance(self, d, coords, isInside)
    class(mesh), intent(in)       :: self
    real(defReal), intent(out)    :: d
    type(coord), intent(inout)    :: coords
    logical(defBool), intent(out) :: isInside
    integer(shortInt)             :: parentIdx

    ! Initialise isInside = .true.
    isInside = .true.
    
    ! If particle is already inside a tetrahedron, simply compute the distance to the next mesh face and return.
    if (coords % elementIdx > 0) then
      call self % distanceToNextFace(d, coords)
      return

    end if

    ! If not, we need to check if the particle enters the mesh. If yes, update localId from index of the parent element and return.
    call self % distanceToBoundary(d, coords, parentIdx)
    if (coords % elementIdx > 0) then
      coords % localId = self % cellZones % findCellZone(parentIdx)
      return

    end if
      
    ! If reached here, the particle does not enter the mesh and CSG tracking resumes.
    isInside = .false.

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
  pure subroutine distanceToBoundary(self, d, coords, parentIdx)
    class(mesh), intent(in)        :: self
    real(defReal), intent(out)     :: d
    type(coord), intent(inout)     :: coords
    integer(shortInt), intent(out) :: parentIdx
    integer(shortInt)              :: i
    real(defReal)                  :: rComponent, rEndComponent, lowerBoundingBoxComponent, upperBoundingBoxComponent

    ! Initialise d = INF and check that the particle's path intersects the mesh's bounding box.
    d = INF
    do i = 1, 3
      ! If the particle does not intersect the mesh's bounding box return early.
      rComponent = coords % r(i)
      rEndComponent = coords % rEnd(i)
      lowerBoundingBoxComponent = self % boundingBox(i)
      upperBoundingBoxComponent = self % boundingBox(i + 3)
      if ((rComponent < lowerBoundingBoxComponent .and. rEndComponent < lowerBoundingBoxComponent) .or. &
          (rComponent > upperBoundingBoxComponent .and. rEndComponent > upperBoundingBoxComponent)) return

    end do

    ! If particle intersects the bounding box, compute the distance to the next intersected boundary face.
    call self % distanceToBoundaryFace(d, coords, parentIdx)

  end subroutine distanceToBoundary

  !! Subroutine 'findOccupiedElementIdx'
  !!
  !! Basic description:
  !!   Returns the index of the element occupied by the particle as well as the localId to which the element belongs.
  !!
  !! Arguments:
  !!   r [in]           -> Position of the particle.
  !!   u [in]           -> Direction of the particle.
  !!   elementIdx [out] -> Index of the element in which the particle is.
  !!   localId [out]    -> Local Id for the given particle.
  !!
  pure subroutine findOccupiedElementIdx(self, r, u, elementIdx, localId)
    class(mesh), intent(in)                 :: self
    real(defReal), dimension(3), intent(in) :: r, u
    integer(shortInt), intent(out)          :: elementIdx, localId
    integer(shortInt)                       :: parentIdx

    ! Initialise localId = 1 (corresponds to the particle being in the CSG cell).
    localId = 1

    ! Find indices of the occupied mesh element and its parent element. Update localId only if particle is not 
    ! outside the mesh.
    call self % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    if (elementIdx > 0) localId = self % cellZones % findCellZone(parentIdx)

  end subroutine findOccupiedElementIdx

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
  !!   nElementZones -> Number of element zones in the mesh.
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