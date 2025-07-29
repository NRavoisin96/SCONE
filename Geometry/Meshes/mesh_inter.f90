module mesh_inter
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use dictionary_class,             only : dictionary
  use genericProcedures,            only : fatalError, numToChar, openToRead
  use numPrecision
  use publicObjects,                only : coordData
  use universalVariables,           only : INF, NUDGE
  
  implicit none
  private
  
  ! Extendable methods.
  public :: distanceToBoundary, findHostElement, kill
  
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
  type, public, abstract         :: mesh
    private
    integer(shortInt)            :: id = 0, nLocalIds = 0
    type(axisAlignedBoundingBox) :: boundingBox
  contains
    ! Build procedures.
    procedure(init), deferred               :: init
    procedure, non_overridable              :: initBoundingBox
    procedure                               :: kill
    procedure, non_overridable              :: setId
    procedure, non_overridable              :: setLocalIdsNumber
    procedure, non_overridable              :: setupBase
    ! Runtime procedures.
    procedure                               :: distance
    procedure(distanceToBoundary), deferred :: distanceToBoundary
    procedure(distanceToNextFace), deferred :: distanceToNextFace
    procedure(findHostElement), deferred    :: findHostElement
    procedure, non_overridable              :: getBoundingBoxPtr
    procedure, non_overridable              :: getId
    procedure, non_overridable              :: getLocalIdsNumber
  end type mesh
  
  abstract interface
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
    subroutine distanceToBoundary(self, data)
      import                         :: coordData, mesh
      class(mesh), intent(in)        :: self
      type(coordData), intent(inout) :: data
    end subroutine distanceToBoundary

    !! Subroutine 'distanceToNextFace'
    !!
    !! Basic description:
    !!   Returns the distance to the next intersected mesh face.
    !!
    !! Arguments:
    !!   d [out]        -> Distance to the next intersected face.
    !!   coords [inout] -> Particle's coordinates.
    !!
    subroutine distanceToNextFace(self, data)
      import                         :: coordData, mesh
      class(mesh), intent(in)        :: self
      type(coordData), intent(inout) :: data
    end subroutine distanceToNextFace

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
    subroutine findHostElement(self, data)
      import                         :: coordData, mesh
      class(mesh), intent(in)        :: self
      type(coordData), intent(inout) :: data
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
  subroutine distance(self, data)
    class(mesh), intent(in)        :: self
    type(coordData), intent(inout) :: data

    ! Initialise isInside = .true.
    data % isInside = .true.
    
    ! If particle is already inside an element, simply compute the distance to the next mesh face and return.
    if (data % elementIdx > 0) then
      call self % distanceToNextFace(data)

    else
      ! If not, we need to check if the particle enters the mesh. If yes, update localId from index of the parent element and return.
      call self % distanceToBoundary(data)
      if (data % elementIdx == 0) data % isInside = .false.

    end if

  end subroutine distance

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

  !!
  !!
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
  pure subroutine initBoundingBox(self, allCoords)
    class(mesh), intent(inout)                 :: self
    real(defReal), dimension(:, :), intent(in) :: allCoords
    
    call self % boundingBox % computeBounds(allCoords)

  end subroutine initBoundingBox

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  subroutine kill(self)
    class(mesh), intent(inout) :: self
   
    self % id = 0
    self % nLocalIds = 0
    call self % boundingBox % kill()

  end subroutine kill
  
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
    character(*), parameter       :: Here = 'setId (mesh_inter.f90)'
    
    ! Catch invalid id and set id.
    if (id < 1) call fatalError(Here, 'Id must be positive. Is: '//numToChar(id)//'.')
    self % id = id

  end subroutine setId

  !!
  !!
  !!
  elemental subroutine setLocalIdsNumber(self, nLocalIds)
    class(mesh), intent(inout)    :: self
    integer(shortInt), intent(in) :: nLocalIds

    self % nLocalIds = nLocalIds

  end subroutine

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
    if (id < 1) call fatalError(Here, 'Mesh Id must be positive. Is: '//numToChar(id)//'.')
    self % id = id

  end subroutine setupBase

end module mesh_inter