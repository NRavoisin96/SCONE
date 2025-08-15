module mesh_inter
  
  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use charMap_class,                only : charMap
  use dictionary_class,             only : dictionary
  use genericProcedures,            only : fatalError, numToChar, openToRead
  use numPrecision
  use publicObjects,                only : coordData
  use RNG_class,                    only : RNG
  use universalVariables
  
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
    integer(shortInt)                            :: id = 0, nLocalIds = 0
    integer(shortInt), dimension(:), allocatable :: localIdsToMaterialIdxs
    type(axisAlignedBoundingBox)                 :: boundingBox
  contains
    ! Build procedures.
    procedure(init), deferred                       :: init
    procedure, non_overridable                      :: initBoundingBox
    procedure                                       :: kill
    procedure, non_overridable                      :: setId
    procedure, non_overridable                      :: setLocalIdsNumber
    procedure, non_overridable                      :: setLocalIdsToMaterialIdxs
    procedure, non_overridable                      :: setupBase
    ! Runtime procedures.
    procedure                                       :: distance
    procedure(distanceToBoundary), deferred         :: distanceToBoundary
    procedure(distanceToNextFace), deferred         :: distanceToNextFace
    procedure(explicitBoundaryConditions), deferred :: explicitBoundaryConditions
    procedure(findHostElement), deferred            :: findHostElement
    procedure, non_overridable                      :: getBoundingBoxBounds
    procedure, non_overridable                      :: getBoundingBoxPtr
    procedure(getFaceBoundaryConditions), deferred  :: getFaceBoundaryConditions
    procedure(getFaceIsBoundary), deferred          :: getFaceIsBoundary
    procedure(getFacesNumber), deferred             :: getFacesNumber
    procedure(getElementIsActive), deferred         :: getElementIsActive
    procedure(getElementsNumber), deferred          :: getElementsNumber
    procedure(getElementVolume), deferred           :: getElementVolume
    procedure, non_overridable                      :: getId
    procedure, non_overridable                      :: getLocalIdsNumber
    procedure, non_overridable                      :: getLocalIdsToMaterialIdxs
    procedure(getUniqueIdOffset), deferred          :: getUniqueIdOffset
    procedure(sampleInitialPosition), deferred      :: sampleInitialPosition
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

    !!
    !!
    !!
    subroutine explicitBoundaryConditions(self, idx, boundaryConditionType, data)
      import                         :: coordData, mesh, shortInt
      class(mesh), intent(in)        :: self
      integer(shortInt), intent(in)  :: idx, boundaryConditionType
      type(coordData), intent(inout) :: data
    end subroutine explicitBoundaryConditions

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

    !!
    !!
    !!
    function getFaceBoundaryConditions(self, idx) result(boundaryConditions)
      import                                   :: mesh, N_BC_TYPES, shortInt
      class(mesh), intent(in)                  :: self
      integer(shortInt), intent(in)            :: idx
      integer(shortInt), dimension(N_BC_TYPES) :: boundaryConditions
    end function getFaceBoundaryConditions

    !!
    !!
    !!
    function getFaceIsBoundary(self, idx) result(isBoundary)
      import                        :: defBool, mesh, shortInt
      class(mesh), intent(in)       :: self
      integer(shortInt), intent(in) :: idx
      logical(defBool)              :: isBoundary
    end function getFaceIsBoundary

    !!
    !!
    !!
    function getFacesNumber(self, activeOnly) result(nFaces)
      import                                 :: defBool, mesh, shortInt
      class(mesh), intent(in)                :: self
      logical(defBool), intent(in), optional :: activeOnly
      integer(shortInt)                      :: nFaces
    end function getFacesNumber

    !!
    !!
    !!
    function getElementIsActive(self, idx) result(isActive)
      import                        :: defBool, mesh, shortInt
      class(mesh), intent(in)       :: self
      integer(shortInt), intent(in) :: idx
      logical(defBool)              :: isActive
    end function getElementIsActive

    !!
    !!
    !!
    function getElementsNumber(self, activeOnly) result(nElements)
      import                                 :: defBool, mesh, shortInt
      class(mesh), intent(in)                :: self
      logical(defBool), intent(in), optional :: activeOnly
      integer(shortInt)                      :: nElements
    end function getElementsNumber

    !!
    !!
    !!
    function getElementVolume(self, idx) result(volume)
      import                        :: defReal, mesh, shortInt
      class(mesh), intent(in)       :: self
      integer(shortInt), intent(in) :: idx
      real(defReal)                 :: volume
    end function getElementVolume

    !!
    !!
    !!
    function getUniqueIdOffset(self) result(uniqueIdOffset)
      import                  :: mesh, shortInt
      class(mesh), intent(in) :: self
      integer(shortInt)       :: uniqueIdOffset
    end function getUniqueIdOffset

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
    subroutine init(self, folderPath, dict, materialsMap)
      import                        :: charMap, dictionary, mesh, shortInt
      class(mesh), intent(inout)    :: self
      character(*), intent(in)      :: folderPath
      class(dictionary), intent(in) :: dict
      type(charMap), intent(in)     :: materialsMap
    end subroutine init

    !!
    !!
    !!
    subroutine sampleInitialPosition(self, elementIdx, rand, localId, r)
      import                                   :: defReal, mesh, RNG, shortInt
      class(mesh), intent(in)                  :: self
      integer(shortInt), intent(in)            :: elementIdx
      class(RNG), intent(inout)                :: rand
      integer(shortInt), intent(out)           :: localId
      real(defReal), dimension(3), intent(out) :: r
    end subroutine sampleInitialPosition

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
    if (0 < data % elementIdx) then
      call self % distanceToNextFace(data)

    else
      ! If not, we need to check if the particle enters the mesh. If yes, update localId from index of the parent element and return.
      call self % distanceToBoundary(data)
      if (data % elementIdx == 0) data % isInside = .false.

    end if

  end subroutine distance

  !!
  !!
  !!
  pure function getBoundingBoxBounds(self) result(boundingBoxBounds)
    class(mesh), intent(in)        :: self
    real(defReal), dimension(3, 2) :: boundingBoxBounds

    boundingBoxBounds = self % boundingBox % getBounds()

  end function getBoundingBoxBounds

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

  !!
  !!
  !!
  pure function getLocalIdsToMaterialIdxs(self) result(localIdsToMaterialIdxs)
    class(mesh), intent(in)                      :: self
    integer(shortInt), dimension(:), allocatable :: localIdsToMaterialIdxs

    localIdsToMaterialIdxs = self % localIdsToMaterialIdxs

  end function getLocalIdsToMaterialIdxs

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
    if (allocated(self % localIdsToMaterialIdxs)) deallocate(self % localIdsToMaterialIdxs)
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

  !!
  !!
  !!
  pure subroutine setLocalIdsToMaterialIdxs(self, localIdsToMaterialIdxs)
    class(mesh), intent(inout)                  :: self
    integer(shortInt), dimension(:), intent(in) :: localIdsToMaterialIdxs

    self % localIdsToMaterialIdxs = localIdsToMaterialIdxs

  end subroutine setLocalIdsToMaterialIdxs

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