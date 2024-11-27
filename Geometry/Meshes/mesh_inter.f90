module mesh_inter
  
  use numPrecision
  use universalVariables, only : NUDGE
  use genericProcedures,  only : fatalError, numToChar, openToRead
  use dictionary_class,   only : dictionary
  use coord_class,        only : coord
  
  implicit none
  private
  
  ! Extendable methods.
  public :: kill
  
  !!
  !! Abstract interface for all meshes.
  !!
  !! Mesh represents a subdivision of the entire space into local elements defined by faces 
  !! themselves divided into vertices.
  !!
  !! Private Members:
  !!   meshId          -> Id of the mesh.
  !!   boundingBox     -> Axis-aligned bounding box (AABB) of the mesh.
  !!
  !! Interface:
  !!   checkMeshFolder -> Checks the existence of folder containing the mesh files.
  !!   getBoundingBox  -> Gets bounding box of the mesh.
  !!   getName         -> Gets name of the mesh.
  !!   getNZones       -> Gets the number of element zones in the mesh.
  !!   id              -> Gets Id of the mesh.
  !!   setBoundingBox  -> Sets bounding box of the mesh.
  !!   setId           -> Sets Id of the mesh.
  !!   setName         -> Sets name of the mesh.
  !!   init            -> Initialises mesh from input files.
  !!   kill         -> Returns to uninitialised state.
  !!   findElement     -> Returns local cellId and cellIdx in cellShelf for a given position.
  !!   distance        -> Calculates the distance to the next face in the mesh.
  !!
  type, public, abstract :: mesh
    private
    integer(shortInt)                :: meshId = 0
    real(defReal), dimension(6)      :: boundingBox = ZERO
  contains
    ! Build procedures.
    procedure, non_overridable       :: getBoundingBox
    procedure(getNZones), deferred   :: getNZones
    procedure, non_overridable       :: id
    procedure, non_overridable       :: setBoundingBox
    procedure, non_overridable       :: setId
    procedure(init), deferred        :: init
    procedure                        :: kill
    ! Runtime procedures.
    procedure(findElement), deferred :: findElement
    procedure(distance), deferred    :: distance
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
    
    !! Function 'getNZones'
    !!
    !! Basic description:
    !!   Returns the number of element zones in the mesh. This can be (for instance) cell zones in
    !!   OpenFOAM meshes, or other element subdivisions.
    !!
    !! Result:
    !!   nZones: an integer corresponding to the number of element zones in the mesh.
    !!
    elemental function getNZones(self) result(nZones)
      import                  :: mesh, shortInt
      class(mesh), intent(in) :: self
      integer(shortInt)       :: nZones

    end function getNZones
    
    !! Subroutine 'findElement'
    !!
    !! Basic description:
    !!   Returns the index of the tetrahedron occupied by the particle as well as the localId to 
    !!   which the tetrahedron belongs.
    !!
    !! Arguments:
    !!   r [in]               -> Position of the particle.
    !!   u [in]               -> Direction of the particle.
    !!   tetrahedronIdx [out] -> Index of the tetrahedron in which the particle is.
    !!   localId [out]        -> Local Id for the given particle.
    !!
    !! TODO: Change tetrahedronIdx to elementIdx.
    !!
    pure subroutine findElement(self, r, u, tetrahedronIdx, localId)
      import                                  :: mesh, shortInt, defReal
      class(mesh), intent(in)                 :: self
      real(defReal), dimension(3), intent(in) :: r, u
      integer(shortInt), intent(out)          :: tetrahedronIdx, localId

    end subroutine findElement
    
    !! Subroutine 'distance'
    !!
    !! Basic description:
    !!   Returns the distance to the next mesh face intersected by a particle's path.
    !!
    !! Arguments:
    !!   dist [out]     -> Distance to the surface intersected by the particle's path.
    !!   coords [inout] -> Coordinates of the particle within the universe (after transformations
    !!                     and with tetrahedronIdx already set).
    !!   isInside [out] -> Is the particle inside / entering the mesh? Needed as if this is false
    !!                     then CSG tracking resumes.
    !!
    pure subroutine distance(self, dist, coords, isInside)
      import                        :: mesh, coord, defBool, defReal
      class(mesh), intent(in)       :: self
      real(defReal), intent(out)    :: dist
      type(coord), intent(inout)    :: coords
      logical(defBool), intent(out) :: isInside

    end subroutine distance

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
  !!   Sets the 'meshId' component of the mesh.
  !!
  !! Arguments:
  !!   meshId [in] -> Id of the mesh.
  !!
  !! Errors:
  !!   fatalError if meshId < 1.
  !!
  subroutine setId(self, meshId)
    class(mesh), intent(inout)    :: self
    integer(shortInt), intent(in) :: meshId
    character(100), parameter     :: Here = 'setId (mesh_inter.f90)'
    
    ! Catch invalid id and set id.
    if (meshId < 1) call fatalError(Here, 'Id must be +ve. Is: '//numToChar(meshId)//'.')
    self % meshId = meshId

  end subroutine setId
  
  !! Function 'id'
  !!
  !! Basic description:
  !!   Returns the Id of the mesh.
  !!
  !! Result:
  !!   'id' component of the 'mesh' structure.
  !!
  elemental function id(self)
    class(mesh), intent(in) :: self
    integer(shortInt)       :: id
    
    id = self % meshId

  end function id
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(mesh), intent(inout) :: self
   
    self % meshId  = 0
    self % boundingBox = ZERO

  end subroutine kill

end module mesh_inter