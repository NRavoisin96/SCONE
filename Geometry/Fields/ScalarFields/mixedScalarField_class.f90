module mixedScalarField_class

  use coordList_class,                                only : coordList
  use dictionary_class,                               only : dictionary
  use errors_mod,                                     only : fatalError
  use geometry_inter,                                 only : geometry
  use geometryReg_mod,                                only : geomNum, geomPtr
  use intMap_class,                                   only : intMap
  use materialMenu_mod,                               only : nMat
  use numPrecision
  use OpenFOAMScalarField_class,                      only : OpenFOAMScalarField
  use scalarField_inter,                              only : kill_super => kill, scalarField
  use universalVariables,                             only : NOT_PRESENT
  use universalVariables,                             only : INF
  use unstructuredPiecewiseConstantScalarField_inter, only : unstructuredPiecewiseConstantScalarField

  implicit none
  private

  !!
  !!
  !!
  type, public :: unstructuredPiecewiseConstantScalarFieldBox
    class(unstructuredPiecewiseConstantScalarField), pointer :: meshFieldPtr => null()
  end type unstructuredPiecewiseConstantScalarFieldBox

  !!
  !!
  !!
  type, public, extends(scalarField) :: mixedScalarField
    private
    class(scalarField), pointer                                                   :: CSGFieldPtr => null()
    class(unstructuredPiecewiseConstantScalarFieldBox), dimension(:), allocatable :: meshFields
    type(intMap)                                                                  :: meshIdxsToMeshFieldIdxs
  contains
    procedure          :: at
    procedure, private :: computeMaximumAndMinimumMaterialValues
    procedure          :: init
    procedure, private :: initMeshFields
    procedure          :: kill
    procedure          :: setValues
  end type mixedScalarField

contains
  !!
  !!
  !!
  function at(self, coords) result(val)
    class(mixedScalarField), intent(in) :: self
    class(coordList), intent(in)        :: coords
    integer(shortInt)                   :: meshIdx, meshFieldIdx
    real(defReal)                       :: val

    val = ZERO
    meshIdx = coords % getLowestMeshIdx()
    if (0 < meshIdx) then
      meshFieldIdx = self % meshIdxsToMeshFieldIdxs % getOrDefault(meshIdx, NOT_PRESENT)
      if (meshFieldIdx /= NOT_PRESENT) val = self % meshFields(meshFieldIdx) % meshFieldPtr % at(coords)

    elseif (associated(self % CSGFieldPtr)) then
      val = self % CSGFieldPtr % at(coords)

    end if

  end function at

  !!
  !!
  !!
  subroutine computeMaximumAndMinimumMaterialValues(self)
    class(mixedScalarField), intent(inout)   :: self
    integer(shortInt)                        :: i, nMaterials
    real(defReal), dimension(:), allocatable :: maximumMaterialValues, minimumMaterialValues

    ! Query all fields and find maximum values for materials.
    nMaterials = nMat()
    allocate(maximumMaterialValues(nMaterials), minimumMaterialValues(nMaterials))
    maximumMaterialValues = -INF
    minimumMaterialValues = INF
    if (associated(self % CSGFieldPtr)) then

    end if

    if (allocated(self % meshFields)) then
      do i = 1, size(self % meshFields)
        maximumMaterialValues = max(maximumMaterialValues, self % meshFields(i) % meshFieldPtr % getMaximumMaterialValues())
        minimumMaterialValues = min(minimumMaterialValues, self % meshFields(i) % meshFieldPtr % getMinimumMaterialValues())

      end do

    end if
    call self % setMaximumMaterialValues(maximumMaterialValues)
    call self % setMinimumMaterialValues(minimumMaterialValues)

  end subroutine computeMaximumAndMinimumMaterialValues

  !!
  !!
  !!
  subroutine init(self, dict)
    class(mixedScalarField), intent(inout)   :: self
    class(dictionary), intent(in)            :: dict
    class(geometry), pointer                 :: geometryPtr
    character(*), parameter                  :: here = 'init (mixedScalarField_class.f90)'

    ! Check that there is only one geometry and call fatalError if not.
    if (1 < geomNum()) call fatalError(here, 'More than one geometry.')
    geometryPtr => geomPtr(1)
    if (.not. associated(geometryPtr)) call fatalError(here, 'Unable to retrieve geometry pointer.')

    ! Initialise CSG field from dictionary.

    ! Initialise mesh fields from dictionary.
    if (dict % isPresent('meshFields')) call self % initMeshFields(dict % getDictPtr('meshFields'), geometryPtr)

    ! Compute maximum and minimum material values.
    call self % computeMaximumAndMinimumMaterialValues()

  end subroutine init

  !!
  !!
  !!
  subroutine initMeshFields(self, dict, geom)
    class(mixedScalarField), intent(inout)        :: self
    class(dictionary), intent(in)                 :: dict
    class(geometry), intent(in)                   :: geom
    class(dictionary), pointer                    :: meshFieldDict
    character(nameLen)                            :: meshFieldType
    character(nameLen), dimension(:), allocatable :: meshNames
    integer(shortInt)                             :: i, nMeshFields
    character(*), parameter :: here = 'initMeshFields (mixedPiecewiseConstantScalarField_class.f90)'

    ! Get the names of all the meshes that have a field defined.
    call dict % keys(meshNames, 'dict')
    
    ! Compute the number of mesh fields and allocate space.
    nMeshFields = 0
    if (allocated(meshNames)) nMeshFields = size(meshNames)

    if (nMeshFields == 0) return
    allocate(self % meshFields(nMeshFields))
    call self % meshIdxsToMeshFieldIdxs % init(nMeshFields)

    ! Loop through each mesh and build its unstructured field.
    do i = 1, nMeshFields
      meshFieldDict => dict % getDictPtr(meshNames(i))
      call meshFieldDict % get(meshFieldType, 'type')

      select case(meshFieldType)
        case('OpenFOAMScalarField')
          allocate(OpenFOAMScalarField :: self % meshFields(i) % meshFieldPtr)

        case default
          call fatalError(here, 'Invalid unstructured field type.')

      end select

      ! Initialise unstructured field.
      call self % meshFields(i) % meshFieldPtr % init(meshFieldDict)

      ! Get index of mesh in the geometry and add it to the map.
      call self % meshIdxsToMeshFieldIdxs % add(geom % getMeshIdxByName(meshNames(i)), i)

    end do

  end subroutine initMeshFields

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(mixedScalarField), intent(inout) :: self
    integer(shortInt)                      :: i

    ! Superclass.
    call kill_super(self)

    ! Local.
    if (associated(self % CSGFieldPtr)) then
      call self % CSGFieldPtr % kill()
      deallocate(self % CSGFieldPtr)

    end if

    if (allocated(self % meshFields)) then
      do i = 1, size(self % meshFields)
        if (associated(self % meshFields(i) % meshFieldPtr)) then
          call self % meshFields(i) % meshFieldPtr % kill()
          deallocate(self % meshFields(i) % meshFieldPtr)

        end if

      end do
      deallocate(self % meshFields)

    end if

    call self % meshIdxsToMeshFieldIdxs % kill()

  end subroutine kill

  !!
  !!
  !!
  subroutine setValues(self, values)
    class(mixedScalarField), intent(inout)  :: self
    real(defReal), dimension(:), intent(in) :: values
    integer(shortInt)                       :: i

    ! Set values in CSG field first.
    if (associated(self % CSGFieldPtr)) then

    end if

    ! Set values in mesh fields.
    if (allocated(self % meshFields)) then
      do i = 1, size(self % meshFields)
        call self % meshFields(i) % meshFieldPtr % setValues(values)

      end do

    end if

    ! Recompute maximum and minimum material values.
    call self % computeMaximumAndMinimumMaterialValues()

  end subroutine setValues

end module mixedScalarField_class