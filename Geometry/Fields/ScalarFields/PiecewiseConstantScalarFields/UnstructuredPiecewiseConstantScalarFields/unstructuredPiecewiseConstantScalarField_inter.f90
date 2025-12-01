module unstructuredPiecewiseConstantScalarField_inter

  use coordList_class,                    only : coordList
  use dictionary_class,                   only : dictionary
  use element_class,                      only : elementBox
  use field_inter,                        only : field
  use genericProcedures,                  only : fatalError, numToChar
  use geometry_inter,                     only : geometry
  use geometryReg_mod,                    only : geomNum, geomPtr
  use intMap_class,                       only : intMap
  use mesh_inter,                         only : mesh
  use numPrecision
  use piecewiseConstantScalarField_inter, only : kill_super => kill, piecewiseConstantScalarField
  use universalVariables,                 only : INF, NOT_PRESENT
  use unstructuredMesh_inter,             only : getCastUnstructuredMeshPtr, unstructuredMesh

  implicit none
  private

  ! Public procedures.
  public :: getCastUnstructuredPiecewiseConstantScalarField, init
  
  !!
  !!
  !!
  type, public, abstract, extends(piecewiseConstantScalarField) :: unstructuredPiecewiseConstantScalarField
    private
    class(unstructuredMesh), pointer :: meshPtr => null()
    type(intMap)                     :: activeElementIdxToParentElementIdxMap
  contains
    procedure                           :: at
    procedure                           :: createExtremalMaterialValues
    procedure                           :: init
    procedure                           :: kill
    procedure(retrieveValues), deferred :: retrieveValues
  end type unstructuredPiecewiseConstantScalarField

  abstract interface
    !!
    !!
    !!
    subroutine retrieveValues(self, dict)
      import                                                         :: dictionary, unstructuredPiecewiseConstantScalarField
      class(unstructuredPiecewiseConstantScalarField), intent(inout) :: self
      class(dictionary), intent(in)                                  :: dict
    end subroutine retrieveValues

  end interface

contains
  !!
  !!
  !!
  function at(self, defaultValue, coords, mult) result(val)
    class(unstructuredPiecewiseConstantScalarField), intent(in) :: self
    real(defReal), intent(in)                                   :: defaultValue
    type(coordList), intent(in)                                 :: coords
    real(defReal), intent(in), optional                         :: mult
    integer(shortInt)                                           :: elementIdx
    real(defReal)                                               :: val
    character(*), parameter                                     :: here = 'at (unstructuredPiecewiseConstantScalarField_inter.f90)'

    if (.not. associated(self % meshPtr)) call fatalError(here, 'Unassociated mesh pointer.')
    
    ! First retrieve the element containing the particle.
    elementIdx = coords % getLowestElementIdx()
    if (elementIdx == 0) then
      val = defaultValue
      return

    end if
    val = self % getValue(self % activeElementIdxToParentElementIdxMap % get(elementIdx))
    if (present(mult)) val = val * mult

  end function at

  !!
  !!
  !!
  subroutine createExtremalMaterialValues(self, maximumMaterialValues, minimumMaterialValues, map)
    class(unstructuredPiecewiseConstantScalarField), intent(in) :: self
    real(defReal), dimension(:), allocatable, intent(out)       :: maximumMaterialValues, minimumMaterialValues
    type(intMap), intent(out)                                   :: map
    integer(shortInt)                                           :: i, idx, materialIdx, nElements, nMaterials
    integer(shortInt), dimension(:), allocatable                :: localIdsToMaterialIdxs
    real(defReal)                                               :: value
    type(elementBox)                                            :: element

    ! Loop through all parent elements and create map linking indices to material indices.
    nElements = self % meshPtr % getElementsNumber()
    nMaterials = 0
    localIdsToMaterialIdxs = self % meshPtr % getLocalIdsToMaterialIdxs()
    do i = 1, nElements
      element = self % meshPtr % getElementBox(i)
      if (element % ptr % getParentIdx() == 0) then
        materialIdx = localIdsToMaterialIdxs(element % ptr % getLocalId())
        if (map % getOrDefault(materialIdx, NOT_PRESENT) == NOT_PRESENT) then
          nMaterials = nMaterials + 1
          call map % add(materialIdx, nMaterials)

        end if

      end if

    end do

    ! Loop through all parent elements and find maximum field value for all materials.
    allocate(maximumMaterialValues(nMaterials), minimumMaterialValues(nMaterials))
    maximumMaterialValues = -INF
    minimumMaterialValues = INF
    do i = 1, nElements
      element = self % meshPtr % getElementBox(i)
      if (element % ptr % getParentIdx() == 0) then
        ! The current element is a parent element. Retrieve its localId and update maximum value for its material.
        idx = map % get(localIdsToMaterialIdxs(element % ptr % getLocalId()))
        value = self % getValue(element % ptr % getIdx())
        maximumMaterialValues(idx) = max(maximumMaterialValues(idx), value)
        minimumMaterialValues(idx) = min(minimumMaterialValues(idx), value)

      end if

    end do

  end subroutine createExtremalMaterialValues

  !!
  !!
  !!
  function getCastUnstructuredPiecewiseConstantScalarField(source) result(ptr)
    class(field), intent(in)                                 :: source
    class(unstructuredPiecewiseConstantScalarField), pointer :: ptr

    select type(temp => source)
      class is(unstructuredPiecewiseConstantScalarField)
        ptr => temp

      class default
        ptr => null()

    end select

  end function getCastUnstructuredPiecewiseConstantScalarField

  !!
  !!
  !!
  subroutine init(self, dict)
    class(unstructuredPiecewiseConstantScalarField), intent(inout) :: self
    class(dictionary), intent(in)                                  :: dict
    class(geometry), pointer                                       :: geom
    class(mesh), pointer                                           :: meshPtr
    integer(shortInt)                                              :: elementIdx, i, meshId, nActiveElements, nElements, &
                                                                      nGeometries, nParentElements, parentElementIdx
    type(elementBox)                                               :: element
    character(*), parameter :: here = 'init (unstructuredPiecewiseConstantScalarField_inter.f90)'

    ! Retrieve geometry pointer. TODO: Make this more modular.

    ! Check that there is only one geometry in the dictionary,
    nGeometries = geomNum()
    if (nGeometries /= 1) call fatalError(here, 'Geometry registry contains: '//numToChar(nGeometries)//'. Should be 1.')

    ! Get pointer to unstructured mesh geometry.
    geom => geomPtr(1)
    call dict % get(meshId, 'meshId')
    meshPtr => geom % getMeshPtr(meshId)
    self % meshPtr => getCastUnstructuredMeshPtr(meshPtr)
    if (.not. associated(self % meshPtr)) call fatalError(here, 'Unable to retrieve unstructured mesh pointer.')

    ! Loop through all elements in the mesh and begin active and parent counts.
    nElements = self % meshPtr % getElementsNumber()
    nActiveElements = 0
    nParentElements = 0
    do i = 1, nElements
      element = self % meshPtr % getElementBox(i)
      nActiveElements = nActiveElements + merge(1, 0, element % ptr % getIsActive())
      nParentElements = nParentElements + merge(1, 0, element % ptr % getParentIdx() == 0)

    end do
    call self % allocateValues(nParentElements)
    call self % activeElementIdxToParentElementIdxMap % init(nActiveElements)

    ! Loop through all elements and build internal map.
    do i = 1, nElements
      element = self % meshPtr % getElementBox(i)
      if (element % ptr % getIsActive()) then
        elementIdx = element % ptr % getIdx()
        parentElementIdx = element % ptr % getParentIdx()
        call self % activeElementIdxToParentElementIdxMap % add(elementIdx, &
                                                                merge(elementIdx, parentElementIdx, parentElementIdx == 0))

      end if

    end do

    ! Retrieve field values.
    call self % retrieveValues(dict)

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(unstructuredPiecewiseConstantScalarField), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % meshPtr => null()
    call self % activeElementIdxToParentElementIdxMap % kill()

  end subroutine kill

end module unstructuredPiecewiseConstantScalarField_inter