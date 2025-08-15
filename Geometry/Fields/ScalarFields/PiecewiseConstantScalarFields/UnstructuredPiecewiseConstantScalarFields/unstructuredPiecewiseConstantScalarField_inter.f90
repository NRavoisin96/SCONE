module unstructuredPiecewiseConstantScalarField_inter

  use coordList_class,                    only : coordList
  use dictionary_class,                   only : dictionary
  use element_class,                      only : elementBox
  use field_inter,                        only : field
  use genericProcedures,                  only : fatalError, numToChar
  use geometry_inter,                     only : geometry
  use geometryReg_mod,                    only : geomNum, geomPtr
  use intMap_class,                       only : intMap
  use materialMenu_mod,                   only : nMat
  use mesh_inter,                         only : mesh
  use numPrecision
  use piecewiseConstantScalarField_inter, only : kill_super => kill, piecewiseConstantScalarField
  use universalVariables,                 only : INF
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
  function at(self, coords) result(val)
    class(unstructuredPiecewiseConstantScalarField), intent(in) :: self
    class(coordList), intent(in)                                :: coords
    integer(shortInt)                                           :: elementIdx
    real(defReal)                                               :: val
    character(*), parameter                                     :: here = 'at (unstructuredPiecewiseConstantScalarField_inter.f90)'

    val = ZERO
    if (.not. associated(self % meshPtr)) call fatalError(here, 'Unassociated mesh pointer.')
    ! First retrieve the element containing the particle.
    elementIdx = coords % getLowestElementIdx()
    if (elementIdx == 0) return
    val = self % getValue(self % activeElementIdxToParentElementIdxMap % get(elementIdx))

  end function at

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
    integer(shortInt)                                              :: elementIdx, i, materialIdx, meshId, nActiveElements, &
                                                                      nElements, nGeometries, nMaterials, nParentElements, &
                                                                      parentElementIdx
    integer(shortInt), dimension(:), allocatable                   :: localIdsToMaterialIdxs
    real(defReal), dimension(:), allocatable                       :: maximumMaterialValues, minimumMaterialValues
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
      if (element % ptr % getIsActive()) nActiveElements = nActiveElements + 1
      if (element % ptr % getParentIdx() == 0) nParentElements = nParentElements + 1

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

    ! Loop through all parent elements and find maximum field value for all materials.
    localIdsToMaterialIdxs = self % meshPtr % getLocalIdsToMaterialIdxs()
    nMaterials = nMat()
    allocate(maximumMaterialValues(nMaterials), minimumMaterialValues(nMaterials))
    maximumMaterialValues = -INF
    minimumMaterialValues = INF
    do i = 1, nElements
      element = self % meshPtr % getElementBox(i)
      if (element % ptr % getParentIdx() == 0) then
        ! The current element is a parent element. Retrieve its localId and update maximum value for its material.
        materialIdx = localIdsToMaterialIdxs(element % ptr % getLocalId())
        maximumMaterialValues(materialIdx) = max(maximumMaterialValues(materialIdx), self % getValue(element % ptr % getIdx()))
        minimumMaterialValues(materialIdx) = min(minimumMaterialValues(materialIdx), self % getValue(element % ptr % getIdx()))

      end if

    end do
    call self % setMaximumMaterialValues(maximumMaterialValues)
    call self % setMinimumMaterialValues(minimumMaterialValues)

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