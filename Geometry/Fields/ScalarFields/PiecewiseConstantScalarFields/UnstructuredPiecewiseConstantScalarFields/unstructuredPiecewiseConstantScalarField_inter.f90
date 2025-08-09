module unstructuredPiecewiseConstantScalarField_inter

  use dictionary_class,                   only : dictionary
  use element_class,                      only : elementBox
  use genericProcedures,                  only : fatalError, numToChar
  use geometry_inter,                     only : geometry
  use geometryReg_mod,                    only : geomNum, geomPtr
  use geometryStd_class,                  only : geometryStd, geometryStd_CptrCast
  use intMap_class,                       only : intMap
  use mesh_inter,                         only : mesh
  use numPrecision
  use particle_class,                     only : particle
  use piecewiseConstantScalarField_inter, only : kill_super => kill, piecewiseConstantScalarField
  use unstructuredMesh_inter,             only : getCastUnstructuredMeshPtr, unstructuredMesh

  implicit none
  private

  ! Public procedures.
  public :: init
  
  !!
  !!
  !!
  type, public, abstract, extends(piecewiseConstantScalarField) :: unstructuredPiecewiseConstantScalarField
    private
    class(unstructuredMesh), pointer :: meshPtr => null()
    type(intMap)                     :: activeElementIdxToParentElementIdxMap
  contains
    procedure :: at
    procedure :: init
    procedure :: kill
  end type unstructuredPiecewiseConstantScalarField

contains
  !!
  !!
  !!
  function at(self, p) result(val)
    class(unstructuredPiecewiseConstantScalarField), intent(in) :: self
    class(particle), intent(inout)                              :: p
    integer(shortInt)                                           :: elementIdx
    real(defReal)                                               :: val
    character(*), parameter                                     :: here = 'at (unstructuredPiecewiseConstantScalarField_inter.f90)'

    val = ZERO
    if (.not. associated(self % meshPtr)) call fatalError(here, 'Unassociated mesh pointer.')
    ! First retrieve the element containing the particle.
    elementIdx = p % coords % getLowestElementIdx()
    if (elementIdx == 0) return
    val = self % getValue(self % activeElementIdxToParentElementIdxMap % get(elementIdx))

  end function at

  !!
  !!
  !!
  subroutine init(self, dict)
    class(unstructuredPiecewiseConstantScalarField), intent(inout) :: self
    class(dictionary), intent(in)                                  :: dict
    class(geometry), pointer                                       :: geom
    class(geometryStd), pointer                                    :: geomStd
    class(mesh), pointer                                           :: meshPtr
    class(unstructuredMesh), pointer                               :: unstructuredMeshPtr
    integer(shortInt)                                              :: elementIdx, i, meshId, nActiveElements, nElements, &
                                                                      nGeometries, nParentElements, parentElementIdx
    type(elementBox)                                               :: element
    character(*), parameter :: here = 'init (unstructuredPiecewiseConstantScalarField_inter.f90)'

    ! Retrieve geometry pointer. TODO: Make this more modular.

    ! Check that there is only one geometry in the dictionary,
    nGeometries = geomNum()
    if (nGeometries /= 1) call fatalError(here, 'Geometry registry contains: '//numToChar(nGeometries)//'. Should be 1.')

    geom => geomPtr(1)
    geomStd => geometryStd_CptrCast(geom)

    if (.not. associated(geomStd)) call fatalError(here, 'Geometry is not of type geometryStd.')

    ! Get pointer to unstructured mesh geometry.
    call dict % get(meshId, 'meshId')
    meshPtr => geomStd % getMeshPtr(meshId)
    unstructuredMeshPtr => getCastUnstructuredMeshPtr(meshPtr)

    ! Loop through all elements in the mesh and begin active and parent counts.
    nElements = unstructuredMeshPtr % getElementsNumber()
    nActiveElements = 0
    nParentElements = 0
    do i = 1, nElements
      element = unstructuredMeshPtr % getElementBox(i)
      if (element % ptr % getIsActive()) nActiveElements = nActiveElements + 1
      if (element % ptr % getParentIdx() == 0) nParentElements = nParentElements + 1

    end do
    call self % allocateValues(nParentElements)
    call self % activeElementIdxToParentElementIdxMap % init(nActiveElements)

    ! Loop through all elements and build internal map.
    do i = 1, nElements
      element = unstructuredMeshPtr % getElementBox(i)
      if (element % ptr % getIsActive()) then
        elementIdx = element % ptr % getIdx()
        parentElementIdx = element % ptr % getParentIdx()
        call self % activeElementIdxToParentElementIdxMap % add(elementIdx, &
                                                                merge(elementIdx, parentElementIdx, parentElementIdx == 0))

      end if

    end do

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