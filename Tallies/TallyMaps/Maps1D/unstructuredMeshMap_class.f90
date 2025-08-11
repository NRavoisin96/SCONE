module unstructuredMeshMap_class

  use dictionary_class,        only : dictionary
  use element_class,           only : elementBox
  use genericProcedures,       only : fatalError, numToChar
  use geometry_inter,          only : geometry
  use geometryReg_mod,         only : geomNum, geomPtr
  use intMap_class,            only : intMap
  use mesh_inter,              only : mesh
  use numPrecision
  use outputFile_class,        only : outputFile
  use particle_class,          only : particleState
  use tallyMap1D_inter,        only : kill_super => kill, tallyMap1D
  use unstructuredMesh_inter,  only : getCastUnstructuredMeshPtr, unstructuredMesh

  implicit none
  private

  type, public, extends(tallyMap1D) :: unstructuredMeshMap
    private
    integer(shortInt)                            :: nBins = 0
    integer(shortInt), dimension(:), allocatable :: parentElementBinsToIdxs
    real(defReal), dimension(:), allocatable     :: binVolumes
    type(intMap)                                 :: activeElementIdxsToParentElementBinsMap
  contains
    procedure :: bins
    procedure :: getAxisName
    procedure :: getBinVolume
    procedure :: init
    procedure :: kill
    procedure :: map
    procedure :: print
  end type unstructuredMeshMap

contains

  !!
  !!
  !!
  elemental function bins(self, D) result(N)
    class(unstructuredMeshMap), intent(in) :: self
    integer(shortInt), intent(in)          :: D
    integer(shortInt)                      :: N

    N = merge(self % nBins, 0, any([0, 1] == D))

  end function bins

  !!
  !!
  !!
  function getAxisName(self) result(name)
    class(unstructuredMeshMap), intent(in) :: self
    character(nameLen)                     :: name

    name = 'UnstructuredMeshElement'

  end function getAxisName

  !!
  !!
  !!
  elemental function getBinVolume(self, idx) result(binVolume)
    class(unstructuredMeshMap), intent(in) :: self
    integer(shortInt), intent(in)          :: idx
    real(defReal)                          :: binVolume

    binVolume = self % binVolumes(idx)

  end function getBinVolume

  !!
  !!
  !!
  subroutine init(self, dict)
    class(unstructuredMeshMap), intent(inout) :: self
    class(dictionary), intent(in)             :: dict
    integer(shortInt)                         :: elementIdx, i, meshId, nActiveElements, nElements, nGeometries, &
                                                 nParentElements, parentElementIdx
    class(geometry), pointer                  :: geom
    class(mesh), pointer                      :: meshPtr
    class(unstructuredMesh), pointer          :: unstructuredMeshPtr
    type(elementBox)                          :: element
    type(intMap)                              :: parentElementIdxToBinMap
    character(*), parameter                   :: here = 'init (unstructuredMeshMap_class.f90)'

    ! Retrieve geometry pointer. TODO: Make this more modular.

    ! Check that there is only one geometry in the dictionary,
    nGeometries = geomNum()
    if (nGeometries /= 1) call fatalError(here, 'Geometry registry contains: '//numToChar(nGeometries)//'. Should be 1.')

    ! Get pointer to unstructured mesh geometry.
    geom => geomPtr(1)
    call dict % get(meshId, 'meshId')
    meshPtr => geom % getMeshPtr(meshId)
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
    self % nBins = nParentElements

    ! Allocate memory then begin assigning indices.
    allocate(self % binVolumes(nParentElements), self % parentElementBinsToIdxs(nParentElements))
    call self % activeElementIdxsToParentElementBinsMap % init(nActiveElements)
    call parentElementIdxToBinMap % init(nParentElements)

    ! Do a pass to populate maps.
    nParentElements = 0
    do i = 1, nElements
      element = unstructuredMeshPtr % getElementBox(i)
      elementIdx = element % ptr % getIdx()
      parentElementIdx = element % ptr % getParentIdx()
      
      if (parentElementIdx == 0) then
        nParentElements = nParentElements + 1
        call parentElementIdxToBinMap % add(elementIdx, nParentElements)
        self % binVolumes(nParentElements) = element % ptr % getVolume()
        self % parentElementBinsToIdxs(nParentElements) = elementIdx

      end if

      if (element % ptr % getIsActive()) then
        call self % activeElementIdxsToParentElementBinsMap % add(elementIdx, &
        parentElementIdxToBinMap % get(merge(elementIdx, parentElementIdx, parentElementIdx == 0)))

      end if

    end do

    ! Kill temporary map.
    call parentElementIdxToBinMap % kill()

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(unstructuredMeshMap), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % nBins = 0
    if (allocated(self % parentElementBinsToIdxs)) deallocate(self % parentElementBinsToIdxs)
    if (allocated(self % binVolumes)) deallocate(self % binVolumes)
    call self % activeElementIdxsToParentElementBinsMap % kill()

  end subroutine kill

  !!
  !!
  !!
  elemental function map(self, state) result(idx)
    class(unstructuredMeshMap), intent(in) :: self
    class(particleState), intent(in)       :: state
    integer(shortInt)                      :: idx

    idx = self % activeElementIdxsToParentElementBinsMap % getOrDefault(state % elementIdx, 0)

  end function map

  !!
  !!
  !!
  subroutine print(self, out)
    class(unstructuredMeshMap), intent(in) :: self
    class(outputFile), intent(inout)       :: out
    character(nameLen)                     :: name
    integer(shortInt)                      :: i

    ! Name the array.
    name = trim(self % getAxisName()) // 'Bins'

    call out % startArray(name, [1, self % nBins])

    ! Print element indices.
    do i = 1, self % nBins
      call out % addValue(numToChar(self % parentElementBinsToIdxs(i)))

    end do

    call out % endArray()

  end subroutine print

end module unstructuredMeshMap_class