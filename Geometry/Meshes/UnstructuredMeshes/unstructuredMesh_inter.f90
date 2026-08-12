module unstructuredMesh_inter

  use accelerationStructure_inter,       only : accelerationStructure
  use accelerationStructureFactory_func, only : newAccelerationStructurePtr
  use axisAlignedBoundingBox_class,      only : axisAlignedBoundingBox
  use charMap_class,                     only : charMap
  use dictionary_class,                  only : dictionary
  use edge_class,                        only : edgeBox
  use element_class,                     only : buildElementPayload, element, elementBox, inclusionTestResult, &
                                                elementIntersectionTestResult, newElementIntersectionTestPayload
  use extentTopologicalObject_inter,     only : buildExtentTopologicalObjectPayload
  use face_class,                        only : buildFacePayload, face, faceBox
  use genericProcedures,                 only : append, fatalError, numToChar, areEqual
  use mesh_inter,                        only : mesh, kill_super => kill
  use numPrecision
  use publicObjects,                     only : basicEdgeInfo, basicElementInfo, basicFaceInfo, basicVertexInfo, &
                                                coordData, intersectionTestPayload, intersectionTestResult, &
                                                meshBoundaryConditionInfo, meshLocalIdInfo, newCoordData, &
                                                newIntersectionTestPayload
  use RNG_class,                         only : RNG
  use topologicalObject_inter,           only : topologicalObjectBox
  use topologicalObjectShelf_class,      only : topologicalObjectShelf
  use triangulationFactory_func,         only : newTriangulationPtr
  use triangulationMethod_inter,         only : triangulationMethod
  use universalVariables
  use vertex_class,                      only : buildVertexPayload, vertexBox
  use ratint
  use face_class

  implicit none
  private

  ! Public procedures.
  public :: distanceToBoundary, distanceToNextFace, findHostElement, getCastUnstructuredMeshPtr, kill
  
  !! Abstract interface to group all unstructured meshes. An unstructured mesh uses a vertex -> face 
  !! -> element representation of space. Each element is composed by a set of faces which are themselves 
  !! composed by a number of vertices. Elements can be grouped together into zones. This is useful to 
  !! assign material filling to mesh elements. Local ids are assigned in the order of the cell zone 
  !! definition.
  !!
  !! Public members:
  !!   cellZones                -> Shelf that stores cell zones.
  !!   edges                    -> Shelf that stores edges.
  !!   elements                 -> Shelf that stores elements.
  !!   faces                    -> Shelf that stores faces.
  !!   vertices                 -> Shelf that stores vertices.
  !!   nVertices                -> Number of vertices in the mesh.
  !!   nFaces                   -> Number of faces in the mesh.
  !!   nEdges                   -> Number of edges in the mesh.
  !!   nElements                -> Number of elements in the mesh.
  !!   nInternalFaces           -> Number of internal faces in the mesh.
  !!
  !! Interface:
  !!   kill                     -> Returns to an unitialised state.
  !!   printComposition         -> Displays mesh composition to the user.
  !!   distanceToBoundaryFace   -> Checks if a particle enters the mesh and returns distance to entry 
  !!                               intersection.
  !!   distanceToNextFace       -> Returns the distance to the next mesh face.
  !!   findElementAndParentIdxs -> Returns the index of the mesh element occupied by a particle. Also
  !!                               returns the index of the parent mesh element containing the occupied
  !!                               element.
  !!
  type, public, abstract, extends(mesh)   :: unstructuredMesh
    private
    class(accelerationStructure), pointer :: acceleration
    class(triangulationMethod), pointer   :: triangulation
    integer(shortInt)                     :: nVertices = 0, nEdges = 0, nInternalFaces = 0
    type(topologicalObjectShelf)          :: edges, elements, faces, vertices
  contains
    ! Build procedures.
    procedure                                        :: assignBoundaryConditions
    procedure                                        :: assignLocalIds
    procedure(buildBoundaryConditionInfos), deferred :: buildBoundaryConditionInfos
    procedure(buildLocalIdInfos), deferred           :: buildLocalIdInfos
    procedure                                        :: importBoundaryConditionsFromFile
    procedure                                        :: importLocalIdsFromFile
    procedure(importMesh), deferred                  :: importMesh
    procedure                                        :: init
    procedure                                        :: initEdgeShelf
    procedure                                        :: initElementShelf
    procedure                                        :: initFaceShelf
    procedure                                        :: initVertexShelf
    procedure                                        :: kill
    procedure, non_overridable                       :: printComposition
    procedure                                        :: setEdgesNumber
    procedure                                        :: setInternalFacesNumber
    procedure                                        :: setVerticesNumber
    ! Runtime procedures.
    procedure                                        :: distanceToBoundary
    procedure                                        :: distanceToNextFace
    procedure                                        :: explicitBoundaryConditions
    procedure                                        :: findHostElement
    procedure                                        :: getEdgesNumber
    procedure                                        :: getElementBox
    procedure                                        :: getElementIsActive
    procedure                                        :: getElementsNumber
    procedure                                        :: getElementVolume
    procedure                                        :: getFaceBoundaryConditions
    procedure                                        :: getFaceIsBoundary
    procedure                                        :: getFacesNumber
    procedure                                        :: getInternalFacesNumber
    procedure                                        :: getParentElementsNumber
    procedure                                        :: getUniqueIdOffset
    procedure                                        :: getVerticesNumber
    procedure                                        :: sampleInitialPosition
  end type unstructuredMesh

  abstract interface
    !!
    !!
    !!
    function buildBoundaryConditionInfos(self, assignmentMethod) result(boundaryConditionInfos)
      import                                                     :: meshBoundaryConditionInfo, nameLen, unstructuredMesh
      class(unstructuredMesh), intent(in)                        :: self
      character(nameLen), intent(in)                             :: assignmentMethod
      type(meshBoundaryConditionInfo), dimension(:), allocatable :: boundaryConditionInfos
    end function buildBoundaryConditionInfos

    !!
    !!
    !!
    function buildLocalIdInfos(self, assignmentMethod) result(localIdInfos)
      import                                           :: meshLocalIdInfo, nameLen, unstructuredMesh
      class(unstructuredMesh), intent(in)              :: self
      character(nameLen), intent(in)                   :: assignmentMethod
      type(meshLocalIdInfo), dimension(:), allocatable :: localIdInfos
    end function buildLocalIdInfos

    !! Subroutine 'distanceToNextFace'
    !!
    !! Basic description:
    !!   Returns the distance to the next intersected mesh face.
    !!
    !! Arguments:
    !!   d [out]        -> Distance to the next intersected face.
    !!   coords [inout] -> Particle's coordinates.
    !!
    subroutine importMesh(self, folderPath)
      import                                 :: unstructuredMesh
      class(unstructuredMesh), intent(inout) :: self
      character(*), intent(in)               :: folderPath
    end subroutine importMesh

  end interface

contains
  !!
  !!
  !!
  subroutine assignBoundaryConditions(self, dict)
    class(unstructuredMesh), intent(inout)                     :: self
    class(dictionary), intent(in)                              :: dict
    character(nameLen)                                         :: assignmentMethod
    character(pathLen)                                         :: path
    integer(shortInt)                                          :: fixedIdx, i, j, nBoundaryInfos, nBoundaryFaces, nFaces, &
                                                                  nFixedTemperatureBoundaryConditions, &
                                                                  nFixedTemperatureValues, nTemperatureBoundaryConditions, &
                                                                  nTransportBoundaryConditions
    integer(shortInt), dimension(:), allocatable               :: temperatureBoundaryConditions, tempFaceIdxs, &
                                                                  transportBoundaryConditions
    real(defReal), dimension(:), allocatable                   :: fixedTemperatureValues
    type(faceBox)                                              :: box
    type(meshBoundaryConditionInfo), dimension(:), allocatable :: boundaryConditionInfos
    character(*), parameter                                    :: here = 'assignBoundaryConditions (unstructuredMesh_inter.f90)'

    ! Get assignment method from dictionary.
    call dict % get(assignmentMethod, 'assignmentMethod')

    ! Apply logic depending on specific assignment method.
    select case(assignmentMethod)
      case('all')
        allocate(boundaryConditionInfos(1))
        if (dict % isPresent('transportBCs')) &
        call dict % get(boundaryConditionInfos(1) % boundaryConditions(TRANSPORT_BCs), 'transportBCs')
        
        if (dict % isPresent('temperatureBCs')) then
          call dict % get(boundaryConditionInfos(1) % boundaryConditions(TEMPERATURE_BCs), 'temperatureBCs')
          if (boundaryConditionInfos(1) % boundaryConditions(TEMPERATURE_BCs) == FIXED_TEMPERATURE_BC) &
          call dict % get(boundaryConditionInfos(1) % boundaryValues(TEMPERATURE_BCs), 'fixedTemperatureValues')

        end if

        ! Count the number of boundary faces in the mesh.
        nBoundaryFaces = 0
        nFaces = self % faces % getObjectsNumber()
        allocate(boundaryConditionInfos(1) % faceIdxs(nFaces))
        do i =  1, nFaces
          box = self % faces % getFaceBox(i)
          if (box % ptr % getIsBoundary()) then
            nBoundaryFaces = nBoundaryFaces + 1
            boundaryConditionInfos(1) % faceIdxs(nBoundaryFaces) = box % ptr % getIdx()

          end if

        end do
        
        ! Resize.
        allocate(tempFaceIdxs(nBoundaryFaces))
        tempFaceIdxs = boundaryConditionInfos(1) % faceIdxs(1:nBoundaryFaces)
        call move_alloc(tempFaceIdxs, boundaryConditionInfos(1) % faceIdxs)

      case('fileBased')
        call dict % get(path, 'path')
        boundaryConditionInfos = self % importBoundaryConditionsFromFile(path)

      case default
        boundaryConditionInfos = self % buildBoundaryConditionInfos(assignmentMethod)
        nBoundaryInfos = size(boundaryConditionInfos)
        
        ! Get boundary from dictionary.
        if (dict % isPresent('transportBCs')) then
          call dict % get(transportBoundaryConditions, 'transportBCs')
          nTransportBoundaryConditions = size(transportBoundaryConditions)

          if (nTransportBoundaryConditions /= nBoundaryInfos) &
          call fatalError(here, 'Number of transport boundary conditions: '//numToChar(nTransportBoundaryConditions)//&
                                ' does not match the number of boundaries: '//numToChar(nBoundaryInfos)//'.')

          do i = 1, nBoundaryInfos
            boundaryConditionInfos(i) % boundaryConditions(TRANSPORT_BCs) = transportBoundaryConditions(i)

          end do

        end if

        if (dict % isPresent('temperatureBCs')) then
          call dict % get(temperatureBoundaryConditions, 'temperatureBCs')
          nTemperatureBoundaryConditions = size(temperatureBoundaryConditions)

          if (dict % isPresent('fixedTemperatureValues')) call dict % get(fixedTemperatureValues, 'fixedTemperatureValues')
          nFixedTemperatureValues = 0
          if (allocated(fixedTemperatureValues)) nFixedTemperatureValues = size(fixedTemperatureValues)

          if (nTemperatureBoundaryConditions /= nBoundaryInfos) &
          call fatalError(here, 'Number of temperature boundary conditions: '//numToChar(nTemperatureBoundaryConditions)//&
                                ' does not match the number of boundaries: '//numToChar(nBoundaryInfos)//'.')

          nFixedTemperatureBoundaryConditions = 0
          do i = 1, nBoundaryInfos
            boundaryConditionInfos(i) % boundaryConditions(TEMPERATURE_BCs) = temperatureBoundaryConditions(i)
            if (boundaryConditionInfos(i) % boundaryConditions(TEMPERATURE_BCs) == FIXED_TEMPERATURE_BC) &
            nFixedTemperatureBoundaryConditions = nFixedTemperatureBoundaryConditions + 1

          end do

          ! Check that the number of Dirichlet boundary conditions matches the number of input values.
          if (nFixedTemperatureBoundaryConditions /= nFixedTemperatureValues) &
          call fatalError(here, 'Number of Dirichlet temperature boundary conditions does not match the number of fixed values.')

          fixedIdx = 0
          do i = 1, nBoundaryInfos
            if (boundaryConditionInfos(i) % boundaryConditions(TEMPERATURE_BCs) == FIXED_TEMPERATURE_BC) then
              fixedIdx = fixedIdx + 1
              boundaryConditionInfos(i) % boundaryValues(TEMPERATURE_BCs) = fixedTemperatureValues(fixedIdx)

            end if

          end do

        end if

    end select

    ! Set boundary conditions for all faces in the shelf.
    do i = 1, size(boundaryConditionInfos)
      do j = 1, size(boundaryConditionInfos(i) % faceIdxs)
        box = self % faces % getFaceBox(boundaryConditionInfos(i) % faceIdxs(j))
        call box % ptr % setBoundaryConditions(boundaryConditionInfos(i))

      end do

    end do

  end subroutine assignBoundaryConditions

  !!
  !!
  !!
  subroutine assignLocalIds(self, dict, materialsMap)
    class(unstructuredMesh), intent(inout)           :: self
    class(dictionary), intent(in)                    :: dict
    type(charMap), intent(in)                        :: materialsMap
    character(nameLen)                               :: assignmentMethod
    character(nameLen), dimension(:), allocatable    :: fillNames
    character(pathLen)                               :: path
    class(dictionary), pointer                       :: localIdsDict
    integer(shortInt)                                :: i, j, materialIdx, nLocalIds
    integer(shortInt), dimension(:), allocatable     :: localIdsToMaterialIdxs
    type(elementBox)                                 :: element
    type(meshLocalIdInfo), dimension(:), allocatable :: localIdInfos
    character(*), parameter                          :: here = 'assignLocalIds (unstructuredMesh_inter.f90)'

    ! Check if localIdsDict is present. Assign only one localIdInfos if not.
    assignmentMethod = 'all'
    if (dict % isPresent('localIds')) then
      localIdsDict => dict % getDictPtr('localIds')
      call localIdsDict % get(assignmentMethod, 'assignmentMethod')

    end if

    ! Apply logic depending on specific assignment method.
    select case(assignmentMethod)
      case('all')
        ! Assign one localIdInfos.
        allocate(localIdInfos(1))
        localIdInfos(1) % localId = 1
        localIdInfos(1) % elementIdxs = [(i, i = 1, self % elements % getObjectsNumber())]

      case('fileBased')
        call localIdsDict % get(path, 'path')
        localIdInfos = self % importLocalIdsFromFile(path)

      case default
        localIdInfos = self % buildLocalIdInfos(assignmentMethod)

    end select

    ! Assign localIds to all elements in the mesh.
    nLocalIds = size(localIdInfos)
    call self % setLocalIdsNumber(nLocalIds)
    do i = 1, nLocalIds
      do j = 1, size(localIdInfos(i) % elementIdxs)
        element = self % elements % getElementBox(localIdInfos(i) % elementIdxs(j))
        call element % ptr % setLocalId(localIdInfos(i) % localId)

      end do

    end do

    ! Now create map linking each global id to a material fill.
    if (.not. dict % isPresent('fills')) call fatalError(here, 'Missing fills.')
    call dict % get(fillNames, 'fills')
    if (size(fillNames) /= nLocalIds) call fatalError(here, 'Mismatch between number of localIds and material fills.')
    allocate(localIdsToMaterialIdxs(nLocalIds))
    
    do i = 1, nLocalIds
      materialIdx = materialsMap % getOrDefault(fillNames(i), NOT_PRESENT)
      if (materialIdx == NOT_PRESENT) call fatalError(here, 'Unknown materal: '//trim(fillNames(i))//'.')
      localIdsToMaterialIdxs(i) = materialIdx

    end do
    call self % setLocalIdsToMaterialIdxs(localIdsToMaterialIdxs)

  end subroutine assignLocalIds

  !! Subroutine 'distanceToBoundaryFace'
  !!
  !! Basic description:
  !!   Returns the distance to the mesh boundary face intersected by a particle's path. Also returns the index
  !!   of the parent element containing the intersected boundary face.
  !!
  !! See mesh_inter for details.
  !!
  subroutine distanceToBoundary(self, data)
    class(unstructuredMesh), intent(in)                   :: self
    type(coordData), intent(inout)                        :: data
    type(axisAlignedBoundingBox), pointer                 :: boundingBoxPtr
    type(intersectionTestResult)                          :: boundingBoxIntersectionResult
    type(faceBox)                                         :: boundaryFace
    type(topologicalObjectBox), dimension(:), allocatable :: faceElements
    type(coordData)                                       :: endData
    type(inclusionTestResult)                             :: testResult
    logical(defBool)                                      :: overshoot
    real(defReal)                                         :: d
    character(*), parameter                               :: here = 'distanceToBoundary (unstructuredMesh_inter.f90)'
    
    ! Initialise parentIdx = 0, edgeIdx = 0 and vertexIdx = 0 then search the tree for the intersected boundary face.
    boundingBoxPtr => self % getBoundingBoxPtr()
    call boundingBoxPtr % intersects(newIntersectionTestPayload(data % r, data % u, data % dMax), boundingBoxIntersectionResult)
    if (.not. boundingBoxIntersectionResult % intersects) return

    call self % acceleration % findEntranceBoundaryFace(self % faces, data, boundaryFace)
    if (.not. associated(boundaryFace % ptr)) return

    ! Retrieve the element associated with the boundary face.
    faceElements = boundaryFace % ptr % getSharingElements()
    ! Downcast elements to correct type.
    select type(ptr => faceElements(1) % ptr)
      type is(element)
        ! Set coords % endPosition to the minimum computed distance plus a slight forward nudge.
        endData = newCoordData(data % r + (data % d + NUDGE) * data % u, data % u)

        ! If the element associated with the intersected face does not contain the end position, begin rescue.
        testResult = ptr % isPointInside(endData % r)
        if (testResult % status == INSIDE_ELEMENT) then
          data % elementIdx = ptr % getIdx()
          data % faceIdx = boundaryFace % ptr % getIdx()
          data % localId = ptr % getLocalId()

        else
          call self % findHostElement(endData)
          ! Update distance.
          data % elementIdx = endData % elementIdx
          data % localId = endData % localId
          
          overshoot = data % elementIdx == 0
          d = norm2(endData % r - data % r)
          data % d = merge(INF, d, overshoot)
          data % u = merge(data % u, (endData % r - data % r) / d, overshoot)

        end if

      class default
        call fatalError(here, 'Element with index: '//numToChar(ptr % getIdx())//' is not an element.')

    end select

  end subroutine distanceToBoundary


  ! subroutine distanceToNextFace(self, data)
  !   class(unstructuredMesh), intent(in)                   :: self
  !   type(coordData), intent(inout)                        :: data
  !   type(elementBox)                                      :: currentElement
  !   type(elementIntersectionTestResult)                   :: intersectionResult
  !   integer(shortInt)                                     :: i, nElements
  !   type(topologicalObjectBox), dimension(:), allocatable :: faceElements
  !   type(faceBox) :: currentFace
  !   type(orientatedFaceBox), dimension(:), allocatable :: elementFaces
  !   type(inclusionTestResult) :: testIn
  !   real(defReal), dimension(3) :: end
  !   character(*), parameter                               :: here = 'distanceToNextFace (unstructuredMesh_inter.f90)'
    
  !   ! Retrieve the element currently occupied by the particle and compute potential 
  !   ! face intersections.
  !   currentElement = self % elements % getElementBox(data % elementIdx)
  !   call currentElement%ptr%&
  !     intersects(newElementIntersectionTestPayload(data%r,data%u,data%dMax,.true., .false.,data%currentFaceIdxs,data%front),&
  !                                          intersectionResult)
  !   if (.not. intersectionResult % intersects) return
  !   data % d = intersectionResult % d
  !   data % faceIdx = intersectionResult % intersectedFace % ptr % getIdx()
    
  !   ! If the intersected face is a boundary face then the particle is leaving the mesh.
  !   if (intersectionResult % intersectedFace % ptr % getIsBoundary()) then
  !     data % elementIdx = 0
  !     data % localId = 0

  !   else
  !     ! Else, retrieve the elements sharing the intersected face from mesh connectivity then
  !     ! update elementIdx and localId.
  !     faceElements = intersectionResult % intersectedFace % ptr % getSharingElements()
  !     nElements = size(faceElements)
  !     if (nElements /= 2) &
  !     call fatalError(here, 'Internal face: '//numToChar(intersectionResult % intersectedFace % ptr % getIdx())// &
  !                           ' is not associated to the correct number of elements.')

  !     do i = 1, 2
  !       ! Downcast element to correct type.
  !       select type(ptr => faceElements(i) % ptr)
  !         type is(element)
  !           if (.not. associated(currentElement % ptr, ptr)) then
  !             ! We have found our new element.
  !             data % elementIdx = ptr % getIdx()
  !             data % localId = ptr % getLocalId()

  !           end if

  !         class default
  !           call fatalError(here, 'Element with index: '//numToChar(ptr % getIdx())//' is not an element.')

  !       end select

  !     end do

  !   end if

  ! end subroutine distanceToNextFace



  !! Subroutine 'distanceToNextFace'
  !!
  !! Basic description:
  !!   Returns the distance to the next face intersected by the particle's path. Returns INF if the particle
  !!   does not intersect any face (i.e., if its path is entirely contained in the element the particle 
  !!   currently is). Algorithm adapted from Macpherson, et al. (2009). DOI: 10.1002/cnm.1128.
  !!
  !! See mesh_inter for details.
  !!
  subroutine distanceToNextFace(self, data)
    class(unstructuredMesh), intent(in)                   :: self
    type(coordData), intent(inout)                        :: data
    type(elementBox)                                      :: currentElement
    type(elementIntersectionTestResult)                   :: intersectionResult, intersectionResultNew
    integer(shortInt)                                     :: i, j, nElements, currentIdx, k
    type(faceBox) :: currentFace
    type(topologicalObjectBox), dimension(:), allocatable :: faceElements, edgeFaces, edgeElements, vertexElements, vertexFaces
    type(orientatedFaceBox), dimension(:), allocatable :: elementFaces, tempElementFaces
    type(edgeBox), dimension(:), allocatable :: faceEdges
    type(vertexBox), dimension(:), allocatable :: faceVertices
    type(edgeBox) :: currentEdge
    type(vertexBox) ::currentVertex
    type(inclusionTestResult) :: testIn
    real(defReal), dimension(3) :: endPt
    integer(shortInt), dimension(:), allocatable :: faceIDs

    logical :: result, validElem
    
    character(*), parameter                               :: here = 'distanceToNextFace (unstructuredMesh_inter.f90)'
    
    ! Retrieve the element currently occupied by the particle and compute potential 
    ! face intersections.
    currentElement = self % elements % getElementBox(data % elementIdx)
    call currentElement%ptr%&
      intersects(newElementIntersectionTestPayload(data%r,data%u,data%dMax,.true., .false.,data%currentFaceIdxs,data%front),&
                                           intersectionResult)

    if (.not. intersectionResult % intersects) return
    data % d = intersectionResult % d
    data % faceIdx = intersectionResult % intersectedFace % ptr % getIdx()
    data % currentFaceIdxs = intersectionResult%currentFaceIdxs
    data%front = intersectionResult%front

    endPt = data%r + (data%u * data%dMax)


    ! If the intersected face is a boundary face then the particle is leaving the mesh.
    if (intersectionResult % intersectedFace % ptr % getIsBoundary()) then
      data % elementIdx = 0
      data % localId = 0
      return

    else

      if (intersectionResult%front <= 1) then 
        ! Else, retrieve the elements sharing the intersected face from mesh connectivity then
        ! update elementIdx and localId.
        faceElements = intersectionResult % intersectedFace % ptr % getSharingElements()
        nElements = size(faceElements)
        if (nElements /= 2) &
        call fatalError(here, 'Internal face: '//numToChar(intersectionResult % intersectedFace % ptr % getIdx())// &
                              ' is not associated to the correct number of elements.')

        do i = 1, 2
          ! Downcast element to correct type.
          select type(ptr => faceElements(i) % ptr)
            type is(element)
              testIn = ptr%hybridIsPointInside(data % r + (data%u * data%dMax))
              if (.not. associated(currentElement % ptr, ptr) .and. &
                  (testIn%status == INSIDE_ELEMENT .or. testIn%status == ON_BOUNDARY_ELEMENT)) then !! comment this out and tests fail
                ! We have found our new element.
                data % elementIdx = ptr % getIdx()
                data % localId = ptr % getLocalId()
                return

              end if

            class default
              call fatalError(here, 'Element with index: '//numToChar(ptr % getIdx())//' is not an element.')

          end select

        end do
      else 
        elementFaces = currentElement%ptr%getOrientatedFaces()
        if (intersectionResult%front == 2) then 

          faceEdges = intersectionResult%intersectedFace%ptr%getEdges()

          do i=1, size(faceEdges)
            if (areEqual(faceEdges(i)%ptr%distanceSquared(intersectionResult%intersectionPt), ZERO)) then 
              currentEdge = faceEdges(i)
              exit 
            end if 
          end do
          if (.not. associated(currentEdge%ptr)) then 
            return 
          end if

          edgeElements = currentEdge%ptr%getSharingElements() 
          edgeFaces = currentEdge%ptr%getSharingFaces() 

          !print *, size(edgeElements)
          elemLoop: do i=1, size(edgeElements)
            validElem = .true.
            select type(ptr => edgeElements(i) % ptr)
              type is(element)
                if (associated(ptr, currentElement%ptr)) then 
                  cycle elemLoop 
                end if
                !tempElementFaces = ptr%getOrientatedFaces()
                testIn = ptr%hybridIsPointInsideTPO(data%r, data%u, edgeFaces, size(edgeFaces))
                !print *, testIn%status == OUTSIDE_ELEMENT
                if (testIn%status == INSIDE_ELEMENT) then 
                  data % elementIdx = ptr % getIdx()
                  data % localId = ptr % getLocalId()
                  return 
                end if


              class default
                call fatalError(here, 'Element with index: '//numToChar(ptr % getIdx())//' is not an element.')
            end select
          end do elemLoop

          
        else 
          faceVertices = intersectionResult%intersectedFace%ptr%getVertices()
          do i=1, size(faceVertices)
            if (areEqual(faceVertices(i)%ptr%distanceSquared(intersectionResult%intersectionPt), ZERO)) then 
              currentVertex = faceVertices(i)
              exit 
            end if 
          end do

          vertexElements = currentVertex%ptr%getSharingElements() 
          vertexFaces = currentVertex%ptr%getSharingFaces() 

          elemLoopV: do i=1, size(vertexElements)
            validElem = .true.
            select type(ptr => vertexElements(i) % ptr)
              type is(element)
                if (associated(ptr, currentElement%ptr)) then 
                  cycle elemLoopV 
                end if


                testIn = ptr%hybridIsPointInsideTPO(data%r, data%u, vertexFaces, size(vertexFaces))

                if (testIn%status == INSIDE_ELEMENT) then 
                  data % elementIdx = ptr % getIdx()
                  data % localId = ptr % getLocalId()
                  return 
                end if

              class default
                call fatalError(here, 'Element with index: '//numToChar(ptr % getIdx())//' is not an element.')
            end select

          end do elemLoopV

        end if
      end if
    end if
      

    !If none of the previous cases have been triggered, the particle is pointing outside of the mesh
    data % elementIdx = 0
    data % localId = 0

  

  end subroutine distanceToNextFace



                ! tempElementFaces = ptr%getOrientatedFaces()

                ! elemFaceLoop: do j=1, size(tempElementFaces)
                !   edgeFaceLoop: do k=1, size(edgeFaces)

                !     select type(ptrF => edgeFaces(k) % ptr)
                !       type is(face)
                !         if (ptrF%getIdx()==tempElementFaces(j)%face%ptr%getIdx()) then

                !           if (dot_product(ptrF%getNormal(), data%u) > 0) then 
                !             validElem = .false.
                !             cycle elemLoop 
                !           end if

                !         end if 
                !     end select

                !   end do edgeFaceLoop
                ! end do elemFaceLoop
                ! if (validElem) then 
                !   data % elementIdx = ptr % getIdx()
                !   data % localId = ptr % getLocalId()
                ! end if


  
  !!
  !!
  !!
  subroutine explicitBoundaryConditions(self, idx, boundaryConditionType, data)
    class(unstructuredMesh), intent(in)                   :: self
    integer(shortInt), intent(in)                         :: idx, boundaryConditionType
    type(coordData), intent(inout)                        :: data
    integer(shortInt), dimension(N_BC_TYPES)              :: faceBoundaryConditions
    type(faceBox)                                         :: box
    type(topologicalObjectBox), dimension(:), allocatable :: sharingElements
    character(*), parameter                               :: here = 'explicitBoundaryConditions (unstructuredMesh_inter.f90)'

    ! Select logic to use based on specific boundary condition.
    box = self % faces % getFaceBox(idx)
    faceBoundaryConditions = box % ptr % getBoundaryConditions()
    select case(boundaryConditionType)
      case(TRANSPORT_BCs)
        select case(faceBoundaryConditions(TRANSPORT_BCs))
          case(REFLECTIVE_BC)
            call box % ptr % flipDirection(data % u)
            sharingElements = box % ptr % getSharingElements()
            if (1 < size(sharingElements)) call fatalError(here, 'Boundary face is associated with more than one element.')

            ! Downcast pointer to correct type.
            select type(ptr => sharingElements(1) % ptr)
              type is(element)
                ! Set index and localId.
                data % elementIdx = ptr % getIdx()
                data % localId = ptr % getLocalId()

              class default
                call fatalError(here, 'Element with index: '//numToChar(ptr % getIdx())//' is not an element.')

            end select

          case default
            call fatalError(here, &
            'Unsupported transport boundary condition: '//numToChar(faceBoundaryConditions(TRANSPORT_BCs))//'.')

        end select

      case(TEMPERATURE_BCs)
        ! Do nothing for now.

      case default
        call fatalError(here, 'Invalid boundary condition type: '//numToChar(boundaryConditionType)//'.')

    end select

  end subroutine explicitBoundaryConditions

  !! Subroutine 'findElementAndParentIdxs'
  !!
  !! Basic description:
  !!   Returns the index of the mesh element occupied by a particle. Also returns the index of the parent mesh
  !!   element containing the occupied element.
  !!
  !! See mesh_inter for details.
  !!
  subroutine findHostElement(self, data)
    class(unstructuredMesh), intent(in)   :: self
    type(coordData), intent(inout)        :: data
    type(axisAlignedBoundingBox), pointer :: boundingBoxPtr
    logical(defBool)                      :: stopSearch
    
    boundingBoxPtr => self % getBoundingBoxPtr()
    searchLoop: do
      if (.not. boundingBoxPtr % contains(data % r, data % u)) return
      call self % acceleration % findHostElement(self % elements, data, stopSearch)
      if (stopSearch) return

    end do searchLoop

  end subroutine findHostElement

  !!
  !!
  !!
  function getCastUnstructuredMeshPtr(source) result(ptr)
    class(mesh), intent(in)          :: source
    class(unstructuredMesh), pointer :: ptr
    character(*), parameter          :: here = 'getCastUnstructuredMeshPtr (unstructuredMesh_inter.f90)'

    ! Downcast to correct type.
    select type(temp => source)
      class is(unstructuredMesh)
        ptr => temp

      class default
        call fatalError(here, 'Mesh is not of class unstructuredMesh.')

    end select

  end function getCastUnstructuredMeshPtr

  !!
  !!
  !!
  elemental function getEdgesNumber(self) result(nEdges)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt)                   :: nEdges

    nEdges = self % nEdges

  end function getEdgesNumber

  !!
  !!
  !!
  function getElementBox(self, idx) result(box)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt), intent(in)       :: idx
    type(elementBox)                    :: box

    box = self % elements % getElementBox(idx)

  end function getElementBox

  !!
  !!
  !!
  function getElementIsActive(self, idx) result(isActive)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt), intent(in)       :: idx
    logical(defBool)                    :: isActive
    type(elementBox)                    :: box
    
    box = self % elements % getElementBox(idx)
    isActive = box % ptr % getIsActive()

  end function getElementIsActive

  !!
  !!
  !!
  function getElementsNumber(self, activeOnly) result(nElements)
    class(unstructuredMesh), intent(in)    :: self
    logical(defBool), intent(in), optional :: activeOnly
    integer(shortInt)                      :: nElements

    nElements = self % elements % getObjectsNumber(activeOnly)

  end function getElementsNumber

  !!
  !!
  !!
  function getElementVolume(self, idx) result(volume)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt), intent(in)       :: idx
    real(defReal)                       :: volume
    type(elementBox)                    :: box

    box = self % elements % getElementBox(idx)
    volume = box % ptr % getVolume()

  end function getElementVolume

  !!
  !!
  !!
  function getFaceBoundaryConditions(self, idx) result(boundaryConditions)
    class(unstructuredMesh), intent(in)      :: self
    integer(shortInt), intent(in)            :: idx
    integer(shortInt), dimension(N_BC_TYPES) :: boundaryConditions
    type(faceBox)                            :: box

    box = self % faces % getFaceBox(idx)
    boundaryConditions = box % ptr % getBoundaryConditions()

  end function getFaceBoundaryConditions

  !!
  !!
  !!
  function getFaceIsBoundary(self, idx) result(isBoundary)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt), intent(in)       :: idx
    logical(defBool)                    :: isBoundary
    type(faceBox)                       :: box

    box = self % faces % getFaceBox(idx)
    isBoundary = box % ptr % getIsBoundary()

  end function getFaceIsBoundary

  !!
  !!
  !!
  function getFacesNumber(self, activeOnly) result(nFaces)
    class(unstructuredMesh), intent(in)    :: self
    logical(defBool), intent(in), optional :: activeOnly
    integer(shortInt)                      :: nFaces

    nFaces = self % faces % getObjectsNumber(activeOnly)

  end function getFacesNumber

  !!
  !!
  !!
  function getInternalFacesNumber(self, activeOnly) result(nInternalFaces)
    class(unstructuredMesh), intent(in)    :: self
    logical(defBool), intent(in), optional :: activeOnly
    integer(shortInt)                      :: i, nInternalFaces
    logical(defBool)                       :: filterActive
    type(faceBox)                          :: box

    filterActive = .false.
    if (present(activeOnly)) filterActive = activeOnly

    nInternalFaces = 0
    do i = 1, self % faces % getObjectsNumber()
      box = self % faces % getFaceBox(i)
      if (filterActive .and. .not. box % ptr % getIsActive()) cycle
      if (.not. box % ptr % getIsBoundary()) nInternalFaces = nInternalFaces + 1

    end do

  end function getInternalFacesNumber

  !!
  !!
  !!
  function getParentElementsNumber(self) result(nParentElements)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt)                   :: i, nParentElements
    type(elementBox)                    :: box

    nParentElements = 0
    do i = 1, self % elements % getObjectsNumber()
      box = self % elements % getElementBox(i)
      if (box % ptr % getParentIdx() == 0) nParentElements = nParentElements + 1

    end do

  end function getParentElementsNumber

  !!
  !!
  !!
  function getUniqueIdOffset(self) result(uniqueIdOffset)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt)                   :: uniqueIdOffset

    uniqueIdOffset = self % getElementsNumber()

  end function getUniqueIdOffset

  !!
  !!
  !!
  elemental function getVerticesNumber(self) result(nVertices)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt)                   :: nVertices

    nVertices = self % nVertices

  end function getVerticesNumber

  !!
  !!
  !!
  function importBoundaryConditionsFromFile(self, path) result(boundaryConditionInfos)
    class(unstructuredMesh), intent(in)                        :: self
    character(pathLen), intent(in)                             :: path
    type(meshBoundaryConditionInfo), dimension(:), allocatable :: boundaryConditionInfos
    character(*), parameter :: here = 'importBoundaryConditionsFromFile (unstructuredMesh_inter.f90)'

    ! Call fatalError for now.
    call fatalError(here, 'STOP.')

  end function importBoundaryConditionsFromFile

  !!
  !!
  !!
  function importLocalIdsFromFile(self, path) result(localIdInfos)
    class(unstructuredMesh), intent(in)              :: self
    character(pathLen), intent(in)                   :: path
    type(meshLocalIdInfo), dimension(:), allocatable :: localIdInfos
    character(*), parameter                          :: here = 'importLocalIdsFromFile (unstructuredMesh_inter.f90)'

    ! Call fatalError for now.
    call fatalError(here, 'STOP.')

  end function importLocalIdsFromFile

  !!
  !!
  !!
  subroutine init(self, folderPath, dict, materialsMap)
    class(unstructuredMesh), intent(inout) :: self
    character(*), intent(in)               :: folderPath
    class(dictionary), intent(in)          :: dict
    type(charMap), intent(in)              :: materialsMap
    type(faceBox)                          :: face

    ! Set up base components.
    call self % setupBase(dict)
    
    ! Import mesh from files.
    call self % importMesh(folderPath)

    ! Assign localIds.
    call self % assignLocalIds(dict, materialsMap)

    ! Assign boundary conditions.
    if (dict % isPresent('boundaryConditions')) call self % assignBoundaryConditions(dict % getDictPtr('boundaryConditions'))

    ! Initialise triangulation method from dictionary then triangulate mesh.
    call newTriangulationPtr(dict, self % triangulation)
    call self % triangulation % triangulate(self % edges, self % elements, self % faces, self % vertices)

    ! Shrink shelves to their correct size after triangulation.
    call self % edges % shrink()
    call self % elements % shrink()
    call self % faces % shrink()
    call self % vertices % shrink()

    ! Update number of edges, elements, faces, internal faces, and vertices.
    self % nEdges = self % edges % getObjectsNumber()
    self % nVertices = self % vertices % getObjectsNumber()

    ! Initialise acceleration method from dictionary.
    call newAccelerationStructurePtr(dict, self % edges, self % elements, self % faces, self % vertices, self % acceleration)

  end subroutine init

  !!
  !!
  !!
  subroutine initEdgeShelf(self, edgeInfos)
    class(unstructuredMesh), intent(inout)                                :: self
    type(basicEdgeInfo), dimension(:), intent(in)                         :: edgeInfos
    type(buildExtentTopologicalObjectPayload), dimension(size(edgeInfos)) :: payloads
    integer(shortInt)                                                     :: i, nEdges

    nEdges = size(edgeInfos)
    self % nEdges = nEdges
    do i = 1, nEdges
      payloads(i) % idx = edgeInfos(i) % idx
      payloads(i) % vertices = self % vertices % getVertexBox(edgeInfos(i) % vertexIdxs)

    end do
    call self % edges % init(payloads)

  end subroutine initEdgeShelf

  !!
  !!
  !!
  subroutine initElementShelf(self, elementInfos)
    class(unstructuredMesh), intent(inout)                   :: self
    type(basicElementInfo), dimension(:), intent(inout)      :: elementInfos
    type(buildElementPayload), dimension(size(elementInfos)) :: payloads
    integer(shortInt)                                        :: i, j, nElements, nFaces
    type(faceBox)                                            :: face

    nElements = size(elementInfos)
    do i = 1, nElements
      payloads(i) % idx = elementInfos(i) % idx
      payloads(i) % parentIdx = elementInfos(i) % parentIdx
      nFaces = size(elementInfos(i) % faceIdxs)
      allocate(payloads(i) % orientatedFaces(nFaces))
      do j = 1, nFaces
        face = self % faces % getFaceBox(abs(elementInfos(i) % faceIdxs(j)))
        payloads(i) % orientatedFaces(j) % face = face
        if (0 < elementInfos(i) % faceIdxs(j)) then
          payloads(i) % orientatedFaces(j) % isOwner = .true.
          payloads(i) % orientatedFaces(j) % outwardNormal = face % ptr % getNormal()
          payloads(i) % orientatedFaces(j) % ratintOutwardNormal = face % ptr % getRatintNormal()

        else
          payloads(i) % orientatedFaces(j) % outwardNormal = -face % ptr % getNormal()
          payloads(i) % orientatedFaces(j) % ratintOutwardNormal = (-1_8)*(face % ptr % getRatintNormal())

        end if

      end do

      ! Check if we need to construct edges.
      if (.not. allocated(elementInfos(i) % edgeIdxs)) call createEdgeIdxs(elementInfos(i))
      payloads(i) % edges = self % edges % getEdgeBox(elementInfos(i) % edgeIdxs)

      ! Check if we need to construct vertices.
      if (.not. allocated(elementInfos(i) % vertexIdxs)) call createVertexIdxs(elementInfos(i))
      payloads(i) % vertices = self % vertices % getVertexBox(elementInfos(i) % vertexIdxs)

    end do
    call self % elements % init(payloads)
  
  contains
    !!
    !!
    !!
    subroutine createEdgeIdxs(info)
      type(basicElementInfo), intent(inout)                 :: info
      logical(defBool), dimension(self % edges % getSize()) :: isPresent
      integer(shortInt)                                     :: currentSize, idx, k, l, nEdges
      type(faceBox)                                         :: fBox
      type(edgeBox), dimension(:), allocatable              :: faceEdges
      integer(shortInt), dimension(:), allocatable          :: tempIdxs

      ! Initialise isPresent = .false. and allocate info % edgeIdxs to an appropriate initial size.
      isPresent = .false.
      allocate(info % edgeIdxs(6))
      nEdges = 0
      do k = 1, size(info % faceIdxs)
        fBox = self % faces % getFaceBox(abs(info % faceIdxs(k)))
        faceEdges = fBox % ptr % getEdges()
        do l = 1, size(faceEdges)
          idx = faceEdges(l) % ptr % getIdx()
          if (.not. isPresent(idx)) then
            nEdges = nEdges + 1
            currentSize = size(info % edgeIdxs)
            if (currentSize < nEdges) then
              allocate(tempIdxs(2 * currentSize))
              tempIdxs(1:currentSize) = info % edgeIdxs
              call move_alloc(tempIdxs, info % edgeIdxs)

            end if
            info % edgeIdxs(nEdges) = idx
            isPresent(idx) = .true.

          end if

        end do

      end do

      ! Resize info % edgeIdxs to correct size if necessary.
      if (nEdges < size(info % edgeIdxs)) then
        allocate(tempIdxs(nEdges))
        tempIdxs = info % edgeIdxs(1:nEdges)
        call move_alloc(tempIdxs, info % edgeIdxs)

      end if

    end subroutine createEdgeIdxs

    !!
    !!
    !!
    subroutine createVertexIdxs(info)
      type(basicElementInfo), intent(inout)                    :: info
      logical(defBool), dimension(self % vertices % getSize()) :: isPresent
      integer(shortInt)                                        :: currentSize, idx, k, l, nVertices
      type(faceBox)                                            :: fBox
      type(vertexBox), dimension(:), allocatable               :: faceVertices
      integer(shortInt), dimension(:), allocatable             :: tempIdxs

      ! Initialise isPresent = .false. and allocate info % vertexIdxs to an appropriate initial size.
      isPresent = .false.
      allocate(info % vertexIdxs(4))
      nVertices = 0
      do k = 1, size(info % faceIdxs)
        fBox = self % faces % getFaceBox(abs(info % faceIdxs(k)))
        faceVertices = fBox % ptr % getVertices()
        do l = 1, size(faceVertices)
          idx = faceVertices(l) % ptr % getIdx()
          if (.not. isPresent(idx)) then
            nVertices = nVertices + 1
            currentSize = size(info % vertexIdxs)
            if (currentSize < nVertices) then
              allocate(tempIdxs(2 * currentSize))
              tempIdxs(1:currentSize) = info % vertexIdxs
              call move_alloc(tempIdxs, info % vertexIdxs)

            end if
            info % vertexIdxs(nVertices) = idx
            isPresent(idx) = .true.

          end if

        end do

      end do

      ! Resize info % vertexIdxs to correct size if necessary.
      if (nVertices < size(info % vertexIdxs)) then
        allocate(tempIdxs(nVertices))
        tempIdxs = info % vertexIdxs(1:nVertices)
        call move_alloc(tempIdxs, info % vertexIdxs)

      end if

    end subroutine createVertexIdxs

  end subroutine initElementShelf

  !!
  !!
  !!
  subroutine initFaceShelf(self, faceInfos)
    class(unstructuredMesh), intent(inout)             :: self
    type(basicFaceInfo), dimension(:), intent(in)      :: faceInfos
    type(buildFacePayload), dimension(size(faceInfos)) :: payloads
    integer(shortInt)                                  :: i, nFaces

    nFaces = size(faceInfos)
    do i = 1, nFaces
      payloads(i) % idx = faceInfos(i) % idx
      payloads(i) % parentIdx = faceInfos(i) % parentIdx
      payloads(i) % isBoundary = faceInfos(i) % isBoundary
      payloads(i) % vertices = self % vertices % getVertexBox(faceInfos(i) % vertexIdxs)
      payloads(i) % edges = self % edges % getEdgeBox(faceInfos(i) % edgeIdxs)

    end do
    call self % faces % init(payloads)

  end subroutine initFaceShelf

  !!
  !!
  !!
  subroutine initVertexShelf(self, vertexInfos)
    class(unstructuredMesh), intent(inout)                 :: self
    type(basicVertexInfo), dimension(:), intent(in)        :: vertexInfos
    type(buildVertexPayload), dimension(size(vertexInfos)) :: payloads
    integer(shortInt)                                      :: i, nVertices
    real(defReal), dimension(3)                            :: coords
    real(defReal), dimension(:, :), allocatable            :: allCoords

    nVertices = size(vertexInfos)
    self % nVertices = nVertices
    allocate(allCoords(3, nVertices))

    do i = 1, nVertices
      payloads(i) % idx = vertexInfos(i) % idx
      coords = vertexInfos(i) % coordinates
      payloads(i) % coordinates = coords
      payloads(i) %ratintCoordinates = vertexInfos(i) %ratintCoordinates
      allCoords(:, i) = coords

    end do
    call self % vertices % init(payloads)
    call self % initBoundingBox(allCoords)

  end subroutine

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an unitialised state.
  !!
  subroutine kill(self)
    class(unstructuredMesh), intent(inout) :: self

    ! Superclass.
    call kill_super(self)
    
    ! Local.
    call self % acceleration % kill()
    deallocate(self % acceleration)
    self % nVertices = 0
    self % nInternalFaces = 0
    self % nEdges = 0
    call self % elements % kill()
    call self % faces % kill()
    call self % edges % kill()
    call self % vertices % kill()
    deallocate(self % triangulation)

  end subroutine kill

  !! Subroutine 'printComposition'
  !!
  !! Basic description:
  !!   Prints the initial polyhedral composition of the mesh.
  !!
  !! Arguments:
  !!   nTetrahedra [out] -> Number of tetrahedra in the mesh.
  !!
  subroutine printComposition(self, nTetrahedra)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt), intent(out)      :: nTetrahedra
    integer(shortInt)                   :: nFaces, nPentahedra, nHexahedra, nOthers, i
    type(elementBox)                    :: element
    
    ! Initialise the numbers of various polyhedra to zero.
    nTetrahedra = 0
    nPentahedra = 0
    nHexahedra = 0
    nOthers = 0
    
    ! Loop over all elements in the mesh.
    do i = 1, self % elements % getObjectsNumber()
      ! Retrieve the number of faces in the current element and increment specific polyhedra
      ! accordingly.
      element = self % elements % getElementBox(i)
      nFaces = size(element % ptr % getOrientatedFaces())
      select case (nFaces)
        case (4)
          nTetrahedra = nTetrahedra + 1
        case (5)
          nPentahedra = nPentahedra + 1
        case (6)
          nHexahedra = nHexahedra + 1
        case default
          nOthers = nOthers + 1

      end select

    end do
    
    ! Print to screen.
    print *, 'Displaying unstructured mesh composition:'
    print *, '  Number of tetrahedra     : '//numToChar(nTetrahedra)//'.'
    print *, '  Number of pentahedra     : '//numToChar(nPentahedra)//'.'
    print *, '  Number of hexahedra      : '//numToChar(nHexahedra)//'.'
    print *, '  Number of other polyhedra: '//numToChar(nOthers)//'.'

  end subroutine printComposition

  !!
  !!
  !!
  subroutine sampleInitialPosition(self, elementIdx, rand, localId, r)
    class(unstructuredMesh), intent(in)      :: self
    integer(shortInt), intent(in)            :: elementIdx
    type(RNG), intent(inout)                :: rand
    integer(shortInt), intent(out)           :: localId
    real(defReal), dimension(3), intent(out) :: r
    type(elementBox)                         :: box

    box = self % elements % getElementBox(elementIdx)
    call box % ptr % sampleInitialPosition(rand, localId, r)

  end subroutine sampleInitialPosition

  !!
  !!
  !!
  elemental subroutine setEdgesNumber(self, nEdges)
    class(unstructuredMesh), intent(inout) :: self
    integer(shortInt), intent(in)          :: nEdges

    self % nEdges = nEdges

  end subroutine setEdgesNumber

  !!
  !!
  !!
  elemental subroutine setInternalFacesNumber(self, nInternalFaces)
    class(unstructuredMesh), intent(inout) :: self
    integer(shortInt), intent(in)          :: nInternalFaces

    self % nInternalFaces = nInternalFaces

  end subroutine setInternalFacesNumber

  !!
  !!
  !!
  elemental subroutine setVerticesNumber(self, nVertices)
    class(unstructuredMesh), intent(inout) :: self
    integer(shortInt), intent(in)          :: nVertices

    self % nVertices = nVertices

  end subroutine setVerticesNumber

end module unstructuredMesh_inter