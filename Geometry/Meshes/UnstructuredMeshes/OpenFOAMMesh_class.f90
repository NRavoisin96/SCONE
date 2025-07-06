module OpenFOAMMesh_class

  use axisAlignedBoundingBox_class, only : axisAlignedBoundingBox
  use coord_class,                  only : coord
  use genericProcedures,            only : append, fatalError, numToChar, openToRead, quickSort
  use iso_fortran_env,              only : int64
  use longIntMap_class,             only : longIntMap
  use numPrecision
  use publicObjects,                only : basicEdgeInfo, basicElementInfo, basicFaceInfo, meshLocalIdInfo
  use universalVariables,           only : INF, NOT_PRESENT
  use unstructuredMesh_inter,       only : unstructuredMesh, &
                                           distanceToBoundaryFace_super => distanceToBoundaryFace, &
                                           distanceToNextFace_super => distanceToNextFace, &
                                           findHostElement_super => findHostElement, &
                                           kill_super => kill

  implicit none
  private

  type, public, extends(unstructuredMesh) :: OpenFOAMMesh
    private
  contains
    ! Superclass procedures.
    procedure :: distanceToBoundaryFace
    procedure :: distanceToNextFace
    procedure :: findHostElement
    ! Local procedures.
    procedure :: checkFiles
    procedure :: getMeshInfo
    procedure :: importElements
    procedure :: importElementZones
    procedure :: importFacesAndEdges
    procedure :: importMesh
    procedure :: importVertices
    procedure :: kill
  end type OpenFOAMMesh

contains
  !! Subroutine 'checkFiles'
  !!
  !! Basic description:
  !!   The subroutine 'checkFiles' checks the existence of vital mesh data files. These files
  !!   are the 'points', 'faces' and 'neighbour' files. If any of them are missing, the subroutine
  !!   calls the procedure 'fatalError'.
  !!
  !! Arguments:
  !!   folderPath [in]     -> Path of the folder containing the files for the mesh.
  !!   nInternalFaces [in] -> Number of internal faces in the mesh. If this number is
  !!                          zero then the 'neighbour' file is not checked.
  !!
  !! Errors:
  !!   - fatalError if the 'points' file is missing.
  !!   - fatalError if the 'faces' file is missing.
  !!   - fatalError if nInternalFaces > 0 and the 'neighbour' file is missing.
  !!
  !! Notes:
  !!   The existence of the 'owner' file is checked in 'getMeshInfo'.
  !!
  subroutine checkFiles(self, folderPath, nInternalFaces, cellZonesFile)
    class(OpenFOAMMesh), intent(inout) :: self
    character(*), intent(in)           :: folderPath
    integer(shortInt), intent(in)      :: nInternalFaces
    logical(defBool), intent(out)      :: cellZonesFile
    logical(defBool)                   :: pointsFile, facesFile, neighbourFile
    character(*), parameter            :: Here = 'checkFiles (OpenFOAMMesh_class.f90)'

    ! Check the existence of the 'points' and 'faces' files and report errors if they are not found.
    inquire(file = folderPath//'points', exist = pointsFile)
    if (.not. pointsFile) call fatalError(Here, "Missing 'points' file for OpenFOAM mesh with Id: "//numToChar(self % getId())//'.')
    
    inquire(file = folderPath//'faces', exist = facesFile)
    if (.not. facesFile) call fatalError(Here, "Missing 'faces' file for OpenFOAM mesh with Id: "//numToChar(self % getId())//'.')
    
    ! If nInternal = 0, return early.
    if (nInternalFaces == 0) return

    ! If reached here check that the 'neighbour' file exists and report error if not.
    inquire(file = folderPath//'neighbour', exist = neighbourFile)
    if (.not. neighbourFile) &
    call fatalError(Here, "Missing 'neighbour' file for OpenFOAM mesh with Id: "//numToChar(self % getId())//'.')

    ! Check that the 'cellZones' file exists.
    inquire(file = folderPath//'cellZones', exist = cellZonesFile)

  end subroutine checkFiles

  !! Subroutine 'distanceToBoundaryFace'
  !!
  !! Basic description:
  !!   Returns the distance to the mesh boundary face intersected by a particle's path. Also returns the index
  !!   of the parent element containing the intersected boundary face.
  !!
  !! See unstructuredMesh_inter for details.
  !!
  subroutine distanceToBoundaryFace(self, d, coords)
    class(OpenFOAMMesh), intent(in) :: self
    real(defReal), intent(out)      :: d
    type(coord), intent(inout)      :: coords

    call distanceToBoundaryFace_super(self, d, coords)

  end subroutine distanceToBoundaryFace

  !! Subroutine 'distanceToNextFace'
  !!
  !! Basic description:
  !!   Returns the distance to the next face intersected by the particle's path.
  !!
  !! See unstructuredMesh_inter for details.
  !!
  subroutine distanceToNextFace(self, d, coords)
    class(OpenFOAMMesh), intent(in) :: self
    real(defReal), intent(out)      :: d
    type(coord), intent(inout)      :: coords

    call distanceToNextFace_super(self, d, coords)

  end subroutine distanceToNextFace

  !! Subroutine 'findElementAndParentIdxs'
  !!
  !! Basic description:
  !!   Returns the index of the mesh element occupied by a particle. Also returns the index of the parent mesh
  !!   element containing the occupied element.
  !!
  !! See unstructuredMesh_inter for details.
  !!
  subroutine findHostElement(self, coords)
    class(OpenFOAMMesh), intent(in)        :: self
    type(coord), intent(inout)             :: coords

    call findHostElement_super(self, coords)

  end subroutine findHostElement

  !! Subroutine 'getMeshInfo'
  !!
  !! Basic description:
  !!   Retrieves some preliminary information about the mesh (number of vertices, faces, etc.) to be 
  !!   imported using the 'owner' data file.
  !!
  !! Detailed description:
  !!   Opens the 'owner' file and reads it until a line with the keyword 'note' is encountered. From
  !!   this line of text, the number of vertices (nVertices), the number of faces (nFaces), the 
  !!   number of elements (nElements) and the number of internal faces (nInternalFaces) are read.
  !!
  !! Arguments:
  !!   folderPath [in] -> Path of the folder containing the files of the mesh.
  !!
  subroutine getMeshInfo(self, folderPath, nVertices, nFaces, nInternalFaces, nElements, hasCellZones)
    class(OpenFOAMMesh), intent(inout) :: self
    character(*), intent(in)           :: folderPath
    integer(shortInt), intent(out)     :: nVertices, nFaces, nInternalFaces, nElements
    logical(defBool), intent(out)      :: hasCellZones
    logical(defBool)                   :: ownerFile
    integer(shortInt)                  :: position
    integer(shortInt), parameter       :: unit = 10
    character(256)                     :: buffer
    character(:), allocatable          :: ownerPath
    character(*), parameter            :: Here = 'getMeshInfo (OpenFOAMMesh_class.f90)'

    ! Initialise ownerPath and perform a first check for the existence of the 'owner' file. Call fatalError if not found.
    ownerPath = folderPath//'owner'
    inquire(file = ownerPath, exist = ownerFile)
    if (.not. ownerFile) call fatalError(Here, "Missing 'owner' file for OpenFOAM mesh with Id: "//numToChar(self % getId())//'.')
    
    ! Open the 'owner' file and read it. Skip lines until the line with the keyword 'note' is encountered.
    call openToRead(unit, ownerPath)
    read(unit, "(a)") buffer
    do while (index(buffer, 'note') == 0)
      read(unit, "(a)") buffer

    end do
    
    ! Find 'nPoints' and read it.
    position = index(buffer, 'nPoints:')
    if (position == 0) call fatalError(here, "Could not find number of vertices in 'owner' file.")
    read(buffer(position + 8:), *) nVertices
    
    ! Find 'nCells' and read it.
    position = index(buffer, 'nCells:')
    if (position == 0) call fatalError(here, "Could not find number of elements in 'owner' file.")
    read(buffer(position + 7:), *) nElements

    ! Find 'nFaces' and read it.
    position = index(buffer, 'nFaces:')
    if (position == 0) call fatalError(here, "Could not find number of faces in 'owner' file.")
    read(buffer(position + 7:), *) nFaces

    ! Find 'nInternalFaces' and read it.
    position = index(buffer, 'nInternalFaces:')
    if (position == 0) call fatalError(here, "Could not find number of internal faces in 'owner' file.")
    read(buffer(position + 15:index(buffer, ';') - 2), *) nInternalFaces
    
    ! Close 'owner' file and check existence of remaining mesh files.
    close(unit)
    call self % checkFiles(folderPath, nInternalFaces, hasCellZones)

  end subroutine getMeshInfo

  !! Subroutine 'initElementShelf'
  !!
  !! Basic description:
  !!   Initialises the shelf from the 'owner' and 'neighbour' files.
  !!
  !! Arguments:
  !!   folderPath [in]     -> Path of the folder containing the OpenFOAM mesh files.
  !!   nFaces [in]         -> Number of faces in the mesh.
  !!   nInternalFaces [in] -> Number of internal faces in the mesh.
  !!
  subroutine importElements(self, folderPath, nVertices, nFaces, nInternalFaces, nElements)
    class(OpenFOAMMesh), intent(inout)           :: self
    character(*), intent(in)                     :: folderPath
    integer(shortInt), intent(in)                :: nVertices, nFaces, nInternalFaces, nElements
    integer(shortInt)                            :: i, elementIdx
    integer(shortInt), parameter                 :: unit = 10
    integer(shortInt), dimension(:), allocatable :: elementIdxs
    character(256)                               :: buffer
    type(basicElementInfo), dimension(nElements) :: elementInfos

    ! Loop over all elements first and assign their indices.
    do i = 1, nElements
      elementInfos(i) % idx = i

    end do

    ! If there is only one element in the mesh simply add all the faces and vertices to this element.
    if (nElements == 1) then
      allocate(elementInfos(1) % faceIdxs(nFaces))
      do i = 1, nFaces
        elementInfos(1) % faceIdxs(i) = i

      end do

      allocate(elementInfos(1) % vertexIdxs(nVertices))
      do i = 1, nVertices
        elementInfos(1) % vertexIdxs(i) = i

      end do

    else
      ! Open the 'owner' file and read it until a line containing the symbol '(' is encountered.
      call openToRead(unit, folderPath//'owner')
      read(unit, "(a)") buffer
      do while (index(buffer(1:len_trim(buffer)), "(") == 0)
        read(unit, "(a)") buffer

      end do

      ! Skip one more line and loop over all faces.
      read(unit, "(a)") buffer
      do i = 1, nFaces
        ! Read the current element index and add the current face to this element. Note
        ! that we add one since Fortran starts indexing at one and not zero.
        read(buffer, *) elementIdx

        elementIdx = elementIdx + 1
        call append(elementInfos(elementIdx) % faceIdxs, i)

        ! Move onto the next line.
        read(unit, "(a)") buffer

      end do

      ! Close the 'owner' file.
      close(unit)

      ! If nInternalFaces = 0 we can return early here. Else we need to repeat the above procedure
      ! for the 'neighbour file.
      if (nInternalFaces > 0) then
        ! Open the 'neighbour' file and read it until a line containing the symbol '(' is encountered.
        call openToRead(unit, folderPath//'neighbour')
        read(unit, "(a)") buffer
        do while (index(buffer(1:len_trim(buffer)), "(") == 0)
          read(unit, "(a)") buffer

        end do

        ! Check if the current line contains the symbol ')'. If it does, then all element indices
        ! are written on a single line.
        if (index(buffer(1:len_trim(buffer)), ")") > 0) then
          ! Allocate the number of entries in the 'elementIndices' array to the number of internal
          ! faces and copy element indices into this array.
          allocate(elementIdxs(nInternalFaces))
          read(buffer(index(buffer(1:len_trim(buffer)), "(") + 1:&
          index(buffer(1:len_trim(buffer)), ")") - 1), *) elementIdxs

        else
          ! Skip one more line.
          read(unit, "(a)") buffer

        end if

        ! Loop over all internal faces.
        do i = 1, nInternalFaces
          if (allocated(elementIdxs)) then
            elementIdx = elementIdxs(i)

          else
            ! Read the current element index and move onto the next line.
            read(buffer, *) elementIdx
            read(unit, "(a)") buffer

          end if
          ! Update connectivity information.
          elementIdx = elementIdx + 1
          call append(elementInfos(elementIdx) % faceIdxs, -i)

        end do

        ! Close the 'neighbour' file.
        close(unit)

      end if

    end if

    ! Initialise elementShelf.
    call self % initElementShelf(elementInfos)

  end subroutine importElements

  !! Subroutine 'initCellZoneShelf'
  !!
  !! Basic description:
  !!   Initialises the shelf from the cellZones file.
  !!
  !! Notes:
  !!   Due to the structure of the cellZones file the subroutine first scans the file a first time
  !!   to obtain the name of each cell zone and then rewinds it to read the indices of the elements
  !!   contained in each zone.
  !!
  !! Arguments:
  !!   folderPath [in] -> Path of the folder containing the OpenFOAM mesh files.
  !!
  subroutine importElementZones(self, folderPath, hasCellZones, nElements)
    class(OpenFOAMMesh), intent(inout)               :: self
    character(*), intent(in)                         :: folderPath
    logical(defBool), intent(in)                     :: hasCellZones
    integer(shortInt), intent(in)                    :: nElements
    integer(shortInt)                                :: i, j, nElementsInZone, nElementZones
    integer(shortInt), parameter                     :: unit = 10
    character(256)                                   :: buffer
    logical(defBool)                                 :: singleLine
    type(meshLocalIdInfo), dimension(:), allocatable :: localIdInfos

    ! Check if there are element zones. If not, simply allocate one element zone and set all elements to it.
    if (.not. hasCellZones) then
      allocate(localIdInfos(1))
      localIdInfos(1) % localId = 2
      localIdInfos(1) % elementIdxs = [(i, i = 1, nElements)]

    else
      ! Open the 'cellZones' data file and read it. Skip lines until a blank line is encountered.
      call openToRead(unit, folderPath//'cellZones')
      read(unit, "(a)") buffer
      do while (len_trim(buffer) > 0)
        read(unit, "(a)") buffer
      end do
      
      ! Read the current line and copy the number of cell zones into the variable 'nCellZones'.
      read(unit, "(a)") buffer
      read(buffer, *) nElementZones
      
      ! Allocate memory and loop through all cell zones.
      allocate(localIdInfos(nElementZones))
      do i = 1, nElementZones
        localIdInfos(i) % localId = i + 1
        ! Initialise singleLine = .false. and skip lines until a '{' is encountered.
        singleLine = .false.
        do while (index(buffer, "{") == 0)
          read(unit, "(a)") buffer

        end do

        ! Go back to the previous line and read the name of the current cell zone.
        backspace(unit)
        read(unit, "(a)") buffer
        
        ! Skip lines until the word 'cellLabels' is encoutered.
        do while (index(buffer, "cellLabels") == 0)
          read(unit, "(a)") buffer

        end do

        ! If the line contains a ')' then all the elements contained in the current cell zones are
        ! written on the same line.
        if (index(buffer, ")") > 0) singleLine = .true.

        ! Read the number of elements in the current cell zone depending on whether they are all
        ! written a single line or not.
        if (singleLine) then
          read(buffer(index(buffer, ">") + 2:index(buffer, "(") - 1), *) nElementsInZone

        else
          ! Skip one line and retrieve the number of elements in the cell zone..
          read(unit, "(a)") buffer
          read(buffer, *) nElementsInZone

        end if

        ! Allocate memory.
        allocate(localIdInfos(i) % elementIdxs(nElementsInZone))

        ! Now read the indices of the elements in the current cell zone. Again this depends on
        ! whether they are all written on a single line or not.
        if (singleLine) then
          ! Retrieve the element indices.
          read(buffer(index(buffer, "(") + 1:index(buffer, ")") - 1), *) localIdInfos(i) % elementIdxs

        else
          ! Skip two lines and retrieve the index of each element in the cell zone line by line..
          do j = 1, 2
            read(unit, "(a)") buffer

          end do
          do j = 1, nElementsInZone
            read(buffer, *) localIdInfos(i) % elementIdxs(j)
            read(unit, "(a)") buffer

          end do

        end if
        localIdInfos(i) % elementIdxs = localIdInfos(i) % elementIdxs + 1

      end do
      
      ! Close the 'cellZones' file.
      close(unit)

    end if

    call self % assignLocalIds(localIdInfos)

  end subroutine importElementZones

  !! Subroutine 'initFaceShelf'
  !!
  !! Basic description:
  !!   Initialises the shelf from the 'faces' file.
  !!
  !! Detailed description:
  !!   Opens the 'faces' file and reads it until a line with the symbol ')' is encountered, which 
  !!   marks the beginning of the faces' data listing. In the 'faces' file, each line corresponds to
  !!   a single face. The subroutine then retrieves the number of vertices in each face and their
  !!   indices sequentially line-by-line.
  !!
  !! Arguments:
  !!   folderPath [in]     -> Path of the folder containing the mesh files.
  !!   nInternalFaces [in] -> Number of internal faces in the mesh.
  !!
  subroutine importFacesAndEdges(self, folderPath, nFaces, nInternalFaces)
    class(OpenFOAMMesh), intent(inout)             :: self
    character(*), intent(in)                       :: folderPath
    integer(shortInt), intent(in)                  :: nFaces, nInternalFaces
    integer(shortInt)                              :: currentSize, edgeIdx, i, idx, ios, j, lastEdgeIdx, &
                                                      maxVertexIdx, minVertexIdx, nVertices
    integer(int64)                                 :: key
    integer(shortInt), parameter                   :: unit = 10
    character(256)                                 :: buffer
    type(basicEdgeInfo), dimension(:), allocatable :: edgeInfos, tempEdgesInfo
    type(basicFaceInfo), dimension(nFaces)         :: faceInfos
    type(longIntMap)                               :: edgeIdxsMap
    character(*), parameter                        :: here = 'importFacesAndEdges (OpenFOAMMesh_class.f90)'

    ! Open the 'faces' data file and read it until a line containing the symbol ')' is encountered.
    call openToRead(unit, folderPath//'faces')
    do
      read(unit, "(a)", iostat = ios) buffer
      if (ios /= 0) call fatalError(here, 'Could not find beginning of faces data block.')
      if (index(trim(buffer), ')') > 0) exit

    end do
    
    ! Loop through all faces in the file.
    do i = 1, nFaces
      ! Check if the vertex list starts on the same line as the count.
      idx = index(buffer, '(')
      if (idx > 0) then
        ! Single-line format, e.g. 4(1 2 3 4)
        read(buffer(1:idx - 1), *) nVertices
        allocate(faceInfos(i) % vertexIdxs(nVertices))
        read(buffer(idx + 1:len_trim(buffer) - 1), *) faceInfos(i) % vertexIdxs

      else
        ! Multi-line format.
        read(buffer, *) nVertices
        allocate(faceInfos(i) % vertexIdxs(nVertices))
        ! Skip the '(' line.
        read(unit, "(a)")
        do j = 1, nVertices
          read(unit, *) faceInfos(i) % vertexIdxs(j)

        end do
        ! Skip the ')' line.
        read(unit, "(a)")

      end if
      
      ! Add one to the vertices indices since Fortran starts indexing at one rather than zero. Then
      ! loop through all vertices in the face.
      faceInfos(i) % idx = i
      faceInfos(i) % isBoundary = nInternalFaces < i
      faceInfos(i) % vertexIdxs = faceInfos(i) % vertexIdxs + 1
      read(unit, "(a)", iostat = ios) buffer

      ! Skip blank lines (can happen for faces with a large number of vertices).
      do while(len_trim(buffer) == 0)
        read(unit, "(a)", iostat = ios) buffer
        if (ios /= 0) exit

      end do
      ! Exit main loop if file ends.
      if (ios /= 0) exit

    end do
    
    ! Close the 'faces' file.
    close(unit)

    ! At this point, we have the information for all the faces. Start by building edges.
    ! Allocate some initial size for edgeInfos and loop through all the faces.
    allocate(edgeInfos(2 * nFaces))
    lastEdgeIdx = 0
    do i = 1, nFaces
      ! Loop through all the vertices in the face.
      nVertices = size(faceInfos(i) % vertexIdxs)
      allocate(faceInfos(i) % edgeIdxs(nVertices))
      do j = 1, nVertices
        ! Compute vertices of minimum and maximum indices.
        minVertexIdx = min(faceInfos(i) % vertexIdxs(j), faceInfos(i) % vertexIdxs(merge(1, j + 1, j == nVertices)))
        maxVertexIdx = max(faceInfos(i) % vertexIdxs(j), faceInfos(i) % vertexIdxs(merge(1, j + 1, j == nVertices)))
        
        ! Create key and find if it already is in the edgeIdxsMap.
        key = ishft(int(minVertexIdx, int64), 32) + int(maxVertexIdx, int64)
        edgeIdx = edgeIdxsMap % getOrDefault(key, NOT_PRESENT)

        ! If not found, create a new edge.
        if (edgeIdx == NOT_PRESENT) then
          lastEdgeIdx = lastEdgeIdx + 1
          
          ! Double size of edgeInfos if necessary.
          currentSize = size(edgeInfos)
          if (currentSize < lastEdgeIdx) then
            allocate(tempEdgesInfo(2 * currentSize))
            tempEdgesInfo(1:currentSize) = edgeInfos
            call move_alloc(tempEdgesInfo, edgeInfos)

          end if

          ! Now add the new edge and insert it in the map.
          edgeInfos(lastEdgeIdx) % idx = lastEdgeIdx
          edgeInfos(lastEdgeIdx) % vertexIdxs = [minVertexIdx, maxVertexIdx]
          call edgeIdxsMap % add(key, lastEdgeIdx)

          ! Add newly created edge to faceInfo.
          faceInfos(i) % edgeIdxs(j) = lastEdgeIdx

        else
          ! Simply add the edge found to faceInfo.
          faceInfos(i) % edgeIdxs(j) = edgeIdx

        end if

      end do

    end do

    ! Resize edgeInfos if necessary.
    if (lastEdgeIdx < size(edgeInfos)) then
      allocate(tempEdgesInfo(lastEdgeIdx))
      tempEdgesInfo = edgeInfos(1:lastEdgeIdx)
      call move_alloc(tempEdgesInfo, edgeInfos)

    end if

    ! Initialise edgeShelf and faceShelf.
    call self % initEdgeShelf(edgeInfos)
    call self % initFaceShelf(faceInfos)

  end subroutine importFacesAndEdges

  !! Subroutine 'importMesh'
  !!
  !! Basic description:
  !!   Imports an OpenFOAM mesh from the path of the folder containing the mesh files and sets the 
  !!   mesh Id from the dictionary.
  !!
  !! Detailed description:
  !!   'init' first sets the mesh Id from the supplied dictionary. It then checks the existence and 
  !!   consistency of files located in the appropriate mesh folder and allocates memory to the 
  !!   'vertices', 'faces' and 'elements' structures of the 'mesh' structure. It then imports data 
  !!   from the 'points' file into the 'vertices' structures and builds a kd-tree for the mesh from 
  !!   the various vertices' coordinates. The subroutine then proceeds to import data from the 
  !!   'faces', 'owner' and 'neighbour' files and stores it into the appropriate 'faces' and 
  !!   'elements' structure. 'init' also computes the area, centroid and normal vector for each 
  !!   face, as well as the centroid and volume of each element. From the 'faces', 'owner' and 
  !!   'neighbour' files, mesh connectivity information is also assigned to the various structures. 
  !!   The subroutine then imports data from the 'cellZones' file (if it exists) into the 
  !!   'cellZones' structures.
  !!
  !! Arguments:
  !!   folderPath [in] -> Path of the folder containing the files for the mesh geometry.
  !!   dict [in]       -> Input dictionary.
  !!
  !! Errors:
  !!   - fatalError if the mesh contains concave elements.
  !!
  subroutine importMesh(self, folderPath)
    class(OpenFOAMMesh), intent(inout) :: self
    character(*), intent(in)           :: folderPath
    integer(shortInt)                  :: nElements, nFaces, nInternalFaces, nVertices
    logical(defBool)                   :: hasCellZones
    character(*), parameter            :: Here = 'importMesh (OpenFOAMMesh_class.f90)'
    
    ! Retrieve preliminary information about the mesh.
    call self % getMeshInfo(folderPath, nVertices, nFaces, nInternalFaces, nElements, hasCellZones)
    
    ! Import vertices.
    call self % importVertices(folderPath, nVertices)

    ! Import faces and edges.
    call self % importFacesAndEdges(folderPath, nFaces, nInternalFaces)
    
    ! Import elements.
    call self % importElements(folderPath, nVertices, nFaces, nInternalFaces, nElements)

    ! Import element zones.
    call self % importElementZones(folderPath, hasCellZones, nElements)

  end subroutine importMesh

  !! Subroutine 'initVertexShelf'
  !!
  !! Basic description:
  !!   Initialises the shelf from the 'points' file.
  !!
  !! Arguments:
  !!   folderPath [in] -> Path to the folder containing the mesh files.
  !!
  subroutine importVertices(self, folderPath, nVertices)
    class(OpenFOAMMesh), intent(inout)          :: self
    character(*), intent(in)                    :: folderPath
    integer(shortInt), intent(in)               :: nVertices
    integer(shortInt)                           :: i, ios
    integer(shortInt), parameter                :: unit = 10
    real(defReal), dimension(:, :), allocatable :: coords
    logical(defBool)                            :: singleLine
    character(:), allocatable                   :: dataBuffer
    character(256)                              :: lineBuffer ! Note: here the string is longer than usual to deal
                                                              ! with cases when all vertices are written on a
                                                              ! single line.
    character(*), parameter                     :: here = 'initVertexShelf (OpenFOAMMEsh_class.f90)'

    ! Open the 'points' file.
    call openToRead(unit, folderPath//'points')

    ! Read until the start line is encountered.
    do
      read(unit, '(a)', iostat = ios) lineBuffer
      if (ios /= 0) call fatalError(here, 'Could not find start of data block.')
      if (index(trim(lineBuffer), '(') > 0) exit ! Exit when we find the start

    end do

    ! Allocate memory to the coords array.
    allocate(coords(3, nVertices))

    ! The buffer now holds either the single data line '((-0.5...)...)' or the opening parenthesis line '('.
    ! Check if the line contains a closing parenthesis to determine the format.
    singleLine = index(trim(lineBuffer), ')') > 0

    ! Retrieve coordinates depending on whether they are all on a single line or not.
    if (singleLine) then
      dataBuffer = trim(lineBuffer)
      ! Replace all parentheses by blank spaces in the buffer.
      do i = 1, len(dataBuffer)
        if (dataBuffer(i:i) == '(' .or. dataBuffer(i:i) == ')') dataBuffer(i:i) = ' '

      end do
      ! Read all coordinates in one go.
      read(dataBuffer, *) coords

    else
      do i = 1, nVertices
        read(unit, '(a)') lineBuffer
        read(lineBuffer(2:len_trim(lineBuffer) - 1), *) coords(:, i)

      end do

    end if

    close(unit)

    ! Initialise vertexShelf.
    call self % initVertexShelf(coords)

  end subroutine importVertices

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an unitialised state.
  !!
  subroutine kill(self)
    class(OpenFOAMMesh), intent(inout) :: self

    ! Call unstructuredMesh procedure.
    call kill_super(self)

  end subroutine kill

end module OpenFOAMMesh_class