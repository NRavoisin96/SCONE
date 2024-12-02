module OpenFOAMMesh_class_backup
  
  use numPrecision
  use universalVariables
  use genericProcedures,      only : fatalError, linFind, targetNotFound, numToChar, openToRead, &
                                     removeDuplicates, findCommon, findDifferent, hasDuplicates, &
                                     append
  use dictionary_class,       only : dictionary
  use face_class,             only : face
  use faceShelf_class,        only : faceShelf
  use edge_class,             only : edge
  use edgeShelf_class,        only : edgeShelf
  use element_class,          only : element
  use elementShelf_class,     only : elementShelf
  use kdTree_class,           only : kdTree
  use cellZone_class,         only : cellZone
  use cellZoneShelf_class,    only : cellZoneShelf
  use tetrahedron_class,      only : tetrahedron
  use tetrahedronShelf_class, only : tetrahedronShelf
  use pyramidShelf_class,     only : pyramidShelf
  use triangle_class,         only : triangle
  use triangleShelf_class,    only : triangleShelf
  use vertex_class,           only : vertex
  use vertexShelf_class,      only : vertexShelf
  use coord_class,            only : coord
  use mesh_inter,             only : mesh, kill_super => kill
  
  implicit none
  private
  
  !! OpenFOAM mesh. Uses a vertex -> face -> element representation of space. Each element is
  !! composed by a set of faces which are themselves composed by a number of vertices. OpenFOAM 
  !! meshes must consist of AT LEAST 3 files: points, faces and owner. If the mesh contains internal
  !! faces then a neighbour file must also be present.
  !!
  !! In OpenFOAM meshes elements can be grouped together in cell zones. This is useful to assign
  !! material filling in mesh elements. Local Ids are assigned in the order of the cell zone 
  !! definition given in the cellZones file. If this file is absent then only a single local Id can
  !! be assigned.
  !!
  !! Sample Input Dictionary:
  !!   mesh {id 7;
  !!         type OpenFOAMMesh;}
  !!
  !! Notes:
  !!   For tracking all OpenFOAM meshes are internally decomposed into tetrahedral meshes.
  !!
  !! Public members:
  !!   cellZones  -> Shelf that stores cell zones.
  !!   elements   -> Shelf that stores elements.
  !!   faces      -> Shelf that stores faces.
  !!   pyramids   -> Shelf that stores pyramids.
  !!   tetrahedra -> Shelf that stores tetrahedra.
  !!   triangles  -> Shelf that stores triangles.
  !!   vertices   -> Shelf that stores vertices.
  !!
  !! Private members:
  !!   tree -> kd-tree used for nearest-neighbour and entry checks.
  !!
  !! Interface:
  !!   Mesh interface.
  !!
  type, public, extends(mesh) :: OpenFOAMMesh
    private
    integer(shortInt)              :: nVertices = 0, nFaces = 0, nEdges = 0, &
                                      nElements = 0, nInternalFaces = 0, nTriangles = 0
    logical(defBool)               :: cellZonesFile = .false.
    type(cellZoneShelf), public    :: cellZones
!    type(edgeShelf), public        :: edges
    type(elementShelf), public     :: elements
    type(faceShelf), public        :: faces
    type(pyramidShelf), public     :: pyramids
    type(tetrahedronShelf), public :: tetrahedra
    type(triangleShelf), public    :: triangles
    type(vertexShelf), public      :: vertices
    type(kdTree)                   :: tree
  contains
    ! Build procedures
    procedure, private             :: checkFiles
    procedure, private             :: getMeshInfo
    procedure                      :: init
    procedure                      :: getNZones
    procedure                      :: kill
    procedure                      :: printComposition
    procedure                      :: split
    procedure                      :: computePrimitives
    procedure                      :: splitElements
    procedure                      :: splitPyramids
    ! Runtime procedures
    procedure                      :: checkForEntry
    procedure                      :: distanceToNextTriangle
    procedure                      :: distance
    procedure                      :: findElement
    procedure                      :: findTetrahedron
    procedure                      :: findTetrahedronFromEdge
    procedure                      :: findTetrahedronFromFace
    procedure                      :: findTetrahedronFromVertex
    procedure                      :: getNVertices
    procedure                      :: getNFaces
    procedure                      :: getNElements
    procedure                      :: getNInternalFaces
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
  subroutine checkFiles(self, folderPath, nInternalFaces)
    class(OpenFOAMMesh), intent(inout) :: self
    character(*), intent(in)           :: folderPath
    integer(shortInt), intent(in)      :: nInternalFaces
    logical(defBool)                   :: pointsFile, facesFile, neighbourFile
    character(100), parameter          :: Here = 'checkFiles (OpenFOAMMesh_class.f90)'

    ! Check the existence of the 'points' and 'faces' files and report errors if they are not found.
    inquire(file = folderPath//'points', exist = pointsFile)
    if (.not. pointsFile) call fatalError(Here, "Missing 'points' file for OpenFOAM mesh with Id: "//numToChar(self % id())//'.')
    
    inquire(file = folderPath//'faces', exist = facesFile)
    if (.not. facesFile) call fatalError(Here, "Missing 'faces' file for OpenFOAM mesh with Id: "//numToChar(self % id())//'.')
    
    ! If nInternal = 0, return early.
    if (nInternalFaces == 0) return

    ! If reached here check that the 'neighbour' file exists and report error if not.
    inquire(file = folderPath//'neighbour', exist = neighbourFile)
    if (.not. neighbourFile) call fatalError(Here, &
    "Missing 'neighbour' file for OpenFOAM mesh with Id: "//numToChar(self % id())//'.')

    ! Check that the 'cellZones' file exists and register information into mesh.
    inquire(file = folderPath//'cellZones', exist = self % cellZonesFile)
  end subroutine checkFiles
  
  !! Subroutine 'distance'
  !!
  !! See mesh_inter for details.
  !!
  pure subroutine distance(self, d, coords, isInside)
    class(OpenFOAMMesh), intent(in) :: self
    real(defReal), intent(out)      :: d
    type(coord), intent(inout)      :: coords
    logical(defBool), intent(out)   :: isInside

    ! Initialise isInside and compute the particle's end position.
    isInside = .true.
    
    ! If particle is already inside a tetrahedron, simply compute the distance to the next mesh triangle and return.
    if (coords % elementIdx > 0) then
      call self % distanceToNextTriangle(d, coords)
      return

    end if

    ! If not, we need to check if the particle enters the mesh.
    call self % checkForEntry(d, coords)

    ! If the particle enters the mesh retrieve the localId based on the element associated with the entry tetrahedron and return.
    if (coords % elementIdx > 0) then
      coords % localId = self % cellZones % findCellZone(self % tetrahedra % shelf(coords % elementIdx) % getElement())
      return

    end if
      
    ! If reached here, the particle does not enter the mesh and CSG tracking resumes.
    isInside = .false.
  end subroutine distance
  
  !! Subroutine 'findElement'
  !!
  !! See mesh_inter for details.
  !!
  pure subroutine findElement(self, r, u, elementIdx, localId)
    class(OpenFOAMMesh), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r, u
    integer(shortInt), intent(out)          :: elementIdx, localId
    
    ! Initialise localId = 1 (corresponds to the particle being in the CSG cell) and find the index of the tetrahedron
    ! in which the particle resides.
    localId = 1
    elementIdx = self % findTetrahedron(r, u)

    ! If tetrahedronIdx = 0 the particle is in the CSG cell. Return early in this case.
    if (elementIdx == 0) return

    ! Update localId based on the element associated with the tetrahedron.
    localId = self % cellZones % findCellZone(self % tetrahedra % shelf(elementIdx) % getElement())

  end subroutine findElement
  
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
  subroutine getMeshInfo(self, folderPath)
    class(OpenFOAMMesh), intent(inout) :: self
    character(*), intent(in)           :: folderPath
    logical(defBool)                   :: ownerFile
    integer(shortInt)                  :: colIdx, nIdx, endQuoteIdx
    integer(shortInt), parameter       :: unit = 10
    character(100)                     :: string, partialString
    character(:), allocatable          :: ownerPath
    character(*), parameter            :: Here = 'getMeshInfo (OpenFOAMMesh_class.f90)'

    ! Initialise ownerPath and perform a first check for the existence of the 'owner' file. Call fatalError if not found.
    ownerPath = folderPath//'owner'
    inquire(file = ownerPath, exist = ownerFile)
    if (.not. ownerFile) call fatalError(Here, "Missing 'owner' file for OpenFOAM mesh with Id: "//numToChar(self % id())//'.')
    
    ! Open the 'owner' file and read it. Skip lines until the line with the keyword 'note' is encountered.
    call openToRead(unit, ownerPath)
    read(unit, "(a)") string
    do while (index(string, 'note') == 0)
      read(unit, "(a)") string
    end do
    
    ! Locate the indices corresponding to the first ':' and ' ' symbols, and retrieve the variable
    ! 'nVertices', which is simply the piece of the text line located between the two aforementioned
    ! symbols.
    colIdx = index(string, ':')
    nIdx = index(string(colIdx:len_trim(string)), 'n')
    partialString = trim(adjustl(string(colIdx + 1:colIdx + nIdx)))
    read(partialString, *) self % nVertices
    
    ! Repeat for nElements.
    string(colIdx:colIdx) = ' '
    colIdx = index(string, ':')
    nIdx = index(string(colIdx:len_trim(string)), 'n')
    partialString = trim(adjustl(string(colIdx + 1:colIdx + nIdx)))
    read(partialString, *) self % nElements
    
    ! Repeat for nFaces.
    string(colIdx:colIdx) = ' '
    colIdx = index(string, ':')
    nIdx = index(string(colIdx:len_trim(string)), 'n')
    partialString = trim(adjustl(string(colIdx + 1:colIdx + nIdx)))
    read(partialString, *) self % nFaces
    
    ! Repeat for nInternalFaces. Close 'owner' file once done.
    string(colIdx:colIdx) = ' '
    colIdx = index(string, ':')
    endQuoteIdx = index(string(colIdx:len_trim(string)), '"')
    partialString = trim(adjustl(string(colIdx + 1:colIdx + endQuoteIdx - 2)))
    read(partialString, *) self % nInternalFaces
    close(unit)
    
    ! Now check existence of remaining mesh files.
    call self % checkFiles(folderPath, self % nInternalFaces)
  end subroutine getMeshInfo
  
  !! Subroutine 'init'
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
  !!   - fatalError if id < 1.
  !!   - fatalError if the mesh contains concave elements.
  !!
  subroutine init(self, folderPath, dict)
    class(OpenFOAMMesh), intent(inout)           :: self
    character(*), intent(in)                     :: folderPath
    class(dictionary), intent(in)                :: dict
    integer(shortInt)                            :: i, j, k, id, nConcaveElements, &
                                                    freeVertexIdx, faceIdx, vertexIdx
    integer(shortInt), dimension(:), allocatable :: vertexIdxs, faceIdxs
    type(vertex), dimension(:), allocatable      :: vertices
    type(face), dimension(:), allocatable        :: faces
    character(100), parameter                    :: Here = 'init (OpenFOAMMesh_class.f90)'
    
    ! Load id from the dictionary. Call fatal error if id is unvalid.
    call dict % get(id, 'id')
    if (id < 1) call fatalError(Here, 'Mesh Id must be +ve. Is: '//numToChar(id)//'.')
    
    ! Set mesh Id and retrieve preliminary information about the mesh.
    call self % setId(id)
    call self % getMeshInfo(folderPath)
    
    ! Allocate memory to the 'vertices', 'faces' and 'elements' structures.
    allocate(self % vertices % shelf(self % nVertices + self % nElements))
    allocate(self % faces % shelf(self % nFaces))
    allocate(self % elements % shelf(self % nElements))
    
    ! Import vertices and set mesh bounding box.
    call self % vertices % init(folderPath, self % nVertices)
    call self % setBoundingBox(self % vertices % getExtremalCoordinates())

    ! Import faces.
    call self % faces % init(folderPath, self % nInternalFaces)
    
    ! Set vertex-to-face connectivity information and compute the area and normal of each face.
    do i = 1, self % nFaces
      vertexIdxs = self % faces % shelf(i) % getVertices()
      do j = 1, size(vertexIdxs)
        call self % vertices % shelf(vertexIdxs(j)) % addFaceIdx(i)

      end do
      call self % faces % shelf(i) % computeAreaAndNormal(self % vertices % shelf(vertexIdxs))

    end do
    
    ! Import elements.
    call self % elements % init(folderPath, self % nElements, self % nFaces, self % nInternalFaces)
    
    ! Print original mesh composition.
    !    call self % printComposition(nTetrahedra)
    
    ! Initialise nConcaveElements = 0 and set vertex-to-element connectivity information.
    nConcaveElements = 0
    do i = 1, self % nElements
      faceIdxs = abs(self % elements % shelf(i) % getFaces())
      do j = 1, size(faceIdxs)
        faceIdx = faceIdxs(j)
        call self % faces % shelf(faceIdx) % addElementToFace(i)
        vertexIdxs = self % faces % shelf(faceIdx) % getVertices()
        do k = 1, size(vertexIdxs)
          vertexIdx = vertexIdxs(k)
          call self % vertices % shelf(vertexIdx) % addElementIdx(i)
          call self % elements % shelf(i) % addVertexToElement(vertexIdx)
        
        end do

      end do
      
      ! Compute the volume and centroid of current element.
      faces = self % faces % shelf(faceIdxs)
      vertices = self % vertices % shelf(self % elements % shelf(i) % getVertices())
      call self % elements % shelf(i) % computeVolumeAndCentroid(vertices, faces)

      ! Check for concave elements. Cycle if the current element is a tetrahedron.
      if (size(vertices) == 4) cycle
      if (.not. self % elements % shelf(i) % isConvex(vertices, faces)) nConcaveElements = nConcaveElements + 1

    end do

    ! If there are concave elements in the mesh call fatalError.
    if (nConcaveElements > 0) call fatalError(Here, numToChar(nConcaveElements)//' elements failed convexity test.')
    
    ! Split non-tetrahedral mesh elements into tetrahedra.
    call self % split(freeVertexIdx)
    
    ! Collapse the 'vertices' structure to save memory and initialise kd-tree for the mesh.
    call self % vertices % collapse(freeVertexIdx)
    call self % tree % init(self % vertices % getAllCoordinates(), .true.)

    ! If the 'cellZones' file is missing, allocate only one 'cellZones' structure and return.
    if (.not. self % cellZonesFile) then
      allocate(self % cellZones % shelf(1))
      call self % cellZones % shelf(1) % init('Default', 1, self % elements % shelf(self % nElements) % getIdx())
      return

    end if
    
    ! If reached here import cell zones.
    call self % cellZones % init(folderPath)
  end subroutine init
  
  !! Function 'getNZones'
  !!
  !! Basic description:
  !!   Returns the number of cell zones in the mesh.
  !!
  !! Result:
  !!   nZones -> Number of cell zones in the mesh.
  !!
  elemental function getNZones(self) result(nZones)
    class(OpenFOAMMesh), intent(in) :: self
    integer(shortInt)               :: nZones
    
    nZones = size(self % cellZones % shelf)
  end function getNZones
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(OpenFOAMMesh), intent(inout) :: self
    ! Superclass.
    call kill_super(self)
    ! Local.
    self % nVertices = 0
    self % nFaces = 0
    self % nInternalFaces = 0
    self % nElements = 0
    self % nEdges = 0
    self % nTriangles = 0
    self % cellZonesFile = .false.
    if (allocated(self % vertices % shelf)) call self % vertices % kill()
    if (allocated(self % faces % shelf)) call self % faces % kill()
    if (allocated(self % elements % shelf)) call self % elements % kill()
    if (allocated(self % cellZones % shelf)) call self % cellZones % kill()
    if (allocated(self % triangles % shelf)) call self % triangles % kill()
    if (allocated(self % pyramids % shelf)) call self % pyramids % kill()
    if (allocated(self % tetrahedra % shelf)) call self % tetrahedra % kill()
    call self % tree % kill()
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
    class(OpenFOAMMesh), intent(in) :: self
    integer(shortInt), intent(out)  :: nTetrahedra
    integer(shortInt)               :: nFaces, nPentahedra, nHexahedra, nOthers, i
    
    ! Initialise the numbers of various polyhedra to zero.
    nTetrahedra = 0
    nPentahedra = 0
    nHexahedra = 0
    nOthers = 0
    
    ! Loop over all elements in the mesh.
    do i = 1, self % nElements
      ! Retrieve the number of faces in the current element and increment specific polyhedra
      ! accordingly.
      nFaces = size(self % elements % shelf(i) % getFaces())
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
    print *, 'Displaying OpenFOAM mesh composition:'
    print *, '  Number of tetrahedra     : '//numToChar(nTetrahedra)//'.'
    print *, '  Number of pentahedra     : '//numToChar(nPentahedra)//'.'
    print *, '  Number of hexahedra      : '//numToChar(nHexahedra)//'.'
    print *, '  Number of other polyhedra: '//numToChar(nOthers)//'.'

  end subroutine printComposition
  
  !! Subroutine 'split'
  !!
  !! Basic description:
  !!   Splits a mesh into tetrahedral elements. If a given element is already a tetrahedron it is
  !!   not split but simply added to the shelf of tetrahedra in the mesh.
  !!
  !! Detailed description:
  !!   'split' starts by computing the number of pyramids, tetrahedra and triangles that will be
  !!   generated in the resulting mesh. Then, each element is split into a set of pyramids, whose
  !!   bases are each of the element's face and whose (common) apex is the element's centroid. This
  !!   apex is also appended to the list of vertices in the mesh in the process. Once this is done,
  !!   each face in the original mesh is subdivided into triangles. Lastly, each pyramid previously
  !!   created is further split into tetrahedra.
  !!
  !! Arguments:
  !!   lastVertexIdx [out] -> Index of the last vertex in the resulting mesh.
  !!
  elemental subroutine split(self, lastVertexIdx)
    class(OpenFOAMMesh), intent(inout) :: self
    integer(shortInt), intent(out)     :: lastVertexIdx
    integer(shortInt)                  :: nTriangles, nPyramids, nTetrahedra, &
                                          lastTriangleIdx, lastTetrahedronIdx, lastPyramidIdx
    
    ! Compute the number of pyramids, triangles and tetrahedra to be created and
    ! allocate memory to the corresponding structures.
    call self % computePrimitives(nPyramids, nTriangles, nTetrahedra)
    allocate(self % pyramids % shelf(nPyramids))
    allocate(self % triangles % shelf(nTriangles))
    allocate(self % tetrahedra % shelf(nTetrahedra))
    
    ! Initialise lastVertexIdx, lastPyramidIdx, lastTetrahedronIdx and lastTriangleIdx then
    ! split all elements into pyramids and all pyramids into tetrahedra.
    lastVertexIdx = self % nVertices
    lastPyramidIdx = 0
    lastTetrahedronIdx = 0
    lastTriangleIdx = 0
    call self % splitElements(lastVertexIdx, lastPyramidIdx, lastTetrahedronIdx, lastTriangleIdx)
    call self % splitPyramids(lastTetrahedronIdx, lastTriangleIdx)

  end subroutine split
  
  !! Subroutine 'computePrimitives'
  !!
  !! Basic description:
  !!   Computes the number of pyramids, triangles and tetrahedra to be created during the mesh
  !!   splitting process.
  !!
  !! Detailed description:
  !!   The number of pyramids is simply given by the sum of the number of faces in each element in
  !!   the original element. The number of triangles is more complex: each pyramid created during
  !!   the splitting process also creates a number of triangles equal to the number of edges (or
  !!   vertices) in the current face. However, since all these triangles are internal they are
  !!   always shared between two pyramids; therefore, the number of triangles to be generated during
  !!   the pyramid creation process is, for a given element, equal to the sum of the number of
  !!   vertices in each of the element's face divided by two. Triangles are also created during the
  !!   splitting of the original mesh's faces: for a given face, the number of triangles to be
  !!   created is simply equal to the number of vertices in the face less two. Lastly, during the
  !!   splitting of pyramids into tetrahedra, additional internal triangles are created, given by
  !!   the number of vertices in a given pyramid's base less three. The number of tetrahedra to be
  !!   generated simply is, for a given face, the number of triangles it is decomposed into.
  !!
  !! Arguments:
  !!   nPyramids [out]   -> Number of pyramids to be generated.
  !!   nTriangles [out]  -> Number of triangles to be generated.
  !!   nTetrahedra [out] -> Number of tetrahedra to be generated.
  !!
  elemental subroutine computePrimitives(self, nPyramids, nTriangles, nTetrahedra)
    class(OpenFOAMMesh), intent(in)              :: self
    integer(shortInt), intent(out)               :: nPyramids, nTriangles, nTetrahedra
    integer(shortInt)                            :: i, j, nVertices, nFaces, nVerticesInFace, &
                                                    absFaceIdx
    integer(shortInt), dimension(:), allocatable :: faceIdxs
    
    ! Initialise nPyramids, nTetrahedra and nTriangles and loop through all elements.
    nPyramids = 0
    nTetrahedra = 0
    nTriangles = 0
    do i = 1, self % nElements
      ! Retrieve the number of vertices and indices of the faces in the current element.
      nVertices = size(self % elements % shelf(i) % getVertices())
      faceIdxs = self % elements % shelf(i) % getFaces()
      
      ! Check if the current element is already a tetrahedron. If yes, increment nTetrahedra by 1
      ! and nTriangles by the number of triangles owned by the tetrahedron then cycle.
      if (nVertices == 4) then
        nTetrahedra = nTetrahedra + 1
        nTriangles = nTriangles + count(faceIdxs > 0)
        cycle
      
      end if   

      ! Compute the number of faces in the current element.
      nFaces = size(faceIdxs)
      
      ! Since there is one pyramid per element face, increase nPyramids by nFaces.
      nPyramids = nPyramids + nFaces
      
      ! Initialise nVertices and loop through all faces.
      nVertices = 0
      do j = 1, nFaces
        absFaceIdx = abs(faceIdxs(j))
        ! Retrieve the number of vertices in the current face and increase the total 
        ! number of vertices by the number of vertices in the current face.
        nVerticesInFace = size(self % faces % shelf(absFaceIdx) % getVertices())
        nVertices = nVertices + nVerticesInFace
        
        ! Increase the number of triangles corresponding to new internal faces by nVerticesInFace - 3.
        nTriangles = nTriangles + nVerticesInFace - 3

        ! If the element owns the current face, increase the number of triangles by nVerticesInFace - 2.
        if (faceIdxs(j) > 0) nTriangles = nTriangles + nVerticesInFace - 2
        
        ! There will be as many tetrahedra as the number of triangles in each face, which is given
        ! by nVerticesInFace - 2.
        nTetrahedra = nTetrahedra + nVerticesInFace - 2
      
      end do
      
      ! The number of pyramids' faces is given by half the total number of vertices.
      nTriangles = nTriangles + nVertices / 2
    
    end do

  end subroutine computePrimitives
  
  !! Subroutine 'splitElements'
  !!
  !! Basic description:
  !!   Splits all elements in the original mesh into pyramids. If a given element is already a
  !!   tetrahedron it is not split but simply appended to the list of existing tetrahedra.
  !!
  !! Arguments:
  !!   lastVertexIdx [inout]      -> Index of the last vertex in the mesh.
  !!   lastPyramidIdx [inout]     -> Index of the last pyramid in the mesh.
  !!   lastTetrahedronIdx [inout] -> Index of the last tetrahedron in the mesh.
  !!   lastTriangleIdx [inout]    -> Index of the last triangle in the mesh.
  !!
  elemental subroutine splitElements(self, lastVertexIdx, lastPyramidIdx, lastTetrahedronIdx, lastTriangleIdx)
    class(OpenFOAMMesh), intent(inout)           :: self
    integer(shortInt), intent(inout)             :: lastVertexIdx, lastPyramidIdx, lastTetrahedronIdx, lastTriangleIdx
    integer(shortInt)                            :: i, j, initialTriangleIdx
    integer(shortInt), dimension(:), allocatable :: faces, vertices, triangles
    real(defReal), dimension(3)                  :: centroid
                                        
    ! Loop through all elements.
    do i = 1, self % nElements
      ! Retrieve the indices of the vertices and the faces in the current element as well as the centroid
      ! of the element.
      vertices = self % elements % shelf(i) % getVertices()
      faces = self % elements % shelf(i) % getFaces()
      centroid = self % elements % shelf(i) % getCentroid()

      ! If the current element is not a tetrahedron we need to split it. Update lastVertexIdx and add
      ! the current element's centroid to the shelf.
      if (size(vertices) > 4) then
        lastVertexIdx = lastVertexIdx + 1
        call self % vertices % shelf(lastVertexIdx) % setIdx(lastVertexIdx)
        call self % vertices % shelf(lastVertexIdx) % setCoordinates(centroid)
        
        ! Set initialTriangleIdx, split the current element into pyramids and cycle to the next element.
        initialTriangleIdx = lastTriangleIdx + 1
        call self % elements % shelf(i) % split(self % faces, self % vertices, self % triangles, &
                                                self % pyramids, lastTriangleIdx, lastPyramidIdx, &
                                                lastVertexIdx)
        cycle

      end if
      
      ! If the current element is already a tetrahedron then there is no need to split it: simply
      ! create a new tetrahedron in the shelf.
      lastTetrahedronIdx = lastTetrahedronIdx + 1
      call self % tetrahedra % shelf(lastTetrahedronIdx) % setIdx(lastTetrahedronIdx)
      call self % tetrahedra % shelf(lastTetrahedronIdx) % setVertices(vertices)
      call self % tetrahedra % shelf(lastTetrahedronIdx) % setElement(i)
      call self % tetrahedra % shelf(lastTetrahedronIdx) % setCentroid(centroid)
      
      ! Loop through all the faces of the new tetrahedron.
      do j = 1, 4
        ! Add the tetrahedron to the vertex.
        call self % vertices % shelf(vertices(j)) % addTetrahedronIdx(lastTetrahedronIdx)
        
        ! If faces(j) > 0 then the tetrahedron owns the current face and we need to add a 
        ! new triangle corresponding to this face to the shelf.
        if (faces(j) > 0) then
          ! Generate a new triangle from the current face and update mesh connectivity information.
          call self % faces % shelf(faces(j)) % split(self % triangles, self % vertices, lastTriangleIdx)
          call self % triangles % shelf(lastTriangleIdx) % addTetrahedronIdx(lastTetrahedronIdx)
          call self % tetrahedra % shelf(lastTetrahedronIdx) % addTriangle(lastTriangleIdx)
          cycle

        end if

        ! If faces(j) < 0 then a new triangle has already been generated for the current face. 
        ! Retrieve the triangle corresponding to the face and set mesh connectivity information.
        triangles = self % faces % shelf(abs(faces(j))) % getTriangles()
        call self % tetrahedra % shelf(lastTetrahedronIdx) % addTriangle(-triangles(1))
        call self % triangles % shelf(triangles(1)) % addTetrahedronIdx(lastTetrahedronIdx)

      end do
      
    end do

  end subroutine splitElements
  
  !! Subroutine 'splitPyramids'
  !!
  !! Basic description:
  !!   Splits all pyramids created during the splitting of original elements into tetrahedra. Also
  !!   splits the base of the pyramid into triangles if the pyramid owns the base.
  !!
  !! Arguments:
  !!   lastTetrahedronIdx [inout] -> Index of the last tetrahedron in the mesh.
  !!   lastTriangleIdx [inout]    -> Index of the last triangle in the mesh.
  !!
  elemental subroutine splitPyramids(self, lastTetrahedronIdx, lastTriangleIdx)
    class(OpenFOAMMesh), intent(inout)           :: self
    integer(shortInt), intent(inout)             :: lastTetrahedronIdx, lastTriangleIdx
    integer(shortInt)                            :: i, j, faceIdx, initialTriangleIdx, &
                                                    initialTetrahedronIdx
    integer(shortInt), dimension(:), allocatable :: triangles
    
    ! Loop over all the pyramids.
    do i = 1, size(self % pyramids % shelf)
      ! Retrieve the face making the pyramid's base.
      faceIdx = self % pyramids % shelf(i) % getFace()
      
      ! If the pyramid owns the face, split it into triangles.
      if (faceIdx > 0) then
        ! Set initialTriangleIdx and split the face.
        initialTriangleIdx = lastTriangleIdx + 1
        call self % faces % shelf(faceIdx) % split(self % triangles, self % vertices, lastTriangleIdx)

        ! Add all new triangles to the current pyramid.
        do j = initialTriangleIdx, lastTriangleIdx
          call self % pyramids % shelf(i) % addTriangle(j)

        end do

      else
        ! If the pyramid does not own the face then it has already been subdivided. Retrieve the 
        ! triangles in the current face and add each triangle to the pyramid.
        triangles = self % faces % shelf(abs(faceIdx)) % getTriangles()
        do j = 1, size(triangles)
          call self % pyramids % shelf(i) % addTriangle(-triangles(j))

        end do

      end if

      ! Update initialTriangleIdx and initialTetrahedronIdx then split the pyramid into tetrahedra.
      initialTriangleIdx = lastTriangleIdx + 1
      initialTetrahedronIdx = lastTetrahedronIdx + 1
      call self % pyramids % shelf(i) % split(self % faces % shelf(abs(faceIdx)), self % vertices, &
                                              self % triangles, self % tetrahedra, &
                                              lastTriangleIdx, lastTetrahedronIdx)

    end do

  end subroutine splitPyramids
  
  !! Subroutine 'checkForEntry'
  !!
  !! Basic description:
  !!   Checks whether a given particle enters from the CSG cell into the mesh. Returns the distance
  !!   to the triangle of intersection in case the particle enters the mesh. Returns INF 
  !!   otherwise.
  !!
  !! Arguments:
  !!   dist [out]     -> Distance to the intersected triangle.
  !!   coords [inout] -> Particle's coordinates.
  !!
  pure subroutine checkForEntry(self, dist, coords)
    class(OpenFOAMMesh), intent(in) :: self
    real(defReal), intent(out)      :: dist
    type(coord), intent(inout)      :: coords
    real(defReal)                   :: rComponent, rEndComponent
    real(defReal), dimension(6)     :: boundingBox
    integer(shortInt)               :: i
    
    ! Initialise distance to INF and check that the particle's path intersects the mesh's bounding box.
    dist = INF
    boundingBox = self % getBoundingBox()
    do i = 1, 3
      ! If both the particle's initial and end locations are on the same side of the plane of the
      ! bounding box return early.
      rComponent = coords % r(i)
      rEndComponent = coords % rEnd(i)
      if ((rComponent < boundingBox(i) .and. rEndComponent < boundingBox(i)) .or. &
          (rComponent > boundingBox(i + 3) .and. rEndComponent > boundingBox(i + 3))) return

    end do
    
    ! If reached here, search the tree for the vertices contained in potential intersected triangles.
    call self % tree % findIntersectedTriangle(dist, coords, self % vertices, self % triangles)

  end subroutine checkForEntry
  
  !! Subroutine 'distanceToNextTriangle'
  !!
  !! Basic description:
  !!   Returns the distance to the next triangle intersected by the particle's path. Returns
  !!   INF if the particle does not intersect any triangle (i.e., if its path is entirely
  !!   contained in the tetrahedron the particle currently is). Algorithm adapted from Macpherson, 
  !!   et al. (2009). DOI: 10.1002/cnm.1128.
  !!
  !! Arguments:
  !!   dist [out]     -> Distance to the next intersected triangle.
  !!   coords [inout] -> Particle's coordinates.
  !!
  !! TODO: Optimise this.
  pure subroutine distanceToNextTriangle(self, dist, coords)
    class(OpenFOAMMesh), intent(in)              :: self
    real(defReal), intent(out)                   :: dist
    type(coord), intent(inout)                   :: coords
    real(defReal), dimension(3)                  :: r, rEnd
    integer(shortInt)                            :: intersectedTriangleIdx, oldElement, newElement
    integer(shortInt), dimension(:), allocatable :: potentialTriangles, triangleToTetrahedra
    real(defReal)                                :: lambda
    type(tetrahedron)                            :: currentTetrahedron
    type(triangle)                               :: intersectedTriangle
    
    ! Initialise dist to INF, retrieve the tetrahedron currently occupied 
    ! by the particle and compute potential triangle intersections.
    dist = INF
    currentTetrahedron = self % tetrahedra % shelf(coords % elementIdx)
    rEnd = coords % rEnd
    call currentTetrahedron % computePotentialTriangles(rEnd, self % triangles, potentialTriangles)
    
    ! If no potential intersections are detected return early.
    if (.not. allocated(potentialTriangles)) return
    
    ! If reached here, compute which triangle is actually intersected and update dist.
    r = coords % r
    call currentTetrahedron % computeIntersectedTriangle(r, rEnd, potentialTriangles, &
                                                         intersectedTriangleIdx, lambda, &
                                                         self % triangles)
    dist = norm2(min(ONE, max(ZERO, lambda)) * (rEnd - r))
    
    ! If the intersected triangle is a boundary triangle then the particle is leaving the mesh.
    intersectedTriangle = self % triangles % shelf(intersectedTriangleIdx)
    if (intersectedTriangle % getIsBoundary()) then
      coords % elementIdx = 0
      coords % localId = 1
      return

    end if

    ! Else, retrieve the tetrahedra sharing the intersected triangle from mesh connectivity then
    ! update tetrahedronIdx and localId (if necessary).
    triangleToTetrahedra = intersectedTriangle % getTetrahedra()
    coords % elementIdx = findDifferent(triangleToTetrahedra, currentTetrahedron % getIdx())
    oldElement = currentTetrahedron % getElement()
    newElement = self % tetrahedra % shelf(coords % elementIdx) % getElement()
    if (newElement /= oldElement) coords % localId = self % cellZones % findCellZone(newElement)

  end subroutine distanceToNextTriangle
  
  !! Function 'findTetrahedron'
  !!
  !! Basic description:
  !!   Returns the index of the tetrahedron occupied by a particle.
  !!
  !! Arguments:
  !!   coords [in] -> Particle's coordinates.
  !!   dir [in]    -> Particle's direction.
  !!
  pure function findTetrahedron(self, coords, dir) result(tetrahedronIdx)
    class(OpenFOAMMesh), intent(in)              :: self
    real(defReal), dimension(3), intent(in)      :: coords, dir
    integer(shortInt)                            :: failedTriangle, i, nearestVertex, &
                                                    tetrahedronIdx
    integer(shortInt), dimension(:), allocatable :: potentialTetrahedra, triangleToTetrahedra, &
                                                    zeroDotProductTriangles
    real(defReal), dimension(6)                  :: boundingBox
    type(tetrahedron)                            :: currentTetrahedron
    
    ! Initialise tetrahedronIdx, retrieve the mesh's bounding box and check that the particle is 
    ! inside the mesh's bounding box. If not we can return early.
    tetrahedronIdx = 0
    boundingBox = self % getBoundingBox()
    do i = 1, 3
      if (coords(i) < boundingBox(i) .or. boundingBox(i + 3) < coords(i)) return

    end do

    ! If the point is inside the bounding box then we need to determine if the point is inside a
    ! pseudoCell. First find the vertex which is nearest to the coordinates.
    nearestVertex = self % tree % findNearestVertex(coords)
    
    ! Retrieve potential tetrahedra occupied by the particle from mesh connectivity information.
    potentialTetrahedra = self % vertices % shelf(nearestVertex) % getVertexToTetrahedra()
    
    ! Initialise the search to the first element in the potential tetrahedra.
    currentTetrahedron = self % tetrahedra % shelf(potentialTetrahedra(1))
    searchLoop: do
      ! Check if the particle is inside the current tetrahedron.
      call currentTetrahedron % testForInclusion(self % triangles, coords, failedTriangle, &
                                                 zeroDotProductTriangles)
      ! If the current tetrahedron is not occupied by the particle, retrieve the tetrahedron sharing
      ! the triangle for which the search failed and update the search.
      if (failedTriangle > 0) then
        triangleToTetrahedra = self % triangles % shelf(failedTriangle) % getTetrahedra()

        ! If the failed triangle is a boundary triangle the particle is outside the mesh.
        if (size(triangleToTetrahedra) == 1) return
        
        ! Update the tetrahedron to be searched and cycle.
        tetrahedronIdx = findDifferent(triangleToTetrahedra, currentTetrahedron % getIdx())
        currentTetrahedron = self % tetrahedra % shelf(tetrahedronIdx)
        cycle searchLoop

      end if
      
      ! If there are faces on which the particle lies we need to employ some more specific
      ! procedures to correctly determine which tetrahedron is actually occupied.
      if (allocated(zeroDotProductTriangles)) then
        select case (size(zeroDotProductTriangles))
          ! Particle is on a tetrahedron's face.
          case(1)
            tetrahedronIdx = self % findTetrahedronFromFace(dir, &
                                                            currentTetrahedron % getIdx(), &
                                                            zeroDotProductTriangles)
          ! Particle is on a tetrahedron's edge.
          case(2)
            tetrahedronIdx = self % findTetrahedronFromEdge(dir, zeroDotProductTriangles)
          ! Particle is on a tetrahedron's vertex.
          case(3)
            tetrahedronIdx = self % findTetrahedronFromVertex(dir, zeroDotProductTriangles)                                             
        end select
        return
      
      end if
      ! If reached this point the particle is in the current tetrahedron. Retrieve its index and 
      ! return.
      tetrahedronIdx = currentTetrahedron % getIdx()
      return

    end do searchLoop

  end function findTetrahedron

  !! Function 'findTetrahedronFromEdge'
  !!
  !! Basic description:
  !!   Returns the index of the tetrahedron occupied by the particle in case the particle lies on
  !!   a tetrahedron edge.
  !!
  !! Detailed description:
  !!   Generalisation of 'findTetrahedronFromFace' applied to an edge. In this case the tetrahedron
  !!   assigned to the particle is that for which the number of negative dot products between the
  !!   tetrahedron's triangles containing the edge and the particle's direction is the greatest.
  !!
  !! Arguments:
  !!   direction [in]               -> Particle's direction.
  !!   zeroDotProductTriangles [in] -> Indices of the triangles on which the particle lies.
  !!
  !! Result:
  !!   tetrahedronIdx -> Index of the tetrahedron occupied by the particle.
  !!
  !! Notes:
  !!   The implementation here is quite ugly and messy. This is because OpenFOAM does not store
  !!   information about edges, so at the moment they have to be constructed on-the-fly. Could be
  !!   made more elegant / faster if there was a procedure to build edges during the importation 
  !!   process, but given how unlikely a particle is to actually sit on an edge, I am not sure this 
  !!   would make much difference...
  !!
  !! TODO: Rework this. Import edges during mesh construction.
  pure function findTetrahedronFromEdge(self, direction, zeroDotProductTriangles) result(tetrahedronIdx)
    class(OpenFOAMMesh), intent(in)                          :: self
    real(defReal), dimension(3), intent(in)                  :: direction
    integer(shortInt), dimension(:), allocatable, intent(in) :: zeroDotProductTriangles
    integer(shortInt)                                        :: i, j, tetrahedronIdx
    integer(shortInt), dimension(:), allocatable             :: absZeroDotProductTriangles, &
                                                                edgeVertices, nTriangles, &
                                                                potentialTetrahedra, &
                                                                absTetrahedronTriangles, &
                                                                tetrahedronTriangles, commonVertices, &
                                                                verticesTriangle1, verticesTriangle2
    type(tetrahedron)                                        :: currentTetrahedron
    type(triangle)                                           :: currentTriangle
    real(defReal), dimension(3)                              :: normal
    ! First create absolute indices.
    absZeroDotProductTriangles = abs(zeroDotProductTriangles)
    ! Now retrieve the vertices making each triangle up.
    verticesTriangle1 = self % triangles % shelf(absZeroDotProductTriangles(1)) % getVertices()
    verticesTriangle2 = self % triangles % shelf(absZeroDotProductTriangles(2)) % getVertices()
    ! Retrieve the vertices making the common edge up.
    edgeVertices = findCommon(verticesTriangle1, verticesTriangle2)
    ! Loop through all triangles.
    do i = 1, size(self % triangles % shelf)
      ! If the index of the 'do' loop corresponds to one of the triangles sharing the common edge
      ! there is no need to search and we cycle.
      if (any(absZeroDotProductTriangles == i)) cycle
      currentTriangle = self % triangles % shelf(i)
      ! Check if the current face contains the common edge. The implementation below is a bit ugly:
      ! conceptually one could do away with the first if-loop and only retain the if (size(...)).
      ! However this causes problems for reasons which are unkill.
      if (hasDuplicates([currentTriangle % getVertices(), edgeVertices])) then
        ! Retrieve the common vertices between the current triangle and the common edge.
        commonVertices = findCommon(currentTriangle % getVertices(), edgeVertices)
        ! Check if the current triangle contains the common edge.
        if (size(findCommon(currentTriangle % getVertices(), edgeVertices)) == 2) then
          ! If so, append the current face's index to the list of faces.
          call append(absZeroDotProductTriangles, i)
        end if
      end if
    end do
    ! Loop through all the triangles that contain the common edge and retrieve their associated
    ! tetrahedra.
    do i = 1, size(absZeroDotProductTriangles)
      currentTriangle = self % triangles % shelf(absZeroDotProductTriangles(i))
      call append(potentialTetrahedra, currentTriangle % getTetrahedra()) 
    end do
    ! We might introduce duplicates in the list of potential tetrahedra so remove them.
    potentialTetrahedra = removeDuplicates(potentialTetrahedra)
    ! Allocate nTriangles and initialise the array to zero.
    allocate(nTriangles(size(potentialTetrahedra)))
    nTriangles = 0
    ! Loop through all potential tetrahedra.
    do i = 1, size(potentialTetrahedra)
      ! Retrieve the triangles making the current tetrahedron and convert them to absolute indices.
      currentTetrahedron = self % tetrahedra % shelf(potentialTetrahedra(i))
      tetrahedronTriangles = currentTetrahedron % getTriangles()
      absTetrahedronTriangles = abs(tetrahedronTriangles)
      ! Loop through all triangles.
      do j = 1, 4
      
        ! Check if the current triangle contains the common edge.
        if(any(absZeroDotProductTriangles == absTetrahedronTriangles(j))) then
        
          currentTriangle = self % triangles % shelf(absTetrahedronTriangles(j))
          ! If so, create the triangle's normal vector.
          normal = currentTriangle % getNormal(tetrahedronTriangles(j))
          ! Check whether the direction of the particle points inwards with respect to the current
          ! triangle.
          if (dot_product(direction, normal) <= ZERO) then
            ! Update nFaces for the current tetrahedron.
            nTriangles(i) = nTriangles(i) + 1
          ! If the current triangle is a boundary triangle, set tetrahedronIdx to zero and return.
          else if (size(currentTriangle % getTetrahedra()) == 1) then
            tetrahedronIdx = 0
            return
          end if
        end if
      end do
    end do
    ! if reached this point, the tetrahedron the particle resides in is the one for which nTriangles 
    ! is the greatest.
    tetrahedronIdx = potentialTetrahedra(maxloc(nTriangles, 1))
  end function findTetrahedronFromEdge
  
  !! Function 'findTetrahedronFromFace'
  !!
  !! Basic description:
  !!   Returns the index of the tetrahedron occupied by the particle in case the particle lies on
  !!   a tetrahedron face.
  !!
  !! Detailed description:
  !!   When a particle is on a tetrahedron face it is not as straightforward to assign a tetrahedron
  !!   to a particle (it is actually in both tetrahedra at the same time). The workaround here is to
  !!   use the particle's direction to determine into which tetrahedron the particle's direction 
  !!   points. This is given by the dot product of the direction and the triangle's normal being
  !!   negative.
  !!
  !! Arguments:
  !!   direction [in]               -> Particle's direction.
  !!   zeroDotProductTriangles [in] -> Indices of the triangles on which the particle lies.
  !!
  !! Result:
  !!   tetrahedronIdx -> Index of the tetrahedron occupied by the particle.
  !!
  pure function findTetrahedronFromFace(self, direction, currentTetrahedronIdx, &
                                        zeroDotProductTriangles) result(tetrahedronIdx)
    class(OpenFOAMMesh), intent(in)                          :: self
    real(defReal), dimension(3), intent(in)                  :: direction
    integer(shortInt), intent(in)                            :: currentTetrahedronIdx
    integer(shortInt), dimension(:), allocatable, intent(in) :: zeroDotProductTriangles
    integer(shortInt)                                        :: absZeroDotProductTriangle, &
                                                                tetrahedronIdx
    integer(shortInt), dimension(:), allocatable             :: triangleToTetrahedra
    real(defReal), dimension(3)                              :: normal
    
    ! Initialise tetrahedronIdx = 0. This corresponds to the particle being outside the mesh.
    tetrahedronIdx = 0
    
    ! Create absolute indices and retrieve the triangle's normal vector.
    absZeroDotProductTriangle = abs(zeroDotProductTriangles(1))
    normal = self % triangles % shelf(absZeroDotProductTriangle) % getNormal(zeroDotProductTriangles(1))
    
    ! Check the sign of the dot product between the particle's direction and the triangle's normal.
    ! If negative it is in the current tetrahedron and we can return early.
    if (dot_product(direction, normal) <= ZERO) then
      tetrahedronIdx = currentTetrahedronIdx
      return

    end if

    ! If not, the particle's in the neighbouring tetrahedron. Retrieve mesh connectivity information.
    triangleToTetrahedra = self % triangles % shelf(absZeroDotProductTriangle) % getTetrahedra()

    ! If the current triangle is not a boundary triangle, update tetrahedronIdx.
    if (size(triangleToTetrahedra) > 1) tetrahedronIdx = findDifferent(triangleToTetrahedra, currentTetrahedronIdx)
  end function findTetrahedronFromFace
  
  !! Function 'findTetrahedronFromVertex'
  !!
  !! Basic description:
  !!   Returns the index of the tetrahedron occupied by the particle in case the particle lies on
  !!   a tetrahedron vertex.
  !!
  !! Detailed description:
  !!   Generalisation of 'findTetrahedronFromFace' applied to a vertex. In this case the tetrahedron
  !!   assigned to the particle is that for which the number of negative dot products between the
  !!   tetrahedron's triangles containing the vertex and the particle's direction is the greatest.
  !!
  !! Arguments:
  !!   direction [in]               -> Particle's direction.
  !!   zeroDotProductTriangles [in] -> Indices of the triangles on which the particle lies.
  !!
  !! Result:
  !!   tetrahedronIdx -> Index of the tetrahedron occupied by the particle.
  !!
  pure function findTetrahedronFromVertex(self, direction, zeroDotProductTriangles) &
                result(tetrahedronIdx)
    class(OpenFOAMMesh), intent(in)                          :: self
    real(defReal), dimension(3), intent(in)                  :: direction
    integer(shortInt), dimension(:), allocatable, intent(in) :: zeroDotProductTriangles
    integer(shortInt)                                        :: commonVertex, i, j, tetrahedronIdx
    integer(shortInt), dimension(:), allocatable             :: absTriangles, &
                                                                absZeroDotProductTriangles, &
                                                                commonVertices, nTriangles, &
                                                                potentialTetrahedra, triangles, &
                                                                verticesTriangle1, &
                                                                verticesTriangle2, &
                                                                verticesTriangle3
    type(triangle)                                           :: currentTriangle
    real(defReal), dimension(3)                              :: normal
    
    ! Initialise tetrahedronIdx = 0. This corresponds to the particle being outside the mesh.
    tetrahedronIdx = 0

    ! First create absolute indices and retrieve the vertices making each of the triangles.
    absZeroDotProductTriangles = abs(zeroDotProductTriangles)
    verticesTriangle1 = self % triangles % shelf(absZeroDotProductTriangles(1)) % getVertices()
    verticesTriangle2 = self % triangles % shelf(absZeroDotProductTriangles(2)) % getVertices()
    verticesTriangle3 = self % triangles % shelf(absZeroDotProductTriangles(3)) % getVertices()
    
    ! Find the common vertex between the three triangles. Retrieve all tetrahedra sharing this vertex.
    commonVertices = findCommon(verticesTriangle1, verticesTriangle2)
    commonVertices = findCommon(commonVertices, verticesTriangle3)
    commonVertex = commonVertices(1)
    potentialTetrahedra = self % vertices % shelf(commonVertex) % getVertexToTetrahedra()
    
    ! Allocate nTriangles and initialise the array to zero.
    allocate(nTriangles(size(potentialTetrahedra)))
    nTriangles = 0
    
    ! Loop through all potential tetrahedra.
    do i = 1, size(potentialTetrahedra)
      ! Retrieve the triangles in the current tetrahedron and create absolute indices.
      triangles = self % tetrahedra % shelf(potentialTetrahedra(i)) % getTriangles()
      absTriangles = abs(triangles)
      
      ! Loop through all triangles.
      do j = 1, 4
        currentTriangle = self % triangles % shelf(absTriangles(j))
        ! Cycle if the current triangle does not contain the common vertex.
        if (.not. any(currentTriangle % getVertices() == commonVertex)) cycle
        
        ! Retrieve the current triangle's normal vector and flip it if necessary.
        normal = currentTriangle % getNormal(triangles(j))
        
        ! If the dot product is negative increment the number of triangles for the current 
        ! tetrahedron and cycle to the next triangle.
        if (dot_product(direction, normal) <= ZERO) then
          nTriangles(i) = nTriangles(i) + 1
          cycle
  
        end if
        
        ! If the test fails and the current triangle is a boundary triangle, the particle is
        ! outside the mesh and we can return early.
        if (size(currentTriangle % getTetrahedra()) == 1) return

      end do

    end do
    
    ! If reached here, the tetrahedron occupied by the particle is that for which nTriangles is the greatest.
    tetrahedronIdx = potentialTetrahedra(maxloc(nTriangles, 1))
  end function findTetrahedronFromVertex
  
  !! Function 'getNVertices'
  !!
  !! Basic description:
  !!   Returns the number of original (unsplit) mesh vertices.
  !!
  !! Result:
  !!   nVertices -> Number of vertices in the original mesh.
  !!
  elemental function getNVertices(self) result(nVertices)
    class(OpenFOAMMesh), intent(in) :: self
    integer(shortInt)               :: nVertices
    
    nVertices = self % nVertices
  end function getNVertices
  
  !! Function 'getNFaces'
  !!
  !! Basic description:
  !!   Returns the number of original (unsplit) mesh faces.
  !!
  !! Result:
  !!   nFaces -> Number of faces in the original mesh.
  !!
  elemental function getNFaces(self) result(nFaces)
    class(OpenFOAMMesh), intent(in) :: self
    integer(shortInt)               :: nFaces
    
    nFaces = self % nFaces
  end function getNFaces  
  
  !! Function 'getNElements'
  !!
  !! Basic description:
  !!   Returns the number of original (unsplit) mesh elements.
  !!
  !! Result:
  !!   nElements -> Number of elements in the original mesh.
  !!
  elemental function getNElements(self) result(nElements)
    class(OpenFOAMMesh), intent(in) :: self
    integer(shortInt)               :: nElements
    
    nElements = self % nElements
  end function getNElements
  
  !! Function 'getNInternalFaces'
  !!
  !! Basic description:
  !!   Returns the number of original (unsplit) mesh internal faces.
  !!
  !! Result:
  !!   nInternalFaces -> Number of internal faces in the original mesh.
  !!
  elemental function getNInternalFaces(self) result(nInternalFaces)
    class(OpenFOAMMesh), intent(in) :: self
    integer(shortInt)               :: nInternalFaces
    
    nInternalFaces = self % nInternalFaces
  end function getNInternalFaces
end module OpenFOAMMesh_class_backup