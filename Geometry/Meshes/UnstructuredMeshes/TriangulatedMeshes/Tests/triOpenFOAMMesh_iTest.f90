module triOpenFOAMMesh_iTest
  
  use coord_class,           only : coord
  use triOpenFOAMMesh_class, only : triOpenFOAMMesh
  use numPrecision
  use genericProcedures
  use dictionary_class,      only : dictionary
  use dictParser_func,       only : charToDict
  use funit
  use universalVariables
  
  implicit none
  
  ! Parameters.
  character(*), parameter :: MESH_DEF = &
  " id 2; type triOpenFOAMMesh; path ./IntegrationTestFiles/Geometry/Meshes/OpenFOAM/testMesh/;"
  ! Variables.
  type(triOpenFOAMMesh)   :: mesh
  type(coord)             :: coords

contains
  
  !!
  !! Import the mesh.
  !!
@Before
  subroutine setUp()
    type(dictionary)   :: dict
    character(pathLen) :: path
    
    call charToDict(dict, MESH_DEF)
    call dict % get(path, 'path')
    call mesh % init(trim(path), dict)
  
  end subroutine setUp
  
  !!
  !! Clean after tests.
  !!
@After
  subroutine cleanUp()
    call mesh % kill()

  end subroutine cleanUp
  
  !!
  !! Test miscellaneous functionality.
  !!
@Test
  subroutine test_misc()
    
    ! Test id.
    @assertEqual(2, mesh % getId())
    call mesh % setId(7)
    @assertEqual(7, mesh % getId())

  end subroutine test_misc
  
  !!
  !! Test mesh information.
  !!
@Test
  subroutine test_info()
    real(defReal) :: TOL = 1.0E-6
    
    ! Test number of vertices.
    @assertEqual(18, mesh % nVertices)
    ! Test number of faces.
    @assertEqual(20, mesh % nFaces)
    ! Test number of elements.
    @assertEqual(4, mesh % nElements)
    ! Test number of internal faces.
    @assertEqual(4, mesh % nInternalFaces)
    ! Test area of two faces.
    @assertEqual(2.0_defReal, mesh % faces % shelf(1) % getArea(), 2.0_defReal * TOL)
    @assertEqual(1.0_defReal, mesh % faces % shelf(11) % getArea(), 1.0_defReal * TOL)
    ! Test volume of one element.
    @assertEqual(2.0_defReal, mesh % elements % shelf(3) % getVolume(), 2.0_defReal * TOL)

  end subroutine test_info
  
  !!
  !! Test decomposition.
  !!
@Test
  subroutine test_decomposition()
    
    ! Test number of vertices after decomposition.
    @assertEqual(22, size(mesh % vertices % shelf))
    ! Test number of triangles after decomposition.
    @assertEqual(112, size(mesh % triangles % shelf))
    ! Test number of tetrahedra after decomposition.
    @assertEqual(48, size(mesh % tetrahedra % shelf))

  end subroutine test_decomposition
  
  !!
  !! Test inside / outside determination.
  !!
@Test
  subroutine test_inside()
    real(defReal), dimension(3) :: r, u
    integer(shortInt)           :: elementIdx, parentIdx
    
    ! Few points inside.
    r = [0.31_defReal, 0.42_defReal, 0.13_defReal]
    u = [ONE, ZERO, ZERO]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(36, elementIdx)
    r = [0.02_defReal, 0.97_defReal, -0.5_defReal]
    u = [ZERO, ONE, ZERO]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(35, elementIdx)
    
    ! Few points outside.
    r = [1.2_defReal, 0.8_defReal, 0.0_defReal]
    u = [ZERO, ONE, ZERO]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)
    r = [0.1_defReal, 0.1_defReal, 1.13_defReal]
    u = [ZERO, ZERO, -ONE]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)
    
    ! Few more difficult points.
    ! A point on a face. Points into the mesh.
    r = [-1.0_defReal, 0.1_defReal, 0.1_defReal]
    u = [ONE, ZERO, ZERO]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(15, elementIdx)
    
    ! Same point but points away from the mesh.
    u = [-ONE, ZERO, ZERO]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)
    
    ! A point on an internal vertex. Different directions.
    r = [0.0_defReal, 0.0_defReal, 0.0_defReal]
    u = [ONE, ZERO, ZERO]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(48, elementIdx)
    u = [-ONE, ZERO, ZERO]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(2, elementIdx)
    u = [ZERO, ONE, ZERO]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(14, elementIdx)
    u = [ZERO, -ONE, ZERO]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(45, elementIdx)
    u = [ZERO, ZERO, ONE]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(45, elementIdx)
    u = [ZERO, ZERO, -ONE]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(45, elementIdx)

  end subroutine test_inside
  
  !!
  !! Test distance calculations.
  !!
@Test
  subroutine test_distance()
    real(defReal)               :: maxDist, distance
    integer(shortInt)           :: parentIdx
    real(defReal), parameter    :: TOL = 1.0E-6
    maxDist = 0.5_defReal
    
    ! Few points inside mesh.
    coords % r = [0.98_defReal, 0.1_defReal, 0.1_defReal]
    coords % dir = [ONE, ZERO, ZERO]
    coords % rEnd = coords % r + coords % dir * maxDist
    
    ! Identify where the point is first.
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToNextFace(distance, coords)
    @assertEqual(0.02_defReal, distance, 0.02_defReal * TOL)
    
    coords % r = [-0.65_defReal, 0.33_defReal, -0.47_defReal]
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToNextFace(distance, coords)
    @assertEqual(0.385_defReal, distance, 0.385_defReal * TOL)    
    
    ! Few points outside the mesh but entering.
    coords % r = [-1.13_defReal, -0.8_defReal, 0.3_defReal]
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertTrue(coords % elementIdx > 0)
    @assertEqual(0.13_defReal, distance, 0.13_defReal * TOL)
    
    coords % r = [0.65_defReal, 0.0_defReal, 1.25_defReal]
    coords % dir = [ZERO, ZERO, -ONE]
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertTrue(coords % elementIdx > 0)
    @assertEqual(0.25_defReal, distance, 0.25_defReal * TOL)
    
    ! Few points outside the mesh and not entering.
    coords % r = [-1.13_defReal, -0.8_defReal, 0.3_defReal]
    coords % dir = [-ONE, ZERO, ZERO]
    coords % rEnd = coords % r + coords % dir * maxDist
    coords % elementIdx = 0
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertTrue(coords % elementIdx == 0)
    @assertEqual(INF, distance)
    
    coords % r = [0.65_defReal, 0.0_defReal, 1.25_defReal]
    coords % dir = [ZERO, ZERO, ONE]
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertTrue(coords % elementIdx == 0)
    @assertEqual(INF, distance)    
  
  end subroutine test_distance

end module triOpenFOAMMesh_iTest