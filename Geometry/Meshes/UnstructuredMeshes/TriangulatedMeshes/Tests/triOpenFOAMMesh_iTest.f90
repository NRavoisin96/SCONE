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
    real(defReal)     :: TOL = 1.0E-6
    integer(shortInt) :: i
    
    ! Test number of vertices.
    @assertEqual(22, mesh % nVertices)
    ! Test number of faces.
    @assertEqual(20, mesh % nFaces)
    ! Test number of elements.
    @assertEqual(4, mesh % nElements)
    ! Test number of internal faces.
    @assertEqual(4, mesh % nInternalFaces)
    ! Test number of edges.
    @assertEqual(85, mesh % nEdges)
    ! Test number of tetrahedra.
    @assertEqual(48, mesh % nTetrahedra)
    ! Test every tetrahedron.
    do i = 1, mesh % nTetrahedra
      ! Test that none of the triangle indices and edge indices are zero.
      @assertTrue(all(mesh % tetrahedra % getTetrahedronTriangleIdxs(i) /= 0))
      @assertTrue(all(mesh % tetrahedra % getTetrahedronEdgeIdxs(i) > 0))

    end do
    ! Test number of triangles.
    @assertEqual(112, mesh % nTriangles)
    ! Test every triangle.
    do i = 1, mesh % nTriangles
      ! Test that none of the triangle vertex indices and edge indices are zero.
      @assertTrue(all(mesh % triangles % getTriangleVertexIdxs(i) > 0))
      @assertTrue(all(mesh % triangles % getTriangleEdgeIdxs(i) > 0))

    end do
    
    ! Test area of two faces.
    @assertEqual(2.0_defReal, mesh % faces % getFaceArea(1), 2.0_defReal * TOL)
    @assertEqual(1.0_defReal, mesh % faces % getFaceArea(11), 1.0_defReal * TOL)
    ! Test volume of one element.
    @assertEqual(2.0_defReal, mesh % elements % getElementVolume(3), 2.0_defReal * TOL)

  end subroutine test_info
  
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
    
    ! A point on a face.
    r = [-1.0_defReal, 0.1_defReal, 0.1_defReal]

    ! Points into the mesh.
    u = [ONE, ZERO, ZERO]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(15, elementIdx)

    ! Points away from the mesh.
    u = [-ONE, ZERO, ZERO]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)
    
    ! A point on an internal edge. Different directions.
    r = ZERO
    u = [2.0_defReal, ONE, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(26, elementIdx)

    u = [ONE, 2.0_defReal, ZERO]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(36, elementIdx)

    u = [-2.0_defReal, ONE, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(24, elementIdx)

    u = [-ONE, 2.0_defReal, ZERO]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(14, elementIdx)
    
    u = [2.0_defReal, -ONE, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(48, elementIdx)

    u = [ONE, -2.0_defReal, ZERO]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(45, elementIdx)
    
    u = [-2.0_defReal, -ONE, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(2, elementIdx)

    u = [-ONE, -2.0_defReal, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(3, elementIdx)

    ! A point on an internal vertex. Different directions.
    r = [0.5_defReal, -0.5_defReal, ZERO]
    u = [2.0_defReal, ONE, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(37, elementIdx)

    u = [ONE, 2.0_defReal, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(47, elementIdx)

    u = [-2.0_defReal, ONE, ZERO]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(45, elementIdx)

    u = [-ONE, 2.0_defReal, ZERO]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(48, elementIdx)

    u = [-2.0_defReal, -ONE, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(46, elementIdx)

    u = [-ONE, -2.0_defReal, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(40, elementIdx)

    u = [2.0_defReal, -ONE, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(38, elementIdx)

    u = [ONE, -2.0_defReal, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(39, elementIdx)

    u = [-ONE, ONE, -3.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(41, elementIdx)

    u = [ONE, -ONE, -3.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(42, elementIdx)

    u = [-ONE, ONE, 3.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(44, elementIdx)

    u = [ONE, -ONE, 3.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(43, elementIdx)

    ! A point on a boundary edge. Different directions.
    r = [ONE, ZERO, 0.5_defReal]

    ! Points inside the mesh.
    u = [-2.0_defReal, ONE, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(25, elementIdx)

    u = [-ONE, 2.0_defReal, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(28, elementIdx)

    u = [-2.0_defReal, -ONE, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(47, elementIdx)

    u = [-ONE, -2.0_defReal, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(37, elementIdx)

    ! Points away from the mesh.
    u = [2.0_defReal, ONE, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [ONE, 2.0_defReal, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [2.0_defReal, -ONE, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [ONE, -2.0_defReal, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    ! A point on a boundary vertex. Different directions.
    r = [ZERO, ZERO, ONE]

    ! Points into the mesh.
    u = [-3.0_defReal, -ONE, -3.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(1, elementIdx)

    u = [-2.0_defReal, -ONE, -5.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(2, elementIdx)

    u = [-ONE, -2.0_defReal, -5.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(3, elementIdx)

    u = [-ONE, -3.0_defReal, -3.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(4, elementIdx)

    u = [-ONE, -2.0_defReal, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(11, elementIdx)

    u = [-2.0_defReal, -ONE, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(12, elementIdx)

    u = [-ONE, 3.0_defReal, -3.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(14, elementIdx)

    u = [-ONE, 2.0_defReal, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(21, elementIdx)

    u = [-3.0_defReal, ONE, -3.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(23, elementIdx)

    u = [-2.0_defReal, ONE, -5.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(24, elementIdx)

    u = [2.0_defReal, ONE, -5.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(26, elementIdx)

    u = [2.0_defReal, ONE, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(33, elementIdx)

    u = [ONE, 2.0_defReal, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(34, elementIdx)

    u = [ONE, 2.0_defReal, -5.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(36, elementIdx)

    u = [ONE, -2.0_defReal, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(44, elementIdx)

    u = [ONE, -2.0_defReal, -5.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(45, elementIdx)

    u = [ONE, -3.0_defReal, -3.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(46, elementIdx)

    u = [3.0_defReal, -ONE, -3.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(48, elementIdx)

    ! Points away from the mesh.
    u = [-3.0_defReal, -ONE, 3.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [-2.0_defReal, -ONE, 5.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [-ONE, -2.0_defReal, 5.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [-ONE, -3.0_defReal, 3.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [-ONE, -2.0_defReal, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [-2.0_defReal, -ONE, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [-ONE, 3.0_defReal, 3.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [-ONE, 2.0_defReal, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [-3.0_defReal, ONE, 3.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [-2.0_defReal, ONE, 5.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [2.0_defReal, ONE, 5.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [2.0_defReal, ONE, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [ONE, 2.0_defReal, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [ONE, 2.0_defReal, 5.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [ONE, -2.0_defReal, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [ONE, -2.0_defReal, 5.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [ONE, -3.0_defReal, 3.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [3.0_defReal, -ONE, 3.0_defReal]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

  end subroutine test_inside
  
  !!
  !! Test distance calculations.
  !!
@Test
  subroutine test_distance()
    real(defReal)               :: distance
    integer(shortInt)           :: parentIdx, i
    real(defReal), parameter    :: TOL = 1.0E-6, maxDist = 2.0_defReal
    integer(shortInt), dimension(:), allocatable :: vertexIdxs
    
    ! Few points inside mesh.
    coords % r = [0.98_defReal, 0.1_defReal, 0.1_defReal]
    coords % dir = [ONE, ZERO, ZERO]
    coords % rEnd = coords % r + coords % dir * maxDist
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
    @assertEqual(5, coords % elementIdx)
    @assertEqual(0.13_defReal, distance, 0.13_defReal * TOL)
    
    coords % r = [0.65_defReal, 0.1_defReal, 1.25_defReal]
    coords % dir = [ZERO, ZERO, -ONE]
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(33, coords % elementIdx)
    @assertEqual(0.25_defReal, distance, 0.25_defReal * TOL)

    ! Few more difficult points entering.
    
    ! Entering through boundary edges.
    coords % r = [-0.6_defReal, 0.2_defReal, -1.25_defReal]
    coords % dir = [6.0_defReal, -4.0_defReal, 5.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(9, coords % elementIdx)
    @assertEqual(sqrt(77.0_defReal) / 20.0_defReal, distance, sqrt(77.0_defReal) / 20.0_defReal * TOL)    

    coords % r = [-0.6_defReal, -0.2_defReal, -1.25_defReal]
    coords % dir = [6.0_defReal, 4.0_defReal, 5.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(20, coords % elementIdx)
    @assertEqual(sqrt(77.0_defReal) / 20.0_defReal, distance, sqrt(77.0_defReal) / 20.0_defReal * TOL)

    coords % r = [1.1_defReal, -0.1_defReal, 1.1_defReal]
    coords % dir = [-5.0_defReal, 2.0_defReal, -2.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(33, coords % elementIdx)
    @assertEqual(sqrt(33.0_defReal) / 20.0_defReal, distance, sqrt(33.0_defReal) / 20.0_defReal * TOL)

    coords % r = [1.1_defReal, 0.1_defReal, 1.1_defReal]
    coords % dir = [-5.0_defReal, -2.0_defReal, -2.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(44, coords % elementIdx)
    @assertEqual(sqrt(33.0_defReal) / 20.0_defReal, distance, sqrt(33.0_defReal) / 20.0_defReal * TOL)

    ! Entering through boundary vertices.
    coords % r = [0.2_defReal, 0.15_defReal, -1.2_defReal]
    coords % dir = [-4.0_defReal, -3.0_defReal, 4.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(9, coords % elementIdx)
    @assertEqual(sqrt(41.0_defReal) / 20.0_defReal, distance, sqrt(41.0_defReal) * TOL / 20.0_defReal)

    coords % r = [0.35_defReal, -0.35_defReal, -1.35_defReal]
    coords % dir = [-ONE, ONE, ONE]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(20, coords % elementIdx)
    @assertEqual(7.0_defReal * sqrt(3.0_defReal) / 20.0_defReal, distance, 7.0_defReal * sqrt(3.0_defReal) * TOL / 20.0_defReal)

    coords % r = [-0.1_defReal, -0.25_defReal, 1.13_defReal]
    coords % dir = [10.0_defReal, 25.0_defReal, -13.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(34, coords % elementIdx)
    @assertEqual(sqrt(894.0_defReal) / 100.0_defReal, distance, sqrt(894.0_defReal) * TOL / 100.0_defReal)

    coords % r = [-0.33_defReal, 0.56_defReal, 1.27_defReal]
    coords % dir = [33.0_defReal, -56.0_defReal, -27.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(44, coords % elementIdx)
    @assertEqual(sqrt(4954.0_defReal) / 100.0_defReal, distance, sqrt(4954.0_defReal) * TOL / 100.0_defReal)
    
    ! Few points outside the mesh and not entering.
    coords % r = [-1.13_defReal, -0.8_defReal, 0.3_defReal]
    coords % dir = [-ONE, ZERO, ZERO]
    coords % rEnd = coords % r + coords % dir * maxDist
    coords % elementIdx = 0
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(0, coords % elementIdx)
    @assertEqual(INF, distance)
    
    coords % r = [0.65_defReal, 0.0_defReal, 1.25_defReal]
    coords % dir = [ZERO, ZERO, ONE]
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(0, coords % elementIdx)
    @assertEqual(INF, distance)

    ! Few more difficult points still not entering.
    ! On a face.
    coords % r = [ONE, 0.1_defReal, 0.3_defReal]
    coords % dir = [ONE, ONE, ONE]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(0, coords % elementIdx)
    @assertEqual(INF, distance)

    ! Going through boundary edges but still pointing outside after intersecting. Use dirty values.
    coords % r = [0.9_defReal, 0.43_defReal, 1.1_defReal]
    coords % dir = [ONE, ZERO, -ONE]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(0, coords % elementIdx)
    @assertEqual(INF, distance)

    coords % r = [0.55_defReal, 0.67_defReal, 1.23_defReal]
    coords % dir = [-22.0_defReal, 33.0_defReal, -23.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(0, coords % elementIdx)
    @assertEqual(INF, distance)

    coords % r = [-0.43_defReal, -0.78_defReal, -1.05_defReal]
    coords % dir = [-57.0_defReal, 3.0_defReal, 5.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(0, coords % elementIdx)
    @assertEqual(INF, distance)

    coords % r = [-0.68_defReal, -0.98_defReal, -1.76_defReal]
    coords % dir = [ZERO, -ONE, 38.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(0, coords % elementIdx)
    @assertEqual(INF, distance)

    ! Going through boundary vertices but still pointing outside after intersecting. Use dirty values.
    coords % r = [1.42_defReal, 0.77_defReal, 0.87_defReal]
    coords % dir = [-42.0_defReal, 23.0_defReal, 13.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(0, coords % elementIdx)
    @assertEqual(INF, distance)

    coords % r = [-0.54_defReal, 0.99_defReal, 1.02_defReal]
    coords % dir = [-46.0_defReal, ONE, -2.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(0, coords % elementIdx)
    @assertEqual(INF, distance)

    coords % r = [0.98_defReal, -0.99_defReal, 2.05_defReal]
    coords % dir = [2.0_defReal, -ONE, -105.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(0, coords % elementIdx)
    @assertEqual(INF, distance)

    coords % r = [-1.1_defReal, -0.9_defReal, 1.1_defReal]
    coords % dir = [ONE, -ONE, -ONE]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(0, coords % elementIdx)
    @assertEqual(INF, distance)

    coords % r = [1.31_defReal, 0.54_defReal, -1.45_defReal]
    coords % dir = [-31.0_defReal, 46.0_defReal, 45.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(0, coords % elementIdx)
    @assertEqual(INF, distance)

    coords % r = [-0.21_defReal, 0.18_defReal, -1.01_defReal]
    coords % dir = [-79.0_defReal, 82.0_defReal, ONE]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(0, coords % elementIdx)
    @assertEqual(INF, distance)

    coords % r = [1.98_defReal, -0.98_defReal, -1.66_defReal]
    coords % dir = [-49.0_defReal, -ONE, 33.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(0, coords % elementIdx)
    @assertEqual(INF, distance)

    coords % r = [-0.99_defReal, -0.99_defReal, -2.67_defReal]
    coords % dir = [-ONE, -ONE, 167.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(0, coords % elementIdx)
    @assertEqual(INF, distance)
  
  end subroutine test_distance

end module triOpenFOAMMesh_iTest