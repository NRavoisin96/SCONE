module OpenFOAMMesh_iTest
  
  use coord_class,        only : coord
  use OpenFOAMMesh_class, only : OpenFOAMMesh
  use numPrecision
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use funit
  use universalVariables
  
  implicit none
  
  ! Parameters.
  character(*), parameter :: MESH_DEF = &
  " id 2; type OpenFOAMMesh; path ./IntegrationTestFiles/Geometry/Meshes/OpenFOAM/testMesh/;"
  ! Variables.
  type(OpenFOAMMesh)      :: mesh
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
    
    u = [ONE, ONE, ONE]
    u = u / norm2(u)

    ! Few points inside.
    r = [-0.32_defReal, -0.65_defReal, 0.73_defReal]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(1, elementIdx)

    r = [-0.02_defReal, 0.34_defReal, -0.56_defReal]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(2, elementIdx)

    r = [0.31_defReal, 0.42_defReal, 0.13_defReal]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(3, elementIdx)

    r = [0.89_defReal, -0.93_defReal, -0.21_defReal]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(4, elementIdx)
    
    ! Few points outside.
    r = [1.2_defReal, 0.8_defReal, 0.0_defReal]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)
    
    r = [0.1_defReal, 0.1_defReal, 1.13_defReal]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)
    
    ! Few more difficult points.
    
    ! A point on a face.
    r = [-1.0_defReal, 0.1_defReal, 0.1_defReal]
    
    ! Points into the mesh.
    u = [ONE, ZERO, ZERO]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(2, elementIdx)
    
    ! Points away from the mesh.
    u = [-ONE, ZERO, ZERO]
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)
    
    ! A point on an internal edge. Different directions.
    r = ZERO
    u = [-ONE, -ONE, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(1, elementIdx)

    u = [-ONE, ONE, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(2, elementIdx)

    u = [ONE, ONE, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(3, elementIdx)
    
    u = [ONE, -ONE, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(4, elementIdx)

    ! A point on a boundary edge. Different directions.
    r = [ONE, ZERO, 0.5_defReal]

    ! Points inside the mesh.
    u = [-ONE, ONE, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(3, elementIdx)

    u = [-ONE, -ONE, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(4, elementIdx)

    ! Points outside the mesh.
    u = [ONE, ONE, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [ONE, -ONE, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    ! A point on a boundary vertex. Different directions.
    r = [ZERO, ZERO, ONE]

    ! Points inside the mesh.
    u = [-ONE, -ONE, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(1, elementIdx)

    u = [-ONE, ONE, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(2, elementIdx)

    u = [ONE, ONE, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(3, elementIdx)

    u = [ONE, -ONE, -ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(4, elementIdx)

    ! Points outside the mesh.
    u = [-ONE, -ONE, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [-ONE, ONE, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [ONE, ONE, ONE]
    u = u / norm2(u)
    call mesh % findElementAndParentIdxs(r, u, elementIdx, parentIdx)
    @assertEqual(0, elementIdx)

    u = [ONE, -ONE, ONE]
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
    integer(shortInt)           :: parentIdx
    real(defReal), parameter    :: TOL = 1.0E-6, maxDist = 2.0_defReal
    
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
    @assertEqual(0.65_defReal, distance, 0.65_defReal * TOL)    
    
    ! Few points outside the mesh but entering.
    coords % r = [-1.13_defReal, -0.8_defReal, 0.3_defReal]
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(1, coords % elementIdx)
    @assertEqual(0.13_defReal, distance, 0.13_defReal * TOL)
    
    coords % r = [0.65_defReal, 0.1_defReal, 1.25_defReal]
    coords % dir = [ZERO, ZERO, -ONE]
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(3, coords % elementIdx)
    @assertEqual(0.25_defReal, distance, 0.25_defReal * TOL)

    ! Few more difficult points entering.
    
    ! Entering through boundary edges.
    coords % r = [-0.6_defReal, 0.2_defReal, -1.25_defReal]
    coords % dir = [6.0_defReal, -4.0_defReal, 5.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(1, coords % elementIdx)
    @assertEqual(sqrt(77.0_defReal) / 20.0_defReal, distance, sqrt(77.0_defReal) / 20.0_defReal * TOL)    

    coords % r = [-0.6_defReal, -0.2_defReal, -1.25_defReal]
    coords % dir = [6.0_defReal, 4.0_defReal, 5.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(2, coords % elementIdx)
    @assertEqual(sqrt(77.0_defReal) / 20.0_defReal, distance, sqrt(77.0_defReal) / 20.0_defReal * TOL)

    coords % r = [1.1_defReal, -0.1_defReal, 1.1_defReal]
    coords % dir = [-5.0_defReal, 2.0_defReal, -2.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(3, coords % elementIdx)
    @assertEqual(sqrt(33.0_defReal) / 20.0_defReal, distance, sqrt(33.0_defReal) / 20.0_defReal * TOL)

    coords % r = [1.1_defReal, 0.1_defReal, 1.1_defReal]
    coords % dir = [-5.0_defReal, -2.0_defReal, -2.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(4, coords % elementIdx)
    @assertEqual(sqrt(33.0_defReal) / 20.0_defReal, distance, sqrt(33.0_defReal) / 20.0_defReal * TOL)

    ! Entering through boundary vertices.
    coords % r = [0.2_defReal, 0.2_defReal, -1.2_defReal]
    coords % dir = [-ONE, -ONE, ONE]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(1, coords % elementIdx)
    @assertEqual(sqrt(3.0_defReal) / 5.0_defReal, distance, sqrt(3.0_defReal) * TOL / 5.0_defReal)

    coords % r = [0.35_defReal, -0.35_defReal, -1.35_defReal]
    coords % dir = [-ONE, ONE, ONE]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(2, coords % elementIdx)
    @assertEqual(7.0_defReal * sqrt(3.0_defReal) / 20.0_defReal, distance, 7.0_defReal * sqrt(3.0_defReal) * TOL / 20.0_defReal)

    coords % r = [-0.1_defReal, -0.25_defReal, 1.13_defReal]
    coords % dir = [10.0_defReal, 25.0_defReal, -13.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(3, coords % elementIdx)
    @assertEqual(sqrt(894.0_defReal) / 100.0_defReal, distance, sqrt(894.0_defReal) * TOL / 100.0_defReal)

    coords % r = [-0.33_defReal, 0.56_defReal, 1.27_defReal]
    coords % dir = [33.0_defReal, -56.0_defReal, -27.0_defReal]
    coords % dir = coords % dir / norm2(coords % dir)
    coords % rEnd = coords % r + coords % dir * maxDist
    call mesh % findElementAndParentIdxs(coords % r, coords % dir, coords % elementIdx, parentIdx)
    call mesh % distanceToBoundary(distance, coords, parentIdx)
    @assertEqual(4, coords % elementIdx)
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
    coords % dir = [ONE, -ONE, ONE]
    coords % dir = coords % dir / norm2(coords % dir)
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

end module OpenFOAMMesh_iTest