module centroidTriangulationMethod_iTest

  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use funit
  use numPrecision
  use OpenFOAMMesh_class, only : OpenFOAMMesh
  use publicObjects,      only : coordData, newCoordData
  use universalVariables
  
  implicit none
  
  ! Parameters.
  character(*), parameter :: MESH_DEF = &
  " id 2; type OpenFOAMMesh; path ./IntegrationTestFiles/Geometry/Meshes/OpenFOAM/testMesh/; triangulationMethod centroidBased;"
  ! Variables.
  type(OpenFOAMMesh) :: mesh

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
    @assertEqual(22, mesh % getVerticesNumber())
    ! Test number of faces.
    @assertEqual(112, mesh % getFacesNumber())
    ! Test number of internal faces.
    @assertEqual(80, mesh % getInternalFacesNumber())
    ! Test number of edges.
    @assertEqual(85, mesh % getEdgesNumber())
    ! Test number of tetrahedra.
    @assertEqual(48, mesh % getElementsNumber(.true.))

  end subroutine test_info
  
  !!
  !! Test inside / outside determination.
  !!
@Test
  subroutine test_inside()
    real(defReal), dimension(3) :: r, u
    type(coordData)             :: data
    
    ! Few points inside.
    r = [0.31_defReal, 0.42_defReal, 0.13_defReal]
    u = [ONE, ZERO, ZERO]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(40, data % elementIdx)

    r = [0.02_defReal, 0.97_defReal, -0.5_defReal]
    u = [ZERO, ONE, ZERO]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(39, data % elementIdx)

    ! Few points outside.
    r = [1.2_defReal, 0.8_defReal, 0.0_defReal]
    u = [ZERO, ONE, ZERO]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)
    
    r = [0.1_defReal, 0.1_defReal, 1.13_defReal]
    u = [ZERO, ZERO, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    ! Few more difficult points.
    
    ! A point on a face.
    r = [-1.0_defReal, 0.1_defReal, 0.1_defReal]

    ! Points into the mesh.
    u = [ONE, ZERO, ZERO]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(19, data % elementIdx)

    ! Points away from the mesh.
    u = [-ONE, ZERO, ZERO]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    ! A point on an internal edge. Different directions.
    r = ZERO
    u = [2.0_defReal, ONE, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(30, data % elementIdx)

    u = [ONE, 2.0_defReal, ZERO]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(40, data % elementIdx)

    u = [-2.0_defReal, ONE, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(28, data % elementIdx)

    u = [-ONE, 2.0_defReal, ZERO]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(18, data % elementIdx)
    
    u = [2.0_defReal, -ONE, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(52, data % elementIdx)

    u = [ONE, -2.0_defReal, ZERO]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(49, data % elementIdx)
    
    u = [-2.0_defReal, -ONE, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(6, data % elementIdx)

    u = [-ONE, -2.0_defReal, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(7, data % elementIdx)

    ! A point on an internal vertex. Different directions.
    r = [0.5_defReal, -0.5_defReal, ZERO]
    u = [2.0_defReal, ONE, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(41, data % elementIdx)

    u = [ONE, 2.0_defReal, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(51, data % elementIdx)

    u = [-2.0_defReal, ONE, ZERO]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(49, data % elementIdx)

    u = [-ONE, 2.0_defReal, ZERO]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(52, data % elementIdx)

    u = [-2.0_defReal, -ONE, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(50, data % elementIdx)

    u = [-ONE, -2.0_defReal, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(44, data % elementIdx)

    u = [2.0_defReal, -ONE, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(42, data % elementIdx)

    u = [ONE, -2.0_defReal, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(43, data % elementIdx)

    u = [-ONE, ONE, -3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(45, data % elementIdx)

    u = [ONE, -ONE, -3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(46, data % elementIdx)

    u = [-ONE, ONE, 3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(48, data % elementIdx)

    u = [ONE, -ONE, 3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(47, data % elementIdx)

    ! A point on a boundary edge. Different directions.
    r = [ONE, ZERO, 0.5_defReal]

    ! Points inside the mesh.
    u = [-2.0_defReal, ONE, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(29, data % elementIdx)

    u = [-ONE, 2.0_defReal, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(32, data % elementIdx)

    u = [-2.0_defReal, -ONE, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(51, data % elementIdx)

    u = [-ONE, -2.0_defReal, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(41, data % elementIdx)

    ! Points away from the mesh.
    u = [2.0_defReal, ONE, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [ONE, 2.0_defReal, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [2.0_defReal, -ONE, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [ONE, -2.0_defReal, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    ! A point on a boundary vertex. Different directions.
    r = [ZERO, ZERO, ONE]

    ! Points into the mesh.
    u = [-3.0_defReal, -ONE, -3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(5, data % elementIdx)

    u = [-2.0_defReal, -ONE, -5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(6, data % elementIdx)

    u = [-ONE, -2.0_defReal, -5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(7, data % elementIdx)

    u = [-ONE, -3.0_defReal, -3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(8, data % elementIdx)

    u = [-ONE, -2.0_defReal, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(15, data % elementIdx)

    u = [-2.0_defReal, -ONE, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(16, data % elementIdx)

    u = [-ONE, 3.0_defReal, -3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(18, data % elementIdx)

    u = [-ONE, 2.0_defReal, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(25, data % elementIdx)

    u = [-3.0_defReal, ONE, -3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(27, data % elementIdx)

    u = [-2.0_defReal, ONE, -5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(28, data % elementIdx)

    u = [2.0_defReal, ONE, -5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(30, data % elementIdx)

    u = [2.0_defReal, ONE, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(37, data % elementIdx)

    u = [ONE, 2.0_defReal, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(38, data % elementIdx)

    u = [ONE, 2.0_defReal, -5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(40, data % elementIdx)

    u = [ONE, -2.0_defReal, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(48, data % elementIdx)

    u = [ONE, -2.0_defReal, -5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(49, data % elementIdx)

    u = [ONE, -3.0_defReal, -3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(50, data % elementIdx)

    u = [3.0_defReal, -ONE, -3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(52, data % elementIdx)

    ! Points away from the mesh.
    u = [-3.0_defReal, -ONE, 3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [-2.0_defReal, -ONE, 5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [-ONE, -2.0_defReal, 5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [-ONE, -3.0_defReal, 3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [-ONE, -2.0_defReal, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [-2.0_defReal, -ONE, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [-ONE, 3.0_defReal, 3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [-ONE, 2.0_defReal, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [-3.0_defReal, ONE, 3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [-2.0_defReal, ONE, 5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [2.0_defReal, ONE, 5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [2.0_defReal, ONE, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [ONE, 2.0_defReal, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [ONE, 2.0_defReal, 5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [ONE, -2.0_defReal, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [ONE, -2.0_defReal, 5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [ONE, -3.0_defReal, 3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    u = [3.0_defReal, -ONE, 3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

  end subroutine test_inside
  
  !!
  !! Test distance calculations.
  !!
@Test
  subroutine test_distance()
    real(defReal), dimension(3) :: r, u
    type(coordData)             :: data
    real(defReal), parameter    :: dMax = TWO, TOL = 1.0E-6

    ! Few points inside mesh.
    data = newCoordData([0.98_defReal, 0.1_defReal, 0.1_defReal], [ONE, ZERO, ZERO], dMax = dMax)
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.02_defReal, data % d, 0.02_defReal * TOL)

    data = newCoordData([-0.65_defReal, 0.33_defReal, -0.47_defReal], [ONE, ZERO, ZERO], dMax = dMax)
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.385_defReal, data % d, 0.385_defReal * TOL)    
    
    ! Few points outside the mesh but entering.
    data = newCoordData([-1.13_defReal, -0.8_defReal, 0.3_defReal], [ONE, ZERO, ZERO], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(9, data % elementIdx)
    @assertEqual(0.13_defReal, data % d, 0.13_defReal * TOL)
    
    data = newCoordData([0.65_defReal, 0.1_defReal, 1.25_defReal], [ZERO, ZERO, -ONE], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(37, data % elementIdx)
    @assertEqual(0.25_defReal, data % d, 0.25_defReal * TOL)

    ! Few more difficult points entering.
    
    ! Entering through boundary edges.
    data = newCoordData([-0.6_defReal, 0.2_defReal, -1.25_defReal], [6.0_defReal, -4.0_defReal, 5.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(13, data % elementIdx)
    @assertEqual(sqrt(77.0_defReal) / 20.0_defReal, data % d, sqrt(77.0_defReal) / 20.0_defReal * TOL)    

    data = newCoordData([-0.6_defReal, -0.2_defReal, -1.25_defReal], [6.0_defReal, 4.0_defReal, 5.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(24, data % elementIdx)
    @assertEqual(sqrt(77.0_defReal) / 20.0_defReal, data % d, sqrt(77.0_defReal) / 20.0_defReal * TOL)

    data = newCoordData([1.1_defReal, -0.1_defReal, 1.1_defReal], [-5.0_defReal, 2.0_defReal, -2.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(37, data % elementIdx)
    @assertEqual(sqrt(33.0_defReal) / 20.0_defReal, data % d, sqrt(33.0_defReal) / 20.0_defReal * TOL)

    data = newCoordData([1.1_defReal, 0.1_defReal, 1.1_defReal], [-5.0_defReal, -2.0_defReal, -2.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(48, data % elementIdx)
    @assertEqual(sqrt(33.0_defReal) / 20.0_defReal, data % d, sqrt(33.0_defReal) / 20.0_defReal * TOL)

    ! Entering through boundary vertices.
    data = newCoordData([0.2_defReal, 0.15_defReal, -1.2_defReal], [-4.0_defReal, -3.0_defReal, 4.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(13, data % elementIdx)
    @assertEqual(sqrt(41.0_defReal) / 20.0_defReal, data % d, sqrt(41.0_defReal) * TOL / 20.0_defReal)

    data = newCoordData([0.35_defReal, -0.35_defReal, -1.35_defReal], [-ONE, ONE, ONE], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(24, data % elementIdx)
    @assertEqual(7.0_defReal * sqrt(3.0_defReal) / 20.0_defReal, data % d, 7.0_defReal * sqrt(3.0_defReal) * TOL / 20.0_defReal)

    data = newCoordData([-0.1_defReal, -0.25_defReal, 1.13_defReal], [10.0_defReal, 25.0_defReal, -13.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(38, data % elementIdx)
    @assertEqual(sqrt(894.0_defReal) / 100.0_defReal, data % d, sqrt(894.0_defReal) * TOL / 100.0_defReal)

    data = newCoordData([-0.33_defReal, 0.56_defReal, 1.27_defReal], [33.0_defReal, -56.0_defReal, -27.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(48, data % elementIdx)
    @assertEqual(sqrt(4954.0_defReal) / 100.0_defReal, data % d, sqrt(4954.0_defReal) * TOL / 100.0_defReal)
    
    ! Few points outside the mesh and not entering.
    data = newCoordData([-1.13_defReal, -0.8_defReal, 0.3_defReal], [-ONE, ZERO, ZERO], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)
    
    data = newCoordData([0.65_defReal, 0.0_defReal, 1.25_defReal], [ZERO, ZERO, ONE], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)

    ! Few more difficult points still not entering.
    ! On a face.
    data = newCoordData([ONE, 0.1_defReal, 0.3_defReal], [ONE, ONE, ONE], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)

    ! Going through boundary edges but still pointing outside after intersecting. Use dirty values.
    data = newCoordData([0.9_defReal, 0.43_defReal, 1.1_defReal], [ONE, ZERO, -ONE], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)

    data = newCoordData([0.55_defReal, 0.67_defReal, 1.23_defReal], [-22.0_defReal, 33.0_defReal, -23.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)

    data = newCoordData([-0.43_defReal, -0.78_defReal, -1.05_defReal], [-57.0_defReal, 3.0_defReal, 5.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)

    data = newCoordData([-0.68_defReal, -0.98_defReal, -1.76_defReal], [ZERO, -ONE, 38.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)

    ! Going through boundary vertices but still pointing outside after intersecting. Use dirty values.
    data = newCoordData([1.42_defReal, 0.77_defReal, 0.87_defReal], [-42.0_defReal, 23.0_defReal, 13.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)

    data = newCoordData([-0.54_defReal, 0.99_defReal, 1.02_defReal], [-46.0_defReal, ONE, -2.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)

    data = newCoordData([0.98_defReal, -0.99_defReal, 2.05_defReal], [2.0_defReal, -ONE, -105.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)

    data = newCoordData([-1.1_defReal, -0.9_defReal, 1.1_defReal], [ONE, -ONE, -ONE], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)

    data = newCoordData([1.31_defReal, 0.54_defReal, -1.45_defReal], [-31.0_defReal, 46.0_defReal, 45.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)

    data = newCoordData([-0.21_defReal, 0.18_defReal, -1.01_defReal], [-79.0_defReal, 82.0_defReal, ONE], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)

    data = newCoordData([1.98_defReal, -0.98_defReal, -1.66_defReal], [-49.0_defReal, -ONE, 33.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)

    data = newCoordData([-0.99_defReal, -0.99_defReal, -2.67_defReal], [-ONE, -ONE, 167.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)
  
  end subroutine test_distance

end module centroidTriangulationMethod_iTest