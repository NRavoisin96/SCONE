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
    type(coord)                 :: coords
    real(defReal), dimension(3) :: r, u
    
    u = [ONE, ONE, ONE]
    u = u / norm2(u)
    call coords % setDirection(u)

    ! Few points inside.
    call coords % setPosition([-0.32_defReal, -0.65_defReal, 0.73_defReal])
    call mesh % findHostElement(coords)
    @assertEqual(1, coords % getElementIdx())

    call coords % setPosition([-0.02_defReal, 0.34_defReal, -0.56_defReal])
    call mesh % findHostElement(coords)
    @assertEqual(2, coords % getElementIdx())

    call coords % setPosition([0.31_defReal, 0.42_defReal, 0.13_defReal])
    call mesh % findHostElement(coords)
    @assertEqual(3, coords % getElementIdx())

    call coords % setPosition([0.89_defReal, -0.93_defReal, -0.21_defReal])
    call mesh % findHostElement(coords)
    @assertEqual(4, coords % getElementIdx())
    
    ! Few points outside.
    call coords % setPosition([1.2_defReal, 0.8_defReal, 0.0_defReal])
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())
    
    call coords % setPosition([0.1_defReal, 0.1_defReal, 1.13_defReal])
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())
    
    ! Few more difficult points.
    
    ! A point on a face.
    r = [-ONE, 0.1_defReal, 0.1_defReal]
    
    ! Points into the mesh.
    call coords % setPosition(r)
    call coords % setDirection([ONE, ZERO, ZERO])
    call mesh % findHostElement(coords)
    @assertEqual(2, coords % getElementIdx())
    
    ! Points away from the mesh.
    call coords % setPosition(r)
    call coords % setDirection([-ONE, ZERO, ZERO])
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())
    
    ! A point on an internal edge. Different directions.
    r = [ZERO, ZERO, ZERO]
    call coords % setPosition(r)
    u = [-ONE, -ONE, -ONE]
    u = u / norm2(u)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(1, coords % getElementIdx())

    call coords % setPosition(r)
    u = [-ONE, ONE, ONE]
    u = u / norm2(u)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(2, coords % getElementIdx())

    call coords % setPosition(r)
    u = [ONE, ONE, ONE]
    u = u / norm2(u)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(3, coords % getElementIdx())
    
    call coords % setPosition(r)
    u = [ONE, -ONE, -ONE]
    u = u / norm2(u)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(4, coords % getElementIdx())

    ! A point on a boundary edge. Different directions.
    r = [ONE, ZERO, 0.5_defReal]

    ! Points inside the mesh.
    call coords % setPosition(r)
    u = [-ONE, ONE, ONE]
    u = u / norm2(u)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(3, coords % getElementIdx())

    call coords % setPosition(r)
    u = [-ONE, -ONE, ONE]
    u = u / norm2(u)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(4, coords % getElementIdx())

    ! Points outside the mesh.
    call coords % setPosition(r)
    u = [ONE, ONE, ONE]
    u = u / norm2(u)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    call coords % setPosition(r)
    u = [ONE, -ONE, -ONE]
    u = u / norm2(u)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    ! A point on a boundary vertex. Different directions.
    r = [ZERO, ZERO, ONE]

    ! Points inside the mesh.
    call coords % setPosition(r)
    u = [-ONE, -ONE, -ONE]
    u = u / norm2(u)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(1, coords % getElementIdx())

    call coords % setPosition(r)
    u = [-ONE, ONE, -ONE]
    u = u / norm2(u)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(2, coords % getElementIdx())

    call coords % setPosition(r)
    u = [ONE, ONE, -ONE]
    u = u / norm2(u)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(3, coords % getElementIdx())

    call coords % setPosition(r)
    u = [ONE, -ONE, -ONE]
    u = u / norm2(u)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(4, coords % getElementIdx())

    ! Points outside the mesh.
    call coords % setPosition(r)
    u = [-ONE, -ONE, ONE]
    u = u / norm2(u)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    call coords % setPosition(r)
    u = [-ONE, ONE, ONE]
    u = u / norm2(u)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    call coords % setPosition(r)
    u = [ONE, ONE, ONE]
    u = u / norm2(u)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    call coords % setPosition(r)
    u = [ONE, -ONE, ONE]
    u = u / norm2(u)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

  end subroutine test_inside
  
  !!
  !! Test distance calculations.
  !!
@Test
  subroutine test_distance()
    real(defReal), dimension(3) :: r, rEnd, u
    type(coord)                 :: coords
    real(defReal)               :: distance
    real(defReal), parameter    :: TOL = 1.0E-6, maxDist = 2.0_defReal
    
    ! Few points inside mesh.
    r = [0.98_defReal, 0.1_defReal, 0.1_defReal]
    u = [ONE, ZERO, ZERO]
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToNextFace(distance, coords)
    @assertEqual(0.02_defReal, distance, 0.02_defReal * TOL)
    
    r = [-0.65_defReal, 0.33_defReal, -0.47_defReal]
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToNextFace(distance, coords)
    @assertEqual(0.65_defReal, distance, 0.65_defReal * TOL)    
    
    ! Few points outside the mesh but entering.
    r = [-1.13_defReal, -0.8_defReal, 0.3_defReal]
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(1, coords % getElementIdx())
    @assertEqual(0.13_defReal, distance, 0.13_defReal * TOL)
    
    r = [0.65_defReal, 0.1_defReal, 1.25_defReal]
    u = [ZERO, ZERO, -ONE]
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(3, coords % getElementIdx())
    @assertEqual(0.25_defReal, distance, 0.25_defReal * TOL)

    ! Few more difficult points entering.
    
    ! Entering through boundary edges.
    r = [-0.6_defReal, 0.2_defReal, -1.25_defReal]
    u = [6.0_defReal, -4.0_defReal, 5.0_defReal]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(1, coords % getElementIdx())
    @assertEqual(sqrt(77.0_defReal) / 20.0_defReal, distance, sqrt(77.0_defReal) / 20.0_defReal * TOL)    

    r = [-0.6_defReal, -0.2_defReal, -1.25_defReal]
    u = [6.0_defReal, 4.0_defReal, 5.0_defReal]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(2, coords % getElementIdx())
    @assertEqual(sqrt(77.0_defReal) / 20.0_defReal, distance, sqrt(77.0_defReal) / 20.0_defReal * TOL)

    r = [1.1_defReal, -0.1_defReal, 1.1_defReal]
    u = [-5.0_defReal, 2.0_defReal, -2.0_defReal]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(3, coords % getElementIdx())
    @assertEqual(sqrt(33.0_defReal) / 20.0_defReal, distance, sqrt(33.0_defReal) / 20.0_defReal * TOL)

    r = [1.1_defReal, 0.1_defReal, 1.1_defReal]
    u = [-5.0_defReal, -2.0_defReal, -2.0_defReal]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(4, coords % getElementIdx())
    @assertEqual(sqrt(33.0_defReal) / 20.0_defReal, distance, sqrt(33.0_defReal) / 20.0_defReal * TOL)

    ! Entering through boundary vertices.
    r = [0.2_defReal, 0.2_defReal, -1.2_defReal]
    u = [-ONE, -ONE, ONE]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(1, coords % getElementIdx())
    @assertEqual(sqrt(3.0_defReal) / 5.0_defReal, distance, sqrt(3.0_defReal) * TOL / 5.0_defReal)

    r = [0.35_defReal, -0.35_defReal, -1.35_defReal]
    u = [-ONE, ONE, ONE]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(2, coords % getElementIdx())
    @assertEqual(7.0_defReal * sqrt(3.0_defReal) / 20.0_defReal, distance, 7.0_defReal * sqrt(3.0_defReal) * TOL / 20.0_defReal)

    r = [-0.1_defReal, -0.25_defReal, 1.13_defReal]
    u = [10.0_defReal, 25.0_defReal, -13.0_defReal]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(3, coords % getElementIdx())
    @assertEqual(sqrt(894.0_defReal) / 100.0_defReal, distance, sqrt(894.0_defReal) * TOL / 100.0_defReal)

    r = [-0.33_defReal, 0.56_defReal, 1.27_defReal]
    u = [33.0_defReal, -56.0_defReal, -27.0_defReal]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(4, coords % getElementIdx())
    @assertEqual(sqrt(4954.0_defReal) / 100.0_defReal, distance, sqrt(4954.0_defReal) * TOL / 100.0_defReal)
    
    ! Few points outside the mesh and not entering.
    r = [-1.13_defReal, -0.8_defReal, 0.3_defReal]
    u = [-ONE, ZERO, ZERO]
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call coords % setElementIdx(0)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(0, coords % getElementIdx())
    @assertEqual(INF, distance)
    
    r = [0.65_defReal, 0.0_defReal, 1.25_defReal]
    u = [ONE, -ONE, ONE]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(0, coords % getElementIdx())
    @assertEqual(INF, distance)

    ! Few more difficult points still not entering.
    ! On a face.
    r = [ONE, 0.1_defReal, 0.3_defReal]
    u = [ONE, ONE, ONE]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(0, coords % getElementIdx())
    @assertEqual(INF, distance)

    ! Going through boundary edges but still pointing outside after intersecting. Use dirty values.
    r = [0.9_defReal, 0.43_defReal, 1.1_defReal]
    u = [ONE, ZERO, -ONE]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(0, coords % getElementIdx())
    @assertEqual(INF, distance)

    r = [0.55_defReal, 0.67_defReal, 1.23_defReal]
    u = [-22.0_defReal, 33.0_defReal, -23.0_defReal]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(0, coords % getElementIdx())
    @assertEqual(INF, distance)

    r = [-0.43_defReal, -0.78_defReal, -1.05_defReal]
    u = [-57.0_defReal, 3.0_defReal, 5.0_defReal]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(0, coords % getElementIdx())
    @assertEqual(INF, distance)

    r = [-0.68_defReal, -0.98_defReal, -1.76_defReal]
    u = [ZERO, -ONE, 38.0_defReal]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(0, coords % getElementIdx())
    @assertEqual(INF, distance)

    ! Going through boundary vertices but still pointing outside after intersecting. Use dirty values.
    r = [1.42_defReal, 0.77_defReal, 0.87_defReal]
    u = [-42.0_defReal, 23.0_defReal, 13.0_defReal]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(0, coords % getElementIdx())
    @assertEqual(INF, distance)

    r = [-0.54_defReal, 0.99_defReal, 1.02_defReal]
    u = [-46.0_defReal, ONE, -2.0_defReal]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(0, coords % getElementIdx())
    @assertEqual(INF, distance)

    r = [0.98_defReal, -0.99_defReal, 2.05_defReal]
    u = [2.0_defReal, -ONE, -105.0_defReal]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(0, coords % getElementIdx())
    @assertEqual(INF, distance)

    r = [-1.1_defReal, -0.9_defReal, 1.1_defReal]
    u = [ONE, -ONE, -ONE]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(0, coords % getElementIdx())
    @assertEqual(INF, distance)

    r = [1.31_defReal, 0.54_defReal, -1.45_defReal]
    u = [-31.0_defReal, 46.0_defReal, 45.0_defReal]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(0, coords % getElementIdx())
    @assertEqual(INF, distance)

    r = [-0.21_defReal, 0.18_defReal, -1.01_defReal]
    u = [-79.0_defReal, 82.0_defReal, ONE]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(0, coords % getElementIdx())
    @assertEqual(INF, distance)

    r = [1.98_defReal, -0.98_defReal, -1.66_defReal]
    u = [-49.0_defReal, -ONE, 33.0_defReal]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(0, coords % getElementIdx())
    @assertEqual(INF, distance)

    r = [-0.99_defReal, -0.99_defReal, -2.67_defReal]
    u = [-ONE, -ONE, 167.0_defReal]
    u = u / norm2(u)
    rEnd = r + u * maxDist
    call coords % setPosition(r)
    call coords % setDirection(u)
    call coords % setEndPosition(rEnd)
    call mesh % findHostElement(coords)
    call mesh % distanceToBoundary(distance, coords)
    @assertEqual(0, coords % getElementIdx())
    @assertEqual(INF, distance)
  
  end subroutine test_distance

end module OpenFOAMMesh_iTest