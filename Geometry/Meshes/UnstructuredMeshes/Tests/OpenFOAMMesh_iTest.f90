module OpenFOAMMesh_iTest
  
  use charMap_class,      only : charMap
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use funit
  use numPrecision
  use OpenFOAMMesh_class, only : OpenFOAMMesh
  use publicObjects,      only : coordData, newCoordData
  use universalVariables
  use ratint
  
  implicit none
  
  ! Parameters.
  character(*), parameter :: MESH_DEF = &
  " id 2; type OpenFOAMMesh; path ./IntegrationTestFiles/Geometry/Meshes/OpenFOAM/testMesh/; fills (fuel);"

  character(*), parameter :: MESH_DEF1 = &
  " id 2; type OpenFOAMMesh; path ./IntegrationTestFiles/Geometry/Meshes/OpenFOAM/testMesh1/; fills (fuel);"
  
  ! Variables.
  type(charMap)      :: mats
  type(charMap)      :: mats1
  type(OpenFOAMMesh) :: mesh
  type(OpenFOAMMesh) :: mesh1

contains
  
  !!
  !! Import the mesh.
  !!
@Before
  subroutine setUp()
    type(dictionary)   :: dict
    type(dictionary)   :: dict1
    character(nameLen) :: name
    character(pathLen) :: path

    name = 'fuel'
    call mats % add(name, 1)
    
    call charToDict(dict, MESH_DEF)
    call dict % get(path, 'path')
    call mesh % init(trim(path), dict, mats)

    call mats1 % add(name, 1)
    
    call charToDict(dict1, MESH_DEF1)
    call dict1 % get(path, 'path')
    call mesh1 % init(trim(path), dict1, mats1)
  
  
  end subroutine setUp
  
  !!
  !! Clean after tests.
  !!
@After
  subroutine cleanUp()

    call mats % kill()
    call mats1 % kill()
    call mesh % kill()
    call mesh1 % kill()

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
    @assertEqual(18, mesh % getVerticesNumber())
    ! Test number of faces.
    @assertEqual(20, mesh % getFacesNumber())
    ! Test number of elements.
    @assertEqual(4, mesh % getElementsNumber())
    ! Test number of internal faces.
    @assertEqual(4, mesh % getInternalFacesNumber())

  end subroutine test_info
  
  !!
  !! Test inside / outside determination.
  !!
@Test
  subroutine test_inside()
    type(coordData)             :: data

    ! Few points inside.
    data = newCoordData([-0.32_defReal, -0.65_defReal, 0.73_defReal], [ONE, ONE, ONE])
    call mesh % findHostElement(data)
    @assertEqual(1, data % elementIdx)

    data = newCoordData([-0.02_defReal, 0.34_defReal, -0.56_defReal], [ONE, ONE, ONE])
    call mesh % findHostElement(data)
    @assertEqual(2, data % elementIdx)

    data = newCoordData([0.31_defReal, 0.42_defReal, 0.13_defReal], [ONE, ONE, ONE])
    call mesh % findHostElement(data)
    @assertEqual(3, data % elementIdx)

    data = newCoordData([0.89_defReal, -0.93_defReal, -0.21_defReal], [ONE, ONE, ONE])
    call mesh % findHostElement(data)
    @assertEqual(4, data % elementIdx)
    
    ! Few points outside.
    data = newCoordData([1.2_defReal, 0.8_defReal, 0.0_defReal], [ONE, ONE, ONE])
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)
    
    data = newCoordData([0.1_defReal, 0.1_defReal, 1.13_defReal], [ONE, ONE, ONE])
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)
    
    ! Few more difficult points.
    
    ! A point on a face. Points into the mesh.
    data = newCoordData([-ONE, 0.1_defReal, 0.1_defReal], [ONE, ZERO, ZERO])
    call mesh % findHostElement(data)
    @assertEqual(2, data % elementIdx)
    
    ! Points away from the mesh.
    data = newCoordData([-ONE, 0.1_defReal, 0.1_defReal], [-ONE, ZERO, ZERO])
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)
    
    ! A point on an internal edge. Different directions.
    data = newCoordData([ZERO, ZERO, ZERO], [-ONE, -ONE, -ONE])
    call mesh % findHostElement(data)
    @assertEqual(1, data % elementIdx)

    data = newCoordData([ZERO, ZERO, ZERO], [-ONE, ONE, ONE])
    call mesh % findHostElement(data)
    @assertEqual(2, data % elementIdx)

    data = newCoordData([ZERO, ZERO, ZERO], [ONE, ONE, ONE])
    call mesh % findHostElement(data)
    @assertEqual(3, data % elementIdx)
    
    data = newCoordData([ZERO, ZERO, ZERO], [ONE, -ONE, -ONE])
    call mesh % findHostElement(data)
    @assertEqual(4, data % elementIdx)

    ! A point on a boundary edge. Different directions.
    ! Points inside the mesh.
    data = newCoordData([ONE, ZERO, 0.5_defReal], [-ONE, ONE, ONE])
    call mesh % findHostElement(data)
    @assertEqual(3, data % elementIdx)

    data = newCoordData([ONE, ZERO, 0.5_defReal], [-ONE, -ONE, ONE])
    call mesh % findHostElement(data)
    @assertEqual(4, data % elementIdx)

    ! Points outside the mesh.
    data = newCoordData([ONE, ZERO, 0.5_defReal], [ONE, ONE, ONE])
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    data = newCoordData([ONE, ZERO, 0.5_defReal], [ONE, -ONE, -ONE])
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    ! A point on a boundary vertex. Different directions.
    ! Points inside the mesh.
    data = newCoordData([ZERO, ZERO, ONE], [-ONE, -ONE, -ONE])
    call mesh % findHostElement(data)
    @assertEqual(1, data % elementIdx)

    data = newCoordData([ZERO, ZERO, ONE], [-ONE, ONE, -ONE])
    call mesh % findHostElement(data)
    @assertEqual(2, data % elementIdx)

    data = newCoordData([ZERO, ZERO, ONE], [ONE, ONE, -ONE])
    call mesh % findHostElement(data)
    @assertEqual(3, data % elementIdx)

    data = newCoordData([ZERO, ZERO, ONE], [ONE, -ONE, -ONE])
    call mesh % findHostElement(data)
    @assertEqual(4, data % elementIdx)

    ! Points outside the mesh.
    data = newCoordData([ZERO, ZERO, ONE], [-ONE, -ONE, ONE])
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    data = newCoordData([ZERO, ZERO, ONE], [-ONE, ONE, ONE])
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    data = newCoordData([ZERO, ZERO, ONE], [ONE, ONE, ONE])
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

    data = newCoordData([ZERO, ZERO, ONE], [ONE, -ONE, ONE])
    call mesh % findHostElement(data)
    @assertEqual(0, data % elementIdx)

  end subroutine test_inside
  
  !!
  !! Test distance calculations.
  !!
@Test
  subroutine test_distance()
    type(coordData)          :: data
    real(defReal), parameter :: dMax = TWO, TOL = 1.0E-6
    real(defReal) :: v1,v2,v3

    ! Few points inside mesh.
    data = newCoordData([0.98_defReal, 0.1_defReal, 0.1_defReal], [ONE, ZERO, ZERO], dMax = dMax)
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.02_defReal, data % d, 0.02_defReal * TOL)

    data = newCoordData([-0.65_defReal, 0.33_defReal, -0.47_defReal], [ONE, ZERO, ZERO], dMax = dMax)
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.65_defReal, data % d, 0.65_defReal * TOL)   


    ! Few points outside the mesh but entering.
    data = newCoordData([-1.13_defReal, -0.8_defReal, 0.3_defReal], [ONE, ZERO, ZERO], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(1, data % elementIdx)
    @assertEqual(0.13_defReal, data % d, 0.13_defReal * TOL)
    
    data = newCoordData([0.65_defReal, 0.1_defReal, 1.25_defReal], [ZERO, ZERO, -ONE], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(3, data % elementIdx)
    @assertEqual(0.25_defReal, data % d, 0.25_defReal * TOL)

    ! Few more difficult points entering.
    
    ! Entering through boundary edges.
    data = newCoordData([-0.6_defReal, 0.2_defReal, -1.25_defReal], [6.0_defReal, -4.0_defReal, 5.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(1, data % elementIdx)
    @assertEqual(sqrt(77.0_defReal) / 20.0_defReal, data % d, sqrt(77.0_defReal) / 20.0_defReal * TOL)    

    data = newCoordData([-0.6_defReal, -0.2_defReal, -1.25_defReal], [6.0_defReal, 4.0_defReal, 5.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(2, data % elementIdx)
    @assertEqual(sqrt(77.0_defReal) / 20.0_defReal, data % d, sqrt(77.0_defReal) / 20.0_defReal * TOL)

    data = newCoordData([1.1_defReal, -0.1_defReal, 1.1_defReal], [-5.0_defReal, 2.0_defReal, -2.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(3, data % elementIdx)
    @assertEqual(sqrt(33.0_defReal) / 20.0_defReal, data % d, sqrt(33.0_defReal) / 20.0_defReal * TOL)

    data = newCoordData([1.1_defReal, 0.1_defReal, 1.1_defReal], [-5.0_defReal, -2.0_defReal, -2.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(4, data % elementIdx)
    @assertEqual(sqrt(33.0_defReal) / 20.0_defReal, data % d, sqrt(33.0_defReal) / 20.0_defReal * TOL)

    ! Entering through boundary vertices.
    data = newCoordData([0.2_defReal, 0.2_defReal, -1.2_defReal], [-ONE, -ONE, ONE], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(1, data % elementIdx)
    @assertEqual(sqrt(3.0_defReal) / 5.0_defReal, data % d, sqrt(3.0_defReal) * TOL / 5.0_defReal)

    data = newCoordData([0.35_defReal, -0.35_defReal, -1.35_defReal], [-ONE, ONE, ONE], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(2, data % elementIdx)
    @assertEqual(7.0_defReal * sqrt(3.0_defReal) / 20.0_defReal, data % d, 7.0_defReal * sqrt(3.0_defReal) * TOL / 20.0_defReal)

    data = newCoordData([-0.1_defReal, -0.25_defReal, 1.13_defReal], [10.0_defReal, 25.0_defReal, -13.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(3, data % elementIdx)
    @assertEqual(sqrt(894.0_defReal) / 100.0_defReal, data % d, sqrt(894.0_defReal) * TOL / 100.0_defReal)

    data = newCoordData([-0.33_defReal, 0.56_defReal, 1.27_defReal], [33.0_defReal, -56.0_defReal, -27.0_defReal], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(4, data % elementIdx)
    @assertEqual(sqrt(4954.0_defReal) / 100.0_defReal, data % d, sqrt(4954.0_defReal) * TOL / 100.0_defReal)
    
    ! Few points outside the mesh and not entering.
    data = newCoordData([-1.13_defReal, -0.8_defReal, 0.3_defReal], [-ONE, ZERO, ZERO], dMax = dMax)
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)
    
    data = newCoordData([0.65_defReal, 0.0_defReal, 1.25_defReal], [ONE, -ONE, ONE], dMax = dMax)
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


  @Test 
  subroutine test_distance_rescue() 
    type(coordData)          :: data
    real(defReal), parameter :: dMax = TWO, TOL = 1.0E-6
    type(ratint_t) :: val1, val2, val3

    !(-0.005 -0.005 -0.005)
    ! data = newCoordData([-0.015_defReal, -0.015_defReal, -0.015_defReal], & 
    !                     [-0.010_defReal, -0.010_defReal, -0.010_defReal], dMax = 1.0_defReal)
    ! call mesh % findHostElement(data)
    ! call mesh % distanceToNextFace(data)
    ! print *, data%elementIdx

    ! print *, 'thisTest'
    ! data = newCoordData([-1.005_defReal, -1.005_defReal, -1.005_defReal], [TWO, TWO, TWO], dMax = 2.0_defReal)
    ! call mesh % findHostElement(data)
    ! call mesh % distanceToNextFace(data)
    ! @assertEqual(1.7320508_defReal, data % d, 1.7320508_defReal * TOL)

    ! data = newCoordData([0.0_defReal, 0.0_defReal, 0.0_defReal],[ONE,ZERO,ZERO], dMax=1.0_defReal)
    ! call mesh1 % findHostElement(data)
    ! call mesh1 % distanceToNextFace(data)
    ! ! print *, data % elementIdx
    ! ! print *, data % d
    ! @assertEqual(0.5_defReal, data % d, 0.5_defReal * TOL) 



    ! data = newCoordData([0.4_defReal, 0.4_defReal, 0.4_defReal],[ONE,ONE,ONE], dMax=0.1_defReal*norm2([ONE,ONE,ONE]))
    ! call mesh1 % findHostElement(data)
    ! call mesh1 % distanceToNextFace(data)
    ! @assertEqual(0.173205_defReal, data % d, 0.173205_defReal * TOL) 

    ! Hasn't crossed
    data=newCoordData([0.4_defReal-1e-12,0.4_defReal-1e-12,0.4_defReal-1e-12],[ONE,ONE,ONE],dMax=0.1_defReal*norm2([ONE,ONE,ONE]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
    @assertEqual(INF, data % d, INF * TOL) 
    @assertEqual(1, data % elementIdx) 


    ! Technically hasn't crossed but here rescue is triggered
    data=newCoordData([0.4_defReal-1e-13,0.4_defReal-1e-13,0.4_defReal-1e-13],[ONE,ONE,ONE],dMax=0.1_defReal*norm2([ONE,ONE,ONE]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
    @assertEqual(INF, data % d, INF * TOL) 
    @assertEqual(1, data % elementIdx) 

    
    ! Hasn't crossed
    data=newCoordData([0.4_defReal,0.4_defReal,0.4_defReal],[ONE,ONE,ONE],dMax=0.099_defReal*norm2([ONE,ONE,ONE]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
    @assertEqual(INF, data % d, INF * TOL) 
    @assertEqual(1, data % elementIdx) 


    ! Crossed near the vertex
    data=newCoordData([0.45_defReal,0.45_defReal,0.45_defReal],[0.05_defReal,0.05_defReal,0.05_defReal], &
                      dMax=0.999999999999999998_defReal*norm2([0.05_defReal,0.05_defReal,0.05_defReal]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
    @assertEqual(0.08660254037_defReal, data % d, 0.08660254037_defReal * TOL) 
    @assertEqual(0, data % elementIdx) 


    ! Stopping on boundary face, pointing outwards on x+
    data=newCoordData([0.25_defReal,0.25_defReal,0.25_defReal],[ONE,ZERO,ZERO],dMax=0.25_defReal*norm2([ONE,ZERO,ZERO]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
  
    @assertEqual(0.25_defReal, data % d, 0.25_defReal * TOL)
    @assertEqual(0, data % elementIdx) 


    ! Stopping on boundary face, pointing outwards on z+
    data=newCoordData([0.25_defReal,0.25_defReal,0.25_defReal],[ZERO,ZERO,ONE],dMax=0.25_defReal*norm2([ZERO,ZERO,ONE]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
    @assertEqual(0.25_defReal, data % d, 0.25_defReal * TOL)
    @assertEqual(0, data % elementIdx) 



    ! Stopping on boundary face, pointing outwards on y+
    data=newCoordData([0.25_defReal,0.25_defReal,0.25_defReal],[ZERO,ONE,ZERO],dMax=0.25_defReal*norm2([ZERO,ONE,ZERO]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
    @assertEqual(0.25_defReal, data % d, 0.25_defReal * TOL)
    @assertEqual(0, data % elementIdx) 


    ! Stopping on boundary face, pointing outwards on x-
    data=newCoordData([0.25_defReal,0.25_defReal,0.25_defReal],[ONE,ZERO,ZERO],dMax=-0.75_defReal*norm2([ONE,ZERO,ZERO]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
    @assertEqual(0.75_defReal, data % d, 0.75_defReal * TOL)
    @assertEqual(0, data % elementIdx) 


    ! Stopping on boundary face, pointing outwards on z-
    data=newCoordData([0.25_defReal,0.25_defReal,0.25_defReal],[ZERO,ZERO,ONE],dMax=-0.75_defReal*norm2([ZERO,ZERO,ONE]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
    @assertEqual(0.75_defReal, data % d, 0.75_defReal * TOL)
    @assertEqual(0, data % elementIdx) 



    ! Stopping on boundary face, pointing outwards on y-
    data=newCoordData([0.25_defReal,0.25_defReal,0.25_defReal],[ZERO,ONE,ZERO],dMax=-0.75_defReal*norm2([ZERO,ONE,ZERO]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
    @assertEqual(0.75_defReal, data % d, 0.75_defReal * TOL)
    @assertEqual(0, data % elementIdx) 


    ! Centre to corner TR+
    data=newCoordData([0.0_defReal,0.0_defReal,0.0_defReal],[ONE,ONE,ONE],dMax=0.5_defReal*norm2([ONE,ONE,ONE]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
    @assertEqual(0.8660254_defReal, data % d, 0.8660254_defReal * TOL)
    @assertEqual(0, data % elementIdx) 


    ! Centre to corner TR-
    data=newCoordData([0.0_defReal,0.0_defReal,0.0_defReal],[ONE,ONE,ONE],dMax=0.5_defReal*norm2([ONE,ONE,-ONE]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
    @assertEqual(0.8660254_defReal, data % d, 0.8660254_defReal * TOL)
    @assertEqual(0, data % elementIdx) 


    ! Centre to corner TL+
    data=newCoordData([0.0_defReal,0.0_defReal,0.0_defReal],[ONE,ONE,ONE],dMax=0.5_defReal*norm2([-ONE,ONE,ONE]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
    @assertEqual(0.8660254_defReal, data % d, 0.8660254_defReal * TOL)
    @assertEqual(0, data % elementIdx) 


    ! Centre to corner TL-
    data=newCoordData([0.0_defReal,0.0_defReal,0.0_defReal],[ONE,ONE,ONE],dMax=0.5_defReal*norm2([-ONE,ONE,-ONE]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
    @assertEqual(0.8660254_defReal, data % d, 0.8660254_defReal * TOL)
    @assertEqual(0, data % elementIdx) 


    ! Centre to corner BL+
    data=newCoordData([0.0_defReal,0.0_defReal,0.0_defReal],[ONE,ONE,ONE],dMax=0.5_defReal*norm2([-ONE,-ONE,ONE]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
    @assertEqual(0.8660254_defReal, data % d, 0.8660254_defReal * TOL)
    @assertEqual(0, data % elementIdx) 


    ! Centre to corner BL-
    data=newCoordData([0.0_defReal,0.0_defReal,0.0_defReal],[ONE,ONE,ONE],dMax=0.5_defReal*norm2([-ONE,-ONE,-ONE]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
    @assertEqual(0.8660254_defReal, data % d, 0.8660254_defReal * TOL)
    @assertEqual(0, data % elementIdx) 



    ! Centre to corner BR+
    data=newCoordData([0.0_defReal,0.0_defReal,0.0_defReal],[ONE,ONE,ONE],dMax=0.5_defReal*norm2([ONE,-ONE,ONE]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
    @assertEqual(0.8660254_defReal, data % d, 0.8660254_defReal * TOL)
    @assertEqual(0, data % elementIdx) 

    ! Centre to corner BR-
    data=newCoordData([0.0_defReal,0.0_defReal,0.0_defReal],[ONE,ONE,ONE],dMax=0.5_defReal*norm2([ONE,-ONE,-ONE]))
    call mesh1 % findHostElement(data)
    call mesh1 % distanceToNextFace(data)
    @assertEqual(0.8660254_defReal, data % d, 0.8660254_defReal * TOL)
    @assertEqual(0, data % elementIdx) 




    ! Tests in other mesh (multiple cubes)

    ! exits through top left edge
    data=newCoordData([0.5_defReal,0.5_defReal,0.5_defReal],[ONE,ONE,ZERO], &
                      dMax=0.5_defReal*norm2([ONE,ONE,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.7071067811_defReal, data % d, 0.7071067811_defReal * TOL) 
    @assertEqual(0, data % elementIdx) 


    !Staying in the same cube, test for element id
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[ZERO,ZERO,ZERO], &
                      dMax=0.5_defReal*norm2([ZERO,ZERO,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(1, data % elementIdx) 

    ! print *, data % elementIdx
    ! print *, data % d


    ! Crossing test for near edge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    !NOTE! THIS DOESNT TRIGGER RESCUE, AND ISN'T HANDLED, REPORTS 3 BUT ACTUALLY EXITS
    data=newCoordData([-0.5_defReal,0.5_defReal,-0.5_defReal],[0.50000000000001_defReal,0.49999999999999_defReal,ZERO], &
                      dMax=2.0_defReal*norm2([0.50000000000001_defReal,0.49999999999999_defReal,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.7071067811_defReal, data % d, 0.7071067811_defReal * TOL) 
    @assertEqual(3, data % elementIdx) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! Moving but staying in same element
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,-0.5_defReal], &
                      dMax=0.01_defReal*norm2([0.5_defReal,0.5_defReal,-0.5_defReal]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(INF, data % d, INF * TOL) 
    @assertEqual(1, data % elementIdx) 

    
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,-0.5_defReal], &
                      dMax=1.5_defReal*norm2([0.5_defReal,0.5_defReal,-0.5_defReal]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)

    ! print *, data % elementIdx
    ! print *, data % d

    
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,-0.5_defReal], &
                      dMax=1.0_defReal*norm2([0.5_defReal,0.5_defReal,-0.5_defReal]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.8660254037_defReal, data % d, 0.8660254037_defReal * TOL) 
    @assertEqual(3, data % elementIdx) 


!!through edge and past it
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,ZERO], &
                      dMax=1.5_defReal*norm2([0.5_defReal,0.5_defReal,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(3, data % elementIdx) 

    ! print *, data % elementIdx
    ! print *, data % d

    !!stopping directly on an edge bordering another element
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[ONE,ONE,ZERO], &
                      dMax=0.5_defReal*norm2([ONE,ONE,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(3, data % elementIdx) 

    ! print *, data % elementIdx
    ! print *, data % d

    @assertEqual(0.7071067811_defReal, data % d, 0.7071067811_defReal * TOL) 

!===================================
    ! Stopping just past the edge
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,ZERO], &
                      dMax=1.1_defReal*norm2([0.5_defReal,0.5_defReal,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.7071067811_defReal, data % d, 0.7071067811_defReal * TOL) 
    @assertEqual(3, data % elementIdx) 


    data=newCoordData([0.5_defReal,-0.5_defReal,0.5_defReal],[ZERO,ZERO,ZERO], &
                      dMax=0.0_defReal*norm2([ZERO,ZERO,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)


    !! stopping directly on a face
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,ZERO,ZERO], &
                      dMax=1.0_defReal*norm2([0.5_defReal,ZERO,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.5_defReal, data % d, 0.5_defReal * TOL) 
    @assertEqual(4, data % elementIdx) 
    !@assertEqual(0.8660254037_defReal, data % elementIdx) 

    

    !Stopping just after a face 
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,ZERO,ZERO], &
                      dMax=1.0_defReal+1e-12*norm2([0.5_defReal,ZERO,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.5_defReal, data % d, 0.5_defReal * TOL) 
    @assertEqual(4, data % elementIdx) 


    ! Through the corner and into the corner cube (diagnonally), stopping at vertex
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,-0.5_defReal], &
                      dMax=1.0_defReal*norm2([0.5_defReal,0.5_defReal,-0.5_defReal]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.8660254037_defReal, data % d, 0.8660254037_defReal * TOL) 
    @assertEqual(3, data % elementIdx) 


    ! Through the corner and into the corner cube (diagnonally), stopping beyond vertex
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,-0.5_defReal], &
                      dMax=1.2_defReal*norm2([0.5_defReal,0.5_defReal,-0.5_defReal]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.8660254037_defReal, data % d, 0.8660254037_defReal * TOL) 
    @assertEqual(3, data % elementIdx) 

  
  
    !Sending particle epsilon close to an edge (before)
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,ZERO], &
                      dMax=(1.0_defReal-1e-12)*norm2([0.5_defReal,0.5_defReal,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(INF, data % d, INF * TOL) 
    @assertEqual(1, data % elementIdx) 


    !Sending particle  just <epsilon close to an edge (before)
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,ZERO], &
                      dMax=(1.0_defReal-1e-13)*norm2([0.5_defReal,0.5_defReal,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(INF, data % d, INF * TOL) 
    @assertEqual(1, data % elementIdx) 



    !Sending particle <epsilon close to an edge (before)
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,ZERO], &
                      dMax=(1.0_defReal-1e-14)*norm2([0.5_defReal,0.5_defReal,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(INF, data % d, INF * TOL) 
    @assertEqual(1, data % elementIdx) 


    !Sending particle <epsilon close to an edge (before)
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,ZERO], &
                      dMax=(1.0_defReal-1e-15)*norm2([0.5_defReal,0.5_defReal,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(INF, data % d, INF * TOL) 
    @assertEqual(1, data % elementIdx) 



    !Sending particle <epsilon close to an edge (before)
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,ZERO], &
                      dMax=(1.0_defReal-1e-16)*norm2([0.5_defReal,0.5_defReal,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(INF, data % d, INF * TOL) 
    @assertEqual(1, data % elementIdx) 



    ! !Sending particle <epsilon close to an edge (before):: NOTE FAILS, assumes crossed at e-17
    ! data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,ZERO], &
    !                   dMax=(1.0_defReal-1e-17)*norm2([0.5_defReal,0.5_defReal,ZERO]))
    ! call mesh % findHostElement(data)
    ! call mesh % distanceToNextFace(data)
    ! @assertEqual(INF, data % d, INF * TOL) 
    ! @assertEqual(1, data % elementIdx) 



    ! !Sending particle <epsilon close to an edge (before)
    ! data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,ZERO], &
    !                   dMax=(1.0_defReal-1e-18)*norm2([0.5_defReal,0.5_defReal,ZERO]))
    ! call mesh % findHostElement(data)
    ! call mesh % distanceToNextFace(data)
    ! @assertEqual(INF, data % d, INF * TOL) 
    ! @assertEqual(1, data % elementIdx) 



    !Sending particle epsilon close to a vertex (after)
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,-0.5_defReal], &
                      dMax=(1.0_defReal)*norm2([0.5_defReal,0.5_defReal,-0.5_defReal]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.8660254037_defReal, data % d, 0.8660254037_defReal * TOL) 
    @assertEqual(3, data % elementIdx)

!=================================
    !Sending particle just >epsilon close to an edge (after)
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[ONE,ONE,ZERO], &
                      dMax=(0.5_defReal+1e-11)*norm2([ONE,ONE,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.7071067811_defReal, data % d, 0.7071067811_defReal * TOL) 
    @assertEqual(3, data % elementIdx)


    !Sending particle epsilon close to an edge (after)
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[ONE,ONE,ZERO], &
                      dMax=(0.5_defReal+1e-12)*norm2([ONE,ONE,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.7071067811_defReal, data % d, 0.7071067811_defReal * TOL) 
    @assertEqual(3, data % elementIdx)


    !Sending particle just <epsilon close to an edge (after)
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[ONE,ONE,ZERO], &
                      dMax=(0.5_defReal+1e-13)*norm2([ONE,ONE,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.7071067811_defReal, data % d, 0.7071067811_defReal * TOL) 
    @assertEqual(3, data % elementIdx)


    !Sending particle <epsilon close to an edge (after)
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[ONE,ONE,ZERO], &
                      dMax=(0.5_defReal+1e-14)*norm2([ONE,ONE,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.7071067811_defReal, data % d, 0.7071067811_defReal * TOL) 
    @assertEqual(3, data % elementIdx)



    !Sending particle <epsilon close to an edge (after)
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[ONE,ONE,ZERO], &
                      dMax=(0.5_defReal+1e-15)*norm2([ONE,ONE,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.7071067811_defReal, data % d, 0.7071067811_defReal * TOL) 
    @assertEqual(3, data % elementIdx)



    !Sending particle <epsilon close to an edge (after)
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[ONE,ONE,ZERO], &
                      dMax=(0.5_defReal+1e-16)*norm2([ONE,ONE,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.7071067811_defReal, data % d, 0.7071067811_defReal * TOL) 
    @assertEqual(3, data % elementIdx)



    !Sending particle <epsilon close to an edge (after)
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[ONE,ONE,ZERO], &
                      dMax=(0.5_defReal+1e-17)*norm2([ONE,ONE,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.7071067811_defReal, data % d, 0.7071067811_defReal * TOL) 
    @assertEqual(3, data % elementIdx)


     !Sending particle <epsilon close to an edge (after)
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[ONE,ONE,ZERO], &
                      dMax=(0.5_defReal+1e-18)*norm2([ONE,ONE,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.7071067811_defReal, data % d, 0.7071067811_defReal * TOL) 
    @assertEqual(3, data % elementIdx)


     !Sending particle <epsilon close to an edge (after)
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[ONE,ONE,ZERO], &
                      dMax=(0.5_defReal+1e-19)*norm2([ONE,ONE,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.7071067811_defReal, data % d, 0.7071067811_defReal * TOL) 
    @assertEqual(3, data % elementIdx)



     !Sending particle <epsilon close to an edge (after)
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[ONE,ONE,ZERO], &
                      dMax=(0.5_defReal+1e-20)*norm2([ONE,ONE,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.7071067811_defReal, data % d, 0.7071067811_defReal * TOL) 
    @assertEqual(3, data % elementIdx)


     !Sending particle <epsilon close to an edge (after)
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[ONE,ONE,ZERO], &
                      dMax=(0.5_defReal+1e-21)*norm2([ONE,ONE,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.7071067811_defReal, data % d, 0.7071067811_defReal * TOL) 
    @assertEqual(3, data % elementIdx)


    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[ONE,ONE,ZERO], &
                      dMax=(0.5_defReal)*norm2([ONE,ONE,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    @assertEqual(0.7071067811_defReal, data % d, 0.7071067811_defReal * TOL) 
    @assertEqual(3, data % elementIdx)


    !print *, '=================================================================='



    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,ZERO], &
                      dMax=(1.0_defReal)*norm2([0.5_defReal,0.5_defReal,ZERO]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    ! print *, data % elementIdx
    ! print *, data % d

    @assertEqual(0.7071067811_defReal, data % d, 0.7071067811_defReal * TOL) 
    @assertEqual(3, data % elementIdx)

    ! epsilon close through the edge and into the corner cube (diagnonally), stopping well inside other cube
    data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,-0.5_defReal], &
                      dMax=1.1_defReal*norm2([0.5_defReal,0.5_defReal,-0.5_defReal]))
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)
    ! print *, data % elementIdx
    ! print *, data % d

    @assertEqual(0.8660254037_defReal, data % d, 0.8660254037_defReal * TOL) 
    @assertEqual(3, data % elementIdx)


    ! data=newCoordData([0.5_defReal,0.5_defReal,-0.5_defReal],[ZERO,ZERO,ZERO], &
    !                   dMax=0.0_defReal*norm2([ZERO,ZERO,ZERO]))
    ! call mesh % findHostElement(data)
    ! call mesh % distanceToNextFace(data)
    ! print *, data % elementIdx
    ! print *, data % d

    ! @assertEqual(0.8660254037_defReal, data % d, 0.8660254037_defReal * TOL) 



    ! data=newCoordData([-0.5_defReal+1e-12,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,-0.5_defReal], &
    !                   dMax=1.2_defReal*norm2([0.5_defReal,0.5_defReal,-0.5_defReal]))
    ! call mesh % findHostElement(data)
    ! call mesh % distanceToNextFace(data)
    ! print *, data % elementIdx
    ! print *, data % d

    ! @assertEqual(0.8660254037_defReal, data % d, 0.8660254037_defReal * TOL) 





    ! only true under wrong function
    ! data=newCoordData([-0.5_defReal,-0.5_defReal,0.5_defReal],[0.5_defReal,0.5_defReal,-0.5_defReal], &
    !                   dMax=4.0_defReal*norm2([0.5_defReal,0.5_defReal,-0.5_defReal]))
    ! call mesh % findHostElement(data)
    ! call mesh % distanceToNextFace(data)
    ! @assertEqual(2.598076_defReal, data % d, 2.598076_defReal * TOL) 


    ! val1 = convert_ieee(0.4_defReal-1e-12+(1/norm2([ONE,ONE,ONE]))*0.1_defReal*norm2([ONE,ONE,ONE]))
    ! ! val2 = convert_ieee(0.4_defReal-1e-12+0.1_defReal*norm2([ONE,ONE,ONE]))

    ! ! call printRatInt(val1)
    ! ! print *, '===='
    ! ! call printRatInt(val2)
    ! val3 = convert_ieee(0.5_defReal)

    ! ! print *, val1 > val2
    ! print *, val3 >= val1 
    ! ! print *, val3 >= val2
    
    ! print *, evaluate(val1)
    ! print *, evaluate(val2)










    ! print *, 'vals!'
    ! print *, data%r 
    ! print *, data%u 
    ! print *, data%dMax



  end subroutine test_distance_rescue

end module OpenFOAMMesh_iTest