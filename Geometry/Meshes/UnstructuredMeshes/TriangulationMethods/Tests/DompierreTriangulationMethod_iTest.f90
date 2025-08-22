module DompierreTriangulationMethod_iTest

  use charMap_class,      only : charMap
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use element_class,      only : elementBox
  use funit
  use numPrecision
  use OpenFOAMMesh_class, only : OpenFOAMMesh
  use publicObjects,      only : coordData, newCoordData
  use universalVariables
  
  implicit none
  
  ! Parameters.
  character(*), parameter :: MESH_DEF = &
  " id 2; type OpenFOAMMesh; path ./IntegrationTestFiles/Geometry/Meshes/OpenFOAM/testMesh/;&
  & accelerationMethod {type none;} triangulationMethod Dompierre; fills (fuel);"
  
  ! Variables.
  type(charMap)      :: mats
  type(OpenFOAMMesh) :: mesh

contains
  
  !!
  !! Import the mesh.
  !!
@Before
  subroutine setUp()
    type(dictionary)   :: dict
    character(nameLen) :: name
    character(pathLen) :: path

    name = 'fuel'
    call mats % add(name, 1)
    
    call charToDict(dict, MESH_DEF)
    call dict % get(path, 'path')
    call mesh % init(trim(path), dict, mats)
  
  end subroutine setUp
  
  !!
  !! Clean after tests.
  !!
@After
  subroutine cleanUp()

    call mats % kill()
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
    
    ! Test number of vertices.
    @assertEqual(18, mesh % getVerticesNumber())
    ! Test number of faces.
    @assertEqual(94, mesh % getFacesNumber(.true.))
    ! Test number of internal faces.
    @assertEqual(62, mesh % getInternalFacesNumber(.true.))
    ! Test number of edges.
    @assertEqual(72, mesh % getEdgesNumber())
    ! Test number of tetrahedra.
    @assertEqual(24, mesh % getElementsNumber(.true.))


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
    @assertEqual(18, data % elementIdx)

    r = [0.02_defReal, 0.97_defReal, -0.5_defReal]
    u = [ZERO, ONE, ZERO]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(19, data % elementIdx)

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
    @assertEqual(12, data % elementIdx)

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
    @assertEqual(17, data % elementIdx)

    u = [ONE, 2.0_defReal, ZERO]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(18, data % elementIdx)

    u = [-2.0_defReal, ONE, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(13, data % elementIdx)

    u = [-ONE, 2.0_defReal, ZERO]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(15, data % elementIdx)
    
    u = [2.0_defReal, -ONE, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(27, data % elementIdx)

    u = [ONE, -2.0_defReal, ZERO]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(25, data % elementIdx)
    
    u = [-2.0_defReal, -ONE, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(10, data % elementIdx)

    u = [-ONE, -2.0_defReal, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(10, data % elementIdx)

    ! A point on a boundary edge. Different directions.
    r = [ONE, ZERO, 0.5_defReal]

    ! Points inside the mesh.
    u = [-2.0_defReal, ONE, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(20, data % elementIdx)

    u = [-ONE, 2.0_defReal, ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(20, data % elementIdx)

    u = [-2.0_defReal, -ONE, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(28, data % elementIdx)

    u = [-ONE, -2.0_defReal, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(28, data % elementIdx)

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
    @assertEqual(10, data % elementIdx)

    u = [-2.0_defReal, -ONE, -5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(10, data % elementIdx)

    u = [-ONE, -2.0_defReal, -5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(10, data % elementIdx)

    u = [-ONE, -3.0_defReal, -3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(10, data % elementIdx)

    u = [-ONE, -2.0_defReal, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(10, data % elementIdx)

    u = [-2.0_defReal, -ONE, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(10, data % elementIdx)

    u = [-ONE, 3.0_defReal, -3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(16, data % elementIdx)

    u = [-ONE, 2.0_defReal, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(16, data % elementIdx)

    u = [-3.0_defReal, ONE, -3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(13, data % elementIdx)

    u = [-2.0_defReal, ONE, -5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(13, data % elementIdx)

    u = [2.0_defReal, ONE, -5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(17, data % elementIdx)

    u = [2.0_defReal, ONE, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(21, data % elementIdx)

    u = [ONE, 2.0_defReal, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(21, data % elementIdx)

    u = [ONE, 2.0_defReal, -5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(18, data % elementIdx)

    u = [ONE, -2.0_defReal, -ONE]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(25, data % elementIdx)

    u = [ONE, -2.0_defReal, -5.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(25, data % elementIdx)

    u = [ONE, -3.0_defReal, -3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(25, data % elementIdx)

    u = [3.0_defReal, -ONE, -3.0_defReal]
    data = newCoordData(r, u)
    call mesh % findHostElement(data)
    @assertEqual(28, data % elementIdx)

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
    real(defReal), parameter    :: TOL = 1.0E-6

  
  end subroutine test_distance

end module DompierreTriangulationMethod_iTest