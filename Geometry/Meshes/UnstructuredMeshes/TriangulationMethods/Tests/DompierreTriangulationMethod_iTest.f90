module DompierreTriangulationMethod_iTest

  use coord_class,        only : coord
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use element_class,      only : elementBox
  use funit
  use numPrecision
  use OpenFOAMMesh_class, only : OpenFOAMMesh
  use universalVariables
  use vertex_class,       only : vertexBox
  
  implicit none
  
  ! Parameters.
  character(*), parameter :: MESH_DEF = &
  " id 2; type OpenFOAMMesh; path ./IntegrationTestFiles/Geometry/Meshes/OpenFOAM/testMesh/; triangulationMethod Dompierre;"
  ! Variables.
  type(OpenFOAMMesh) :: mesh
  type(coord)        :: coords

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
    
    ! Test number of vertices.
    @assertEqual(18, mesh % getVerticesNumber())
    ! Test number of faces.
    @assertEqual(94, mesh % getFacesNumber())
    ! Test number of internal faces.
    @assertEqual(62, mesh % getInternalFacesNumber())
    ! Test number of edges.
    @assertEqual(72, mesh % getEdgesNumber())
    ! Test number of tetrahedra.
    @assertEqual(24, mesh % getElementsNumber())


  end subroutine test_info
  
  !!
  !! Test inside / outside determination.
  !!
@Test
  subroutine test_inside()
    real(defReal), dimension(3) :: r, u
    type(coord)                 :: coords
    type(vertexBox), dimension(4) :: vertices
    integer(shortInt)             :: i
    
    ! Few points inside.
    r = [0.31_defReal, 0.42_defReal, 0.13_defReal]
    u = [ONE, ZERO, ZERO]
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(18, coords % getElementIdx())

    r = [0.02_defReal, 0.97_defReal, -0.5_defReal]
    u = [ZERO, ONE, ZERO]
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(19, coords % getElementIdx())

    ! Few points outside.
    r = [1.2_defReal, 0.8_defReal, 0.0_defReal]
    u = [ZERO, ONE, ZERO]
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())
    
    r = [0.1_defReal, 0.1_defReal, 1.13_defReal]
    u = [ZERO, ZERO, -ONE]
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    ! Few more difficult points.
    
    ! A point on a face.
    r = [-1.0_defReal, 0.1_defReal, 0.1_defReal]

    ! Points into the mesh.
    call coords % setPosition(r)
    u = [ONE, ZERO, ZERO]
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(12, coords % getElementIdx())

    ! Points away from the mesh.
    call coords % setPosition(r)
    u = [-ONE, ZERO, ZERO]
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    ! A point on an internal edge. Different directions.
    r = ZERO
    u = [2.0_defReal, ONE, ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(17, coords % getElementIdx())

    u = [ONE, 2.0_defReal, ZERO]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(18, coords % getElementIdx())

    u = [-2.0_defReal, ONE, -ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(13, coords % getElementIdx())

    u = [-ONE, 2.0_defReal, ZERO]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(15, coords % getElementIdx())
    
    u = [2.0_defReal, -ONE, -ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(27, coords % getElementIdx())

    u = [ONE, -2.0_defReal, ZERO]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(25, coords % getElementIdx())
    
    u = [-2.0_defReal, -ONE, ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(10, coords % getElementIdx())

    u = [-ONE, -2.0_defReal, -ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(10, coords % getElementIdx())

    ! A point on a boundary edge. Different directions.
    r = [ONE, ZERO, 0.5_defReal]

    ! Points inside the mesh.
    u = [-2.0_defReal, ONE, ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(20, coords % getElementIdx())

    u = [-ONE, 2.0_defReal, ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(20, coords % getElementIdx())

    u = [-2.0_defReal, -ONE, -ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(28, coords % getElementIdx())

    u = [-ONE, -2.0_defReal, -ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(28, coords % getElementIdx())

    ! Points away from the mesh.
    u = [2.0_defReal, ONE, ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [ONE, 2.0_defReal, ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [2.0_defReal, -ONE, -ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [ONE, -2.0_defReal, -ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    ! A point on a boundary vertex. Different directions.
    r = [ZERO, ZERO, ONE]

    ! Points into the mesh.
    u = [-3.0_defReal, -ONE, -3.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(10, coords % getElementIdx())

    u = [-2.0_defReal, -ONE, -5.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(10, coords % getElementIdx())

    u = [-ONE, -2.0_defReal, -5.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(10, coords % getElementIdx())

    u = [-ONE, -3.0_defReal, -3.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(10, coords % getElementIdx())

    u = [-ONE, -2.0_defReal, -ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(10, coords % getElementIdx())

    u = [-2.0_defReal, -ONE, -ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(10, coords % getElementIdx())

    u = [-ONE, 3.0_defReal, -3.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(16, coords % getElementIdx())

    u = [-ONE, 2.0_defReal, -ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(16, coords % getElementIdx())

    u = [-3.0_defReal, ONE, -3.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(13, coords % getElementIdx())

    u = [-2.0_defReal, ONE, -5.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(13, coords % getElementIdx())

    u = [2.0_defReal, ONE, -5.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(17, coords % getElementIdx())

    u = [2.0_defReal, ONE, -ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(21, coords % getElementIdx())

    u = [ONE, 2.0_defReal, -ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(21, coords % getElementIdx())

    u = [ONE, 2.0_defReal, -5.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(18, coords % getElementIdx())

    u = [ONE, -2.0_defReal, -ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(25, coords % getElementIdx())

    u = [ONE, -2.0_defReal, -5.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(25, coords % getElementIdx())

    u = [ONE, -3.0_defReal, -3.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(25, coords % getElementIdx())

    u = [3.0_defReal, -ONE, -3.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(28, coords % getElementIdx())

    ! Points away from the mesh.
    u = [-3.0_defReal, -ONE, 3.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [-2.0_defReal, -ONE, 5.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [-ONE, -2.0_defReal, 5.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [-ONE, -3.0_defReal, 3.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [-ONE, -2.0_defReal, ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [-2.0_defReal, -ONE, ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [-ONE, 3.0_defReal, 3.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [-ONE, 2.0_defReal, ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [-3.0_defReal, ONE, 3.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [-2.0_defReal, ONE, 5.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [2.0_defReal, ONE, 5.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [2.0_defReal, ONE, ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [ONE, 2.0_defReal, ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [ONE, 2.0_defReal, 5.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [ONE, -2.0_defReal, ONE]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [ONE, -2.0_defReal, 5.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [ONE, -3.0_defReal, 3.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
    call coords % setDirection(u)
    call mesh % findHostElement(coords)
    @assertEqual(0, coords % getElementIdx())

    u = [3.0_defReal, -ONE, 3.0_defReal]
    u = u / norm2(u)
    call coords % setPosition(r)
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

    
  
  end subroutine test_distance

end module DompierreTriangulationMethod_iTest