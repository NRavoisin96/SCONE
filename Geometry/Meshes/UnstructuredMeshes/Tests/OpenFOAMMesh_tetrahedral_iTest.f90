module OpenFOAMMesh_tetrahedral_iTest

  use charMap_class,      only : charMap
  use coord_class,        only : coord
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use funit
  use numPrecision
  use OpenFOAMMesh_class, only : OpenFOAMMesh
  use universalVariables
  
  implicit none
  
  ! Parameters.
  character(*), parameter :: MESH_DEF = &
  "id 15; type OpenFOAMMesh; path ./IntegrationTestFiles/Geometry/Meshes/OpenFOAM/StanfordBunny_LowPoly/; fills (water);"
  
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

    name = 'water'
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
    @assertEqual(15, mesh % getId())
    call mesh % setId(3)
    @assertEqual(3, mesh % getId())

  end subroutine test_misc
  
  !!
  !! Test mesh information.
  !!
@Test
  subroutine test_info()
    real(defReal)     :: TOL = 1.0E-6
    integer(shortInt) :: i
    
    ! Test number of vertices.
    @assertEqual(154, mesh % getVerticesNumber())
    ! Test number of faces.
    @assertEqual(904, mesh % getFacesNumber())
    ! Test number of internal faces.
    @assertEqual(612, mesh % getInternalFacesNumber())
    ! Test number of edges.
    @assertEqual(678, mesh % getEdgesNumber())
    ! Test number of tetrahedra.
    @assertEqual(379, mesh % getElementsNumber())

  end subroutine test_info

end module OpenFOAMMesh_tetrahedral_iTest