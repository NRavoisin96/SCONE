module meshShelf_iTest
  
  use charMap_class,    only : charMap
  use dictionary_class, only : dictionary
  use dictParser_func,  only : charToDict
  use funit
  use mesh_inter,       only : mesh
  use meshShelf_class,  only : meshShelf
  use numPrecision
  
  implicit none
  
  ! Parameters.
  character(*), parameter :: MESHES_DEF = &
  " testMesh1 { id 11; type OpenFOAMMesh; path ./IntegrationTestFiles/Geometry/Meshes/OpenFOAM/testMesh1/; fills (fuel);} &
   &testMesh2 { id 21; type OpenFOAMMesh; path ./IntegrationTestFiles/Geometry/Meshes/OpenFOAM/testMesh2/; fills (water);}"
  
  ! Variables.
  type(charMap)   :: mats
  type(meshShelf) :: meshes
  
contains
  
  !!
  !! Build shelf.
  !!
@Before
  subroutine setUp()
    character(nameLen) :: materialName
    type(dictionary)   :: dict
    
    materialName = 'fuel'
    call mats % add(materialName, 1)

    materialName = 'water'
    call mats % add(materialName, 2)

    call charToDict(dict, MESHES_DEF)
    call meshes % init(dict, mats)

  end subroutine setUp
  
  !!
  !! Clean after tests.
  !!
@After
  subroutine cleanUp()

    call mats % kill()
    call meshes % kill()

  end subroutine cleanUp
  
  !!
  !! Test shelf.
  !!
@Test
  subroutine test_get()
    class(mesh), pointer :: ptr
    integer(shortInt)    :: idx

    ! Mesh ID 11.
    idx = meshes % getMeshIdx(11)
    ptr => meshes % getMeshPtr(idx)
    @assertEqual(11, ptr % getId())
    @assertEqual(11, meshes % getMeshId(idx))

    ! Mesh ID 21.
    idx = meshes % getMeshIdx(21)
    ptr => meshes % getMeshPtr(idx)
    @assertEqual(21, ptr % getId())
    @assertEqual(21, meshes % getMeshId(idx))   

  end subroutine test_get

end module meshShelf_iTest