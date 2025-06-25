module meshUniverse_iTest
  
  use numPrecision
  use universalVariables
  use genericProcedures
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use coord_class,        only : coord
  use surfaceShelf_class, only : surfaceShelf
  use meshShelf_class,    only : meshShelf
  use cellShelf_class,    only : cellShelf
  use meshUniverse_class, only : meshUniverse
  use funit
  
  implicit none
  
  ! Parameters.
  character(*), parameter :: SURFS_DEF = &
  " surf1 { id 1; type sphere; origin (0.0 0.0 0.0); radius 2;}&
  & surf2 { id 2; type box; origin (0.0 0.0 0.0); halfwidth (1.5 1.5 1.5);}"
  character(*), parameter :: CELL_DEF = &
  " cell1 {id 1; type simpleCell; surfaces (-1); filltype mat; material water;}"
  
  character(*), parameter :: MESH_DEF = &
  " testMesh {id 8; type OpenFOAMMesh; path ./IntegrationTestFiles/Geometry/Meshes/OpenFOAM/testMesh/;}"
  ! 
  ! Note that rotation is such that following axis transformation applies:
  !   x -> z
  !   y -> -y
  !   z -> x
  ! 
  character(*), parameter :: UNI_DEF = &
  "id 1; type meshUniverse; origin (2.0 0.0 0.0); rotation (90.0 90.0 90.0); cell 1; mesh 8; &
   fills (fuel water fuel water);"
  ! Variables.
  type(surfaceShelf) :: surfs
  type(meshShelf)    :: meshes
  type(cellShelf)    :: cells
  type(charMap)      :: mats
  type(meshUniverse) :: uni

contains
  !!
  !! Setup environment.
  !!
@Before
  subroutine setUp()
    character(nameLen)                           :: name
    integer(shortInt), dimension(:), allocatable :: fill
    type(dictionary)                             :: dict

    ! Add materials.
    name = 'water'
    call mats % add(name, 1)
    name = 'fuel'
    call mats % add(name, 2)
    ! Build surfaces.
    call charToDict(dict, SURFS_DEF)
    call surfs % init(dict)
    call dict % kill()
    ! Build cells.
    call charToDict(dict, CELL_DEF)
    call cells % init(dict, surfs, mats)
    call dict % kill()
    
    ! Build meshes.
    call charToDict(dict, MESH_DEF)
    call meshes % init(dict)
    call dict % kill()
    ! Build universe.
    call charToDict(dict, UNI_DEF)
    call uni % init(dict, mats, fill, cells, surfs, meshes)
    call dict % kill()
    ! Set index.
    call uni % setIdx(18)
    ! Verify fills.
    @assertEqual([1, 2, 1, 2, 1], fill)

  end subroutine setUp
  
  !!
  !! Clean environment.
  !!
@After
  subroutine cleanUp()

    call surfs % kill()
    call meshes % kill()
    call cells % kill()
    call mats % kill()
    call uni % kill()

  end subroutine cleanUp
  !!
  !! Test miscellaneous functionality (of generic universe).
  !!
@Test
  subroutine test_misc()

    ! Get ID.
    @assertEqual(1, uni % id())
    ! Set ID.
    call uni % setId(71)
    @assertEqual(71, uni % id())

  end subroutine test_misc
  !!
  !! Test entering the universe.
  !!
@Test
  subroutine test_enter()
    type(coord)                    :: new
    real(defReal), dimension(3)    :: r_ref, u_ref, r, u
    real(defReal), dimension(3, 3) :: rotationMatrix
    real(defReal), parameter       :: TOL = 1.0E-7_defReal

    ! ** Enter into local cell 1.
    r = [ZERO, ZERO, 3.0_defReal]
    u = [ZERO, ZERO, ONE]
    call uni % enter(r, u, new)
    ! Verify location. Note that particle is at the mesh boundary and is going out due to its 
    ! direction.
    r_ref = [ONE, ZERO, ZERO]
    u_ref = [ONE, ZERO, ZERO]
    @assertEqual(r_ref, new % getPosition(), TOL)
    @assertEqual(u_ref, new % getDirection(), TOL)
    @assertEqual(18, new % getUniIdx())
    @assertEqual(1, new % getLocalId())
    @assertEqual(cells % getIdx(1), new % getCellIdx())
    
    ! Change direction and check that the particle is in the mesh.
    r = [ZERO, ZERO, 3.0_defReal]
    u = [ZERO, ZERO, -ONE]
    call uni % enter(r, u, new)
    u_ref = [-ONE, ZERO, ZERO]
    @assertEqual(r_ref, new % getPosition(), TOL)
    @assertEqual(u_ref, new % getDirection(), TOL)
    @assertEqual(18, new % getUniIdx())
    @assertEqual(4, new % getLocalId())
    @assertEqual(cells % getIdx(1), new % getCellIdx())
    ! Verify rotation settings in coord.
    ! * Do it only once.
    @assertTrue(new % getIsRotated())
    rotationMatrix = new % getRotationMatrix()
    @assertEqual([ZERO, ZERO,  ONE], rotationMatrix(1, :), TOL)
    @assertEqual([ZERO, -ONE, ZERO], rotationMatrix(2, :), TOL)
    @assertEqual([ONE , ZERO, ZERO], rotationMatrix(3, :), TOL)

  end subroutine test_enter
  !!
  !! Test distance calculations.
  !!
@Test
  subroutine test_distance()
    real(defReal)            :: d, ref, maxDist
    integer(shortInt)        :: surfIdx
    type(coord)              :: pos
    real(defReal), parameter :: TOL = 1.0E-7_defReal
    
    ! Initialise maxDist.
    maxDist = ONE
    
    ! ** In local cell 1 distance to cell boundary.
    call pos % setPosition([-1.2_defReal, ZERO, ZERO])
    call pos % setDirection([-ONE, ZERO, ZERO])
    call pos % setEndPosition(pos % getPosition() + pos % getDirection() * maxDist)
    call pos % setUniIdx(18)
    call pos % setCellIdx(cells % getIdx(1))
    call pos % setLocalId(1)
    call uni % distance(pos, d, surfIdx)
    ref = 0.8_defReal
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(surfs % getIdx(1), surfIdx)
    
    ! ** In local cell 1 distance to mesh boundary.
    call pos % setPosition([-1.2_defReal, ZERO, ZERO])
    call pos % setDirection([ONE, ZERO, ZERO])
    call pos % setEndPosition(pos % getPosition() + pos % getDirection() * maxDist)
    call pos % setUniIdx(18)
    call pos % setCellIdx(cells % getIdx(1))
    call pos % setLocalId(1)
    call pos % setElementIdx(0)
    call uni % distance(pos, d, surfIdx)
    ref = 0.2_defReal
    @assertEqual(ref, d, TOL * ref)
    @assertEqual(0, surfIdx)
    
    ! ** In local cell 1. Leaving mesh and on mesh boundary.
    call pos % setPosition([ZERO, ONE, ZERO])
    call pos % setDirection([ZERO, ONE, ZERO])
    call pos % setEndPosition(pos % getPosition() + pos % getDirection() * maxDist)
    call pos % setUniIdx(18)
    call pos % setCellIdx(cells % getIdx(1))
    call pos % setLocalId(1)
    call pos % setElementIdx(0)
    call uni % distance(pos, d, surfIdx)
    ref = ONE
    @assertEqual(ref, d, TOL * d)
    
  end subroutine test_distance
  
end module meshUniverse_iTest