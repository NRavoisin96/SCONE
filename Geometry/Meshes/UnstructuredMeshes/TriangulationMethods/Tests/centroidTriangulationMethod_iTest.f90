module centroidTriangulationMethod_iTest

  !! Module 'centroidTriangulationMethod_iTest'
  !!
  !! Basic description:
  !!   Integration tests for the centroid-based triangulation method.
  !!
  !! Detailed description:
  !!   A tetrahedralised mesh describes the same geometry as the mesh it was built from
  !!   and only the decomposition differs, so the physically meaningful outputs of
  !!   tracking are identical between the two. The tests assert that rather than any
  !!   hand-computed tetrahedron indices:
  !!
  !!     d       -> distance to entry must agree.
  !!     localId -> tetrahedra inherit their parent element's local ID, so the local ID
  !!                must agree even though the element indices do not.
  !!
  !!   elementIdx is deliberately not compared. The two meshes number their cells
  !!   differently and the tetrahedron indices carry no independent meaning. Asserting
  !!   them was the previous design and it produced a long table of numbers that had to
  !!   be recomputed by hand whenever the decomposition changed, with no way to tell a
  !!   wrong expectation from a wrong answer.
  !!
  !!   Using the undecomposed mesh as the oracle is also stronger than a fixed table.
  !!   Any defect specific to the tetrahedralised path, such as missing rational face
  !!   normals or missing coordinates on the manufactured centroid vertices, shows up as
  !!   a disagreement instead of a plausible-looking number.
  !!
  !!   NOT COMPARED:
  !!     A single tracking step. One crossing means different things in the two meshes,
  !!     since the decomposition adds internal faces: after one step the undecomposed
  !!     mesh may have left the mesh while the tetrahedralised one is still inside it.
  !!     Neither the step distance nor the state after one step is comparable. Comparing
  !!     the total distance travelled through the mesh would be valid but needs the real
  !!     transport loop rather than a hand-rolled walk, and is left as a follow-up.
  !!
  !!   Rays are built from exactly representable coordinates and integer direction
  !!   patterns. See OpenFOAMMesh_iTest for why: decimals such as 0.2 do not round the
  !!   same way as 1.2, and directions with unequal component magnitudes cannot place a
  !!   crossing exactly on a feature.
  !!

  use charMap_class,      only : charMap
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use funit
  use numPrecision
  use OpenFOAMMesh_class, only : OpenFOAMMesh
  use publicObjects,      only : coordData, newCoordData
  use universalVariables

  implicit none

  ! Same mesh, imported with and without the triangulation method.
  character(*), parameter :: MESH_DEF_TET = &
  " id 2; type OpenFOAMMesh; path ./IntegrationTestFiles/Geometry/Meshes/OpenFOAM/testMesh/; triangulationMethod centroidBased;&
  & fills (water);"

  character(*), parameter :: MESH_DEF_PLAIN = &
  " id 3; type OpenFOAMMesh; path ./IntegrationTestFiles/Geometry/Meshes/OpenFOAM/testMesh/; fills (water);"

  ! Geometry constants. All exactly representable. The mesh spans [-1, 1] in x and
  ! y with an internal plane at x = 0 and y = 0, and [-1, 1] in z as a single layer.
  real(defReal), parameter :: EXTENT = ONE, HALF_CELL = HALF, QUARTER_CELL = FOURTH
  real(defReal), parameter :: FAR = TWO, TOL = 1.0E-6_defReal

  type(charMap)      :: mats, matsPlain
  type(OpenFOAMMesh) :: mesh, meshPlain

contains

  !!
  !! Import the mesh twice, once tetrahedralised and once as it comes.
  !!
@Before
  subroutine setUp()
    type(dictionary)   :: dict, dictPlain
    character(nameLen) :: name
    character(pathLen) :: path

    name = 'water'
    call mats % add(name, 1)
    call charToDict(dict, MESH_DEF_TET)
    call dict % get(path, 'path')
    call mesh % init(trim(path), dict, mats)

    call matsPlain % add(name, 1)
    call charToDict(dictPlain, MESH_DEF_PLAIN)
    call dictPlain % get(path, 'path')
    call meshPlain % init(trim(path), dictPlain, matsPlain)

  end subroutine setUp

  !!
  !! Clean after tests.
  !!
@After
  subroutine cleanUp()

    call mats % kill()
    call matsPlain % kill()
    call mesh % kill()
    call meshPlain % kill()

  end subroutine cleanUp

  !! ------------------------------------------------------------------------------
  !! Helpers
  !! ------------------------------------------------------------------------------

  !! Function 'newRay'
  !!
  !! Basic description:
  !!   Builds a ray starting at r and travelling along the integer direction pattern
  !!   uPattern for at most dMax.
  !!
  !! Arguments:
  !!   r [in]        -> Starting coordinates of the ray.
  !!   uPattern [in] -> Direction pattern. Components must be in {-1, 0, 1}.
  !!   dMax [in]     -> Optional. Maximum distance travelled. Defaults to FAR.
  !!
  !! Result:
  !!   data          -> Coordinate data for the ray.
  !!
  function newRay(r, uPattern, dMax) result(data)
    real(defReal), dimension(3), intent(in)     :: r
    integer(shortInt), dimension(3), intent(in) :: uPattern
    real(defReal), intent(in), optional         :: dMax
    real(defReal)                               :: maximumDistance
    type(coordData)                             :: data

    maximumDistance = FAR
    if (present(dMax)) maximumDistance = dMax

    data = newCoordData(r, real(uPattern, defReal), dMax = maximumDistance)

  end function newRay

  !! Subroutine 'reportDisagreement'
  !!
  !! Basic description:
  !!   Prints the case name when the two meshes disagree.
  !!
  !! Detailed description:
  !!   This pFUnit build takes no message argument on the assertions used here, and both
  !!   helpers are called from several tests, so without this a failure only gives a line
  !!   number inside the helper.
  !!
  !! Arguments:
  !!   dataPlain [in] -> Coordinate data obtained from the undecomposed mesh.
  !!   dataTet [in]   -> Coordinate data obtained from the tetrahedralised mesh.
  !!   message [in]   -> Name of the case, printed on disagreement.
  !!
  subroutine reportDisagreement(dataPlain, dataTet, message)
    type(coordData), intent(in) :: dataPlain, dataTet
    character(*), intent(in)    :: message
    logical(defBool)            :: disagrees

    disagrees = (dataPlain % elementIdx == 0) .neqv. (dataTet % elementIdx == 0)
    if (dataPlain % localId /= dataTet % localId) disagrees = .true.

    if (disagrees) then
      print *, 'DISAGREEMENT: '//message
      print *, '  undecomposed: d =', dataPlain % d, ' elementIdx =', dataPlain % elementIdx, &
               ' localId =', dataPlain % localId
      print *, '  tetrahedral : d =', dataTet % d, ' elementIdx =', dataTet % elementIdx, &
               ' localId =', dataTet % localId

    end if

  end subroutine reportDisagreement

  !! Subroutine 'assertEntryAgrees'
  !!
  !! Basic description:
  !!   Enters the mesh from outside along the same ray in both meshes and checks that the
  !!   decomposition made no difference.
  !!
  !! Arguments:
  !!   r [in]        -> Starting coordinates of the ray.
  !!   uPattern [in] -> Direction pattern. Components must be in {-1, 0, 1}.
  !!   message [in]  -> Name of the case, printed on disagreement.
  !!
  subroutine assertEntryAgrees(r, uPattern, message)
    real(defReal), dimension(3), intent(in)     :: r
    integer(shortInt), dimension(3), intent(in) :: uPattern
    character(*), intent(in)                    :: message
    type(coordData)                             :: dataPlain, dataTet

    dataPlain = newRay(r, uPattern)
    call meshPlain % distanceToBoundary(dataPlain)

    dataTet = newRay(r, uPattern)
    call mesh % distanceToBoundary(dataTet)

    call reportDisagreement(dataPlain, dataTet, message)

    @assertTrue((dataPlain % elementIdx == 0) .eqv. (dataTet % elementIdx == 0))
    @assertEqual(dataPlain % d, dataTet % d, max(abs(dataPlain % d), ONE) * TOL)
    @assertEqual(dataPlain % localId, dataTet % localId)

  end subroutine assertEntryAgrees

  !! ------------------------------------------------------------------------------
  !! Mesh information
  !! ------------------------------------------------------------------------------

  !!
  !! Test miscellaneous functionality.
  !!
@Test
  subroutine test_misc()

    @assertEqual(2, mesh % getId())
    call mesh % setId(7)
    @assertEqual(7, mesh % getId())

  end subroutine test_misc

  !!
  !! Test mesh information. These are properties of the decomposition itself, so this is
  !! the one place where fixed numbers belong.
  !!
@Test
  subroutine test_info()

    @assertEqual(22, mesh % getVerticesNumber())
    @assertEqual(112, mesh % getFacesNumber(.true.))
    @assertEqual(80, mesh % getInternalFacesNumber(.true.))
    @assertEqual(85, mesh % getEdgesNumber())
    @assertEqual(48, mesh % getElementsNumber(.true.))

  end subroutine test_info

  !!
  !! The tetrahedralisation must cover exactly the same region as the mesh it came from.
  !! A point is inside one if and only if it is inside the other.
  !!
@Test
  subroutine test_inside_agrees_with_undecomposed_mesh()
    integer(shortInt)           :: i, j, k
    real(defReal), dimension(3) :: r
    type(coordData)             :: dataPlain, dataTet

    do i = -3, 3
      do j = -3, 3
        do k = -3, 3
          ! Sweep a lattice of exactly representable points covering the mesh and its
          ! surroundings, including points on the boundary and on the internal planes.
          r = HALF_CELL * real([i, j, k], defReal)

          dataPlain = newRay(r, [1, 0, 0])
          call meshPlain % findHostElement(dataPlain)

          dataTet = newRay(r, [1, 0, 0])
          call mesh % findHostElement(dataTet)

          @assertTrue((dataPlain % elementIdx == 0) .eqv. (dataTet % elementIdx == 0))
          if (0 < dataPlain % elementIdx) then
            @assertEqual(dataPlain % localId, dataTet % localId)

          end if

        end do

      end do

    end do

  end subroutine test_inside_agrees_with_undecomposed_mesh

  !! ------------------------------------------------------------------------------
  !! Entry into the mesh from outside
  !! ------------------------------------------------------------------------------

  !!
  !! Entry through the interior of a boundary face.
  !!
@Test
  subroutine test_entry_through_face_interiors()

    call assertEntryAgrees([-TWO, -HALF_CELL, ZERO], [1, 0, 0], 'entering along +x')
    call assertEntryAgrees([HALF_CELL, TWO, ZERO], [0, -1, 0], 'entering along -y')

  end subroutine test_entry_through_face_interiors

  !!
  !! Entry aimed exactly at boundary edges and corners.
  !!
@Test
  subroutine test_entry_through_features()

    ! Aimed at the boundary edge x = -1, y = 0.
    call assertEntryAgrees([-TWO, -EXTENT, ZERO], [1, 1, 0], 'entering at boundary edge')

    ! Aimed at the boundary corner (-1, -1, -1).
    call assertEntryAgrees([-TWO, -TWO, -TWO], [1, 1, 1], 'entering at boundary corner')

    ! Aimed at the centre of the mesh from outside, through a boundary face interior and
    ! on towards the internal features.
    call assertEntryAgrees([-TWO, -TWO, ZERO], [1, 1, 0], 'entering towards mesh centre')

  end subroutine test_entry_through_features

  !!
  !! Rays that do not enter the mesh must not enter either version of it.
  !!
@Test
  subroutine test_no_entry_agrees()

    ! Points away from the mesh.
    call assertEntryAgrees([-TWO, ZERO, ZERO], [-1, 0, 0], 'pointing away')

    ! Passes the mesh entirely.
    call assertEntryAgrees([-TWO, -TWO, ZERO], [1, 0, 0], 'missing the mesh')

    ! Grazes the boundary corner from outside without entering.
    call assertEntryAgrees([-TWO, -TWO, -TWO], [1, -1, 0], 'grazing the corner')

  end subroutine test_no_entry_agrees

end module centroidTriangulationMethod_iTest