module OpenFOAMMesh_iTest

  !! Module 'OpenFOAMMeshExact_iTest'
  !!
  !! Basic description:
  !!   Integration tests for the exact (rational) tracking path on the 8-cube mesh.
  !!
  !! Detailed description:
  !!   A case only tests what it says it tests if the crossing point lands exactly on
  !!   the feature it names. Two rules make that true by construction and both must be
  !!   respected when adding cases:
  !!
  !!   1. Every coordinate is an exact binary fraction (0.5, 0.25, 0.125, ...). Decimals
  !!      such as 0.2 or 1.2 are not exactly representable and do not round the same
  !!      way, so a case built from them misses the intended feature by ~1e-17 and ends
  !!      up testing something else entirely.
  !!
  !!   2. Direction components are taken from {-1, 0, 1} only. Once normalised the
  !!      non-zero components are the same double up to sign, so u_i / u_j is exactly
  !!      +-1 and the crossing point p_j = r_j - (u_j / u_i) * r_i is an exact
  !!      difference of binary fractions. A direction such as (6, -4, 5) has no exact
  !!      component ratio and can never land on a feature.
  !!
  !!   NOTE: the end point r + u * dMax is never exact, even under these rules, since
  !!         u = v / norm2(v) followed by u * norm2(v) does not give v back. No case
  !!         relies on the end point landing on a feature. dMax is either comfortably
  !!         large or set to a deliberate multiple of a known crossing distance.
  !!
  !! Mesh:
  !!   ./IntegrationTestFiles/Geometry/Meshes/OpenFOAM/cubes/ -> eight cubes spanning
  !!   [-0.5, 0.5]^3 with internal planes at x = y = z = 0 and cell centres at
  !!   (+-0.25, +-0.25, +-0.25). Cells are numbered by octant:
  !!
  !!     1 = (-,-,-)   2 = (+,-,-)   3 = (-,+,-)   4 = (+,+,-)
  !!     5 = (-,-,+)   6 = (+,-,+)   7 = (-,+,+)   8 = (+,+,+)
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

  character(*), parameter :: MESH_DEF = &
  " id 2; type OpenFOAMMesh; path ./IntegrationTestFiles/Geometry/Meshes/OpenFOAM/cubes/; fills (fuel);"

  ! Geometry constants. All exactly representable.
  real(defReal), parameter :: HALF_EXTENT = HALF, EIGHTH = 0.125_defReal
  real(defReal), parameter :: FAR = TWO, TOL = 1.0E-6_defReal

  ! Cell indices ordered by octant, flattened as in the mesh file.
  integer(shortInt), dimension(8), parameter :: OCTANT_CELL = [1, 2, 3, 4, 5, 6, 7, 8]

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

  !! ------------------------------------------------------------------------------
  !! Helpers
  !! ------------------------------------------------------------------------------

  !! Function 'newRay'
  !!
  !! Basic description:
  !!   Builds a ray starting at r and travelling along the integer direction pattern
  !!   uPattern for at most dMax. See the module header for why the direction pattern
  !!   is restricted to components in {-1, 0, 1}.
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
    real(defReal), dimension(3), intent(in) :: r
    integer(shortInt), dimension(3), intent(in) :: uPattern
    real(defReal), intent(in), optional      :: dMax
    real(defReal)                            :: maximumDistance
    type(coordData)                          :: data

    maximumDistance = FAR
    if (present(dMax)) maximumDistance = dMax

    data = newCoordData(r, real(uPattern, defReal), dMax = maximumDistance)

  end function newRay

  !! Function 'getOctantCellIdx'
  !!
  !! Basic description:
  !!   Returns the index of the cell occupying a given octant.
  !!
  !! Arguments:
  !!   octant [in] -> Octant of the cell. Each component is +1 on the positive side of
  !!                  the corresponding internal plane and -1 on the negative side.
  !!
  !! Result:
  !!   cellIdx     -> Index of the cell.
  !!
  pure function getOctantCellIdx(octant) result(cellIdx)
    integer(shortInt), dimension(3), intent(in) :: octant
    integer(shortInt)                           :: cellIdx, i, j, k

    i = merge(1, 0, octant(1) > 0)
    j = merge(1, 0, octant(2) > 0)
    k = merge(1, 0, octant(3) > 0)
    cellIdx = OCTANT_CELL(1 + i + 2 * j + 4 * k)

  end function getOctantCellIdx

  !! ------------------------------------------------------------------------------
  !! Crossings from inside the mesh
  !! ------------------------------------------------------------------------------

  !!
  !! Crossing through the interior of an internal face. No feature is involved here, so
  !! the fast path should resolve it on its own.
  !!
@Test
  subroutine test_crossing_through_face_interior()
    type(coordData) :: data

    ! From the centre of cell 1 towards +x. Crosses x = 0 at (0, -0.25, -0.25), well
    ! inside the face, and enters cell 2.
    data = newRay([-FOURTH, -FOURTH, -FOURTH], [1, 0, 0])
    call mesh % findHostElement(data)
    @assertEqual(getOctantCellIdx([-1, -1, -1]), data % elementIdx)

    call mesh % distanceToNextFace(data)
    @assertEqual(FOURTH, data % d, FOURTH * TOL)
    @assertEqual(getOctantCellIdx([1, -1, -1]), data % elementIdx)
    @assertEqual(1, data % front)

  end subroutine test_crossing_through_face_interior

  !!
  !! Crossing exactly along an internal edge. Two faces are crossed at the same distance
  !! and the host is the cell the direction points into.
  !!
@Test
  subroutine test_crossing_through_internal_edge()
    type(coordData) :: data

    ! From the centre of cell 1 towards (+x, +y). Crosses x = 0 and y = 0 together at
    ! (0, 0, -0.25), on the internal edge, and enters cell 4.
    data = newRay([-FOURTH, -FOURTH, -FOURTH], [1, 1, 0])
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)

    @assertEqual(FOURTH * sqrt(TWO), data % d, FOURTH * sqrt(TWO) * TOL)
    @assertEqual(2, data % front)
    @assertEqual(getOctantCellIdx([1, 1, -1]), data % elementIdx)

  end subroutine test_crossing_through_internal_edge

  !!
  !! Crossing exactly through the centre vertex. Three faces are crossed at the same
  !! distance and the host is the diagonally opposite cell.
  !!
@Test
  subroutine test_crossing_through_centre_vertex()
    type(coordData) :: data

    ! From the centre of cell 1 towards (+x, +y, +z). Crosses all three internal planes
    ! together at the origin and enters cell 8.
    data = newRay([-FOURTH, -FOURTH, -FOURTH], [1, 1, 1])
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)

    @assertEqual(FOURTH * sqrt(3.0_defReal), data % d, FOURTH * sqrt(3.0_defReal) * TOL)
    @assertEqual(3, data % front)
    @assertEqual(getOctantCellIdx([1, 1, 1]), data % elementIdx)

  end subroutine test_crossing_through_centre_vertex

  !!
  !! At the centre vertex the host is decided by direction alone. Sweeping the eight
  !! diagonal directions from the eight cells must land in the diagonally opposite cell
  !! every time.
  !!
@Test
  subroutine test_centre_vertex_host_all_directions()
    integer(shortInt)               :: i, j, k
    integer(shortInt), dimension(3) :: octant
    type(coordData)                 :: data

    do i = -1, 1, 2
      do j = -1, 1, 2
        do k = -1, 1, 2
          octant = [i, j, k]

          ! Start at the centre of the cell in this octant and head for the origin.
          data = newRay(FOURTH * real(octant, defReal), -octant)
          call mesh % findHostElement(data)
          @assertEqual(getOctantCellIdx(octant), data % elementIdx)

          call mesh % distanceToNextFace(data)
          @assertEqual(3, data % front)
          @assertEqual(getOctantCellIdx(-octant), data % elementIdx)

        end do

      end do

    end do

  end subroutine test_centre_vertex_host_all_directions

  !! ------------------------------------------------------------------------------
  !! Stopping short of, exactly on, and past a feature
  !! ------------------------------------------------------------------------------

  !!
  !! A segment stopping short of the centre vertex does not cross it. dMax is taken from
  !! the exact crossing distance so that the comparison is unambiguous.
  !!
@Test
  subroutine test_stopping_short_of_centre_vertex()
    real(defReal)   :: crossingDistance
    type(coordData) :: data

    crossingDistance = FOURTH * sqrt(3.0_defReal)

    ! Stop halfway to the vertex. No crossing, host unchanged.
    data = newRay([-FOURTH, -FOURTH, -FOURTH], [1, 1, 1], dMax = HALF * crossingDistance)
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)

    @assertEqual(INF, data % d)
    @assertEqual(getOctantCellIdx([-1, -1, -1]), data % elementIdx)

  end subroutine test_stopping_short_of_centre_vertex

  !!
  !! A segment reaching well past the centre vertex crosses it, and the distance
  !! returned is the distance to the vertex rather than to the end of the segment.
  !!
@Test
  subroutine test_passing_beyond_centre_vertex()
    real(defReal)   :: crossingDistance
    type(coordData) :: data

    crossingDistance = FOURTH * sqrt(3.0_defReal)

    data = newRay([-FOURTH, -FOURTH, -FOURTH], [1, 1, 1], dMax = TWO * crossingDistance)
    call mesh % findHostElement(data)
    call mesh % distanceToNextFace(data)

    @assertEqual(crossingDistance, data % d, crossingDistance * TOL)
    @assertEqual(getOctantCellIdx([1, 1, 1]), data % elementIdx)

  end subroutine test_passing_beyond_centre_vertex

  !! ------------------------------------------------------------------------------
  !! Sideways offsets near the centre vertex
  !! ------------------------------------------------------------------------------

  !!
  !! Sweeps a lateral offset applied to the start point. The offset only exists if it
  !! survives the rounding of (-0.25 - offset): below half an ulp it vanishes and the ray
  !! passes exactly through the vertex, above it the ray passes beside the vertex and
  !! crosses an internal edge instead.
  !!
  !! NOTE: the expectation is derived from whether the offset survives rather than
  !!       tabulated, so the test states the rule instead of its consequences.
  !!
@Test
  subroutine test_sideways_offsets_near_centre_vertex()
    integer(shortInt)                        :: i
    logical(defBool)                         :: offsetSurvives
    real(defReal)                            :: offset, startCoordinate
    real(defReal), dimension(6), parameter   :: OFFSETS = [1.0E-9_defReal, 1.0E-12_defReal, 1.0E-14_defReal, &
                                                           1.0E-16_defReal, 1.0E-17_defReal, 1.0E-20_defReal]
    type(coordData)                          :: data

    do i = 1, size(OFFSETS)
      offset = OFFSETS(i)
      startCoordinate = -FOURTH - offset
      offsetSurvives = startCoordinate /= -FOURTH

      ! Offset applied to z only. With the offset present the ray reaches x = 0 and
      ! y = 0 together while z = 0 is still ahead, so an internal edge is crossed. Once
      ! the offset rounds away all three planes coincide.
      data = newRay([-FOURTH, -FOURTH, startCoordinate], [1, 1, 1])
      call mesh % findHostElement(data)
      call mesh % distanceToNextFace(data)

      if (offsetSurvives) then
        ! Edge crossing into the cell below the diagonal one.
        @assertEqual(2, data % front)
        @assertEqual(getOctantCellIdx([1, 1, -1]), data % elementIdx)

      else
        ! Offset rounded away, so this is the exact vertex crossing.
        @assertEqual(3, data % front)
        @assertEqual(getOctantCellIdx([1, 1, 1]), data % elementIdx)

      end if

    end do

  end subroutine test_sideways_offsets_near_centre_vertex

  !! ------------------------------------------------------------------------------
  !! Entry into the mesh from outside
  !! ------------------------------------------------------------------------------

  !!
  !! Entry through the interior of a boundary face.
  !!
@Test
  subroutine test_entry_through_boundary_face()
    type(coordData) :: data

    ! Approaches the x = -0.5 face along +x and enters the interior of cell 1's boundary
    ! face at (-0.5, -0.25, -0.25).
    data = newRay([-ONE, -FOURTH, -FOURTH], [1, 0, 0])
    call mesh % distanceToBoundary(data)

    @assertEqual(HALF_EXTENT, data % d, HALF_EXTENT * TOL)
    @assertEqual(getOctantCellIdx([-1, -1, -1]), data % elementIdx)
    @assertEqual(1, data % front)

  end subroutine test_entry_through_boundary_face

  !!
  !! Entry exactly through a boundary edge. Two boundary faces are hit at the same
  !! distance and the direction decides which cell is entered.
  !!
@Test
  subroutine test_entry_through_boundary_edge()
    type(coordData) :: data

    ! Aims at the boundary edge x = -0.5, y = 0, arriving at (-0.5, 0, -0.25). The +y
    ! component means cell 3's side is entered.
    data = newRay([-ONE, -HALF_EXTENT, -FOURTH], [1, 1, 0])
    call mesh % distanceToBoundary(data)

    @assertEqual(2, data % front)
    @assertEqual(getOctantCellIdx([-1, 1, -1]), data % elementIdx)

  end subroutine test_entry_through_boundary_edge

  !!
  !! Entry exactly through a boundary vertex. Three boundary faces are hit together.
  !!
@Test
  subroutine test_entry_through_boundary_vertex()
    type(coordData) :: data

    ! Aims at the corner (-0.5, -0.5, -0.5) from outside along (+1, +1, +1).
    data = newRay([-ONE, -ONE, -ONE], [1, 1, 1])
    call mesh % distanceToBoundary(data)

    @assertEqual(3, data % front)
    @assertEqual(getOctantCellIdx([-1, -1, -1]), data % elementIdx)

  end subroutine test_entry_through_boundary_vertex

  !!
  !! Rays reaching the mesh boundary but travelling away from it, and rays missing the
  !! mesh altogether, must report no entry.
  !!
@Test
  subroutine test_no_entry()
    type(coordData) :: data

    ! Points directly away from the mesh.
    data = newRay([-ONE, -FOURTH, -FOURTH], [-1, 0, 0])
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)

    ! Passes the mesh entirely, offset in y beyond its extent.
    data = newRay([-ONE, -ONE, -FOURTH], [1, 0, 0])
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)

    ! Aimed at the mesh but stops short of it.
    data = newRay([-ONE, -FOURTH, -FOURTH], [1, 0, 0], dMax = FOURTH)
    call mesh % distanceToBoundary(data)
    @assertEqual(0, data % elementIdx)
    @assertEqual(INF, data % d)

  end subroutine test_no_entry

  !! ------------------------------------------------------------------------------
  !! Invariants
  !! ------------------------------------------------------------------------------

  !!
  !! Walks a ray across the mesh one crossing at a time and checks the invariants that
  !! must hold at every step, without any hand-computed expectations: the distance is
  !! positive and finite until the mesh is left, the number of faces reported is between
  !! 1 and 3, and the particle ends up outside after a bounded number of crossings.
  !!
@Test
  subroutine test_traversal_invariants()
    integer(shortInt)           :: nCrossings
    real(defReal), dimension(3) :: u
    type(coordData)             :: data

    data = newRay([-ONE, -EIGHTH, EIGHTH], [1, 0, 0])
    call mesh % distanceToBoundary(data)
    @assertTrue(0 < data % elementIdx)

    u = data % u
    nCrossings = 0
    do while (0 < data % elementIdx .and. nCrossings < 10)
      ! Move to the crossing point and take the next step from there.
      data = newCoordData(data % r + data % d * u, u, dMax = FAR)
      call mesh % findHostElement(data)
      if (data % elementIdx == 0) exit

      call mesh % distanceToNextFace(data)
      nCrossings = nCrossings + 1

      if (0 < data % elementIdx) then
        @assertTrue(ZERO < data % d)
        @assertTrue(data % d < INF)
        @assertTrue(1 <= data % front .and. data % front <= 3)

      end if

    end do

    ! A straight ray through a 2 x 2 x 2 block crosses at most two internal planes.
    @assertTrue(nCrossings <= 3)

  end subroutine test_traversal_invariants

end module OpenFOAMMesh_iTest