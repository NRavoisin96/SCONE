module patchSearchStatistics_mod
  !!
  !! Runtime counters for the Patch-Search mapping-trigger statistics
  !! reported in Section 4.2 of the JCP paper. Enabled via the
  !! PATCH_SEARCH_STATS build option; the default build leaves the query
  !! path uninstrumented so that timing measurements are unaffected.
  !!
  !! Counters are module-level because findHostElementIdx takes self as
  !! intent(in). They are not thread-safe: the build refuses to configure
  !! with OpenMP enabled when PATCH_SEARCH_STATS is ON.
  !!
  !! Call resetPatchSearchStats() before a measurement run and
  !! reportPatchSearchStats() after it.
  !!
  use numPrecision

  implicit none
  private

  integer(longInt), public, save :: nAngularSearch = 0_longInt, nDirectElement = 0_longInt, nOutside = 0_longInt, &
                                    nQueries = 0_longInt, nSingleFace = 0_longInt, nVertexDisplacement = 0_longInt
  real(defReal), public, save    :: edgeMappingVolume = ZERO, elementMappingVolume = ZERO, outsideVolume = ZERO, &
                                    singleFaceVolume = ZERO, totalVolume = ZERO, vertexMappingVolume = ZERO

  public :: resetPatchSearchStats, reportPatchSearchStats

contains

  subroutine resetPatchSearchStats()
    nQueries = 0_longInt
    nOutside = 0_longInt
    nDirectElement = 0_longInt
    nSingleFace = 0_longInt
    nAngularSearch = 0_longInt
    nVertexDisplacement = 0_longInt
    edgeMappingVolume = ZERO
    elementMappingVolume = ZERO
    outsideVolume = ZERO
    singleFaceVolume = ZERO
    vertexMappingVolume = ZERO
    totalVolume = ZERO
    
  end subroutine resetPatchSearchStats

  subroutine reportPatchSearchStats()
    real(defReal) :: q

    q = real(max(1_longInt, nQueries), defReal)
    print '(A)', ''
    print '(A)', ' Patch-Search mapping trigger statistics'
    print '(A)', ' ---------------------------------------'
    print '(A,I14)',        'Queries                 : ', nQueries
    print '(A,I14,F9.4,A)', 'Outside mesh            : ', nOutside, 100.0_defReal * nOutside / q, ' %'
    print '(A,I14,F9.4,A)', 'Direct element (psi)    : ', nDirectElement, 100.0_defReal * nDirectElement / q, ' %'
    print '(A,I14,F9.4,A)', 'Single-face branch      : ', nSingleFace, 100.0_defReal * nSingleFace / q, ' %'
    print '(A,I14,F9.4,A)', 'Angular search (phi)    : ', nAngularSearch, 100.0_defReal * nAngularSearch / q, ' %'
    print '(A,I14,F9.4,A)', 'Vertex displ. (varphi)  : ', nVertexDisplacement, 100.0_defReal * nVertexDisplacement / q, ' %'
    print '(A,I14)',        'Total grid descents     : ', nQueries + nVertexDisplacement
    print '(A)', ''

  end subroutine reportPatchSearchStats

end module patchSearchStatistics_mod