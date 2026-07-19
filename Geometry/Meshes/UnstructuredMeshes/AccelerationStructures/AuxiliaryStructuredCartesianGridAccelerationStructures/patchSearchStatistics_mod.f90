module patchSearchStatistics_mod
  !!
  !! Runtime mapping-trigger counters for the Patch-Search acceleration
  !! structure (revision experiment, R1.4). Module-level counters are used
  !! deliberately: findHostElementIdx takes self as intent(in), so the
  !! statistics cannot live on the object. Single-threaded use only -- if
  !! OpenMP is ever enabled, guard the increments with !$omp atomic.
  !!
  !! Usage: call resetPatchSearchStats() before a measurement run and
  !! reportPatchSearchStats() after it. Keep this on a dedicated stats
  !! branch/build -- never in a timing build.
  !!
  use numPrecision

  implicit none
  private

  integer(longInt), public, save :: nQueries            = 0
  integer(longInt), public, save :: nOutside            = 0
  integer(longInt), public, save :: nDirectElement      = 0  ! psi(T) > 0
  integer(longInt), public, save :: nSingleFace         = 0  ! single-face branch
  integer(longInt), public, save :: nAngularSearch      = 0  ! phi(T) > 0
  integer(longInt), public, save :: nVertexDisplacement = 0  ! varphi path (re-entrant)
  real(defReal), public, save    :: edgeMappingVolume = ZERO, elementMappingVolume = ZERO, outsideVolume = ZERO, &
                                    singleFaceVolume = ZERO, vertexMappingVolume = ZERO, totalVolume = ZERO

  public :: resetPatchSearchStats
  public :: reportPatchSearchStats

contains

  subroutine resetPatchSearchStats()
    nQueries            = 0
    nOutside            = 0
    nDirectElement      = 0
    nSingleFace         = 0
    nAngularSearch      = 0
    nVertexDisplacement = 0
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
    print '(A,I14)',          '  queries                 : ', nQueries
    print '(A,I14,F9.4,A)',   '  outside mesh            : ', nOutside, &
                              100.0_defReal * nOutside / q, ' %'
    print '(A,I14,F9.4,A)',   '  direct element (psi)    : ', nDirectElement, &
                              100.0_defReal * nDirectElement / q, ' %'
    print '(A,I14,F9.4,A)',   '  single-face branch      : ', nSingleFace, &
                              100.0_defReal * nSingleFace / q, ' %'
    print '(A,I14,F9.4,A)',   '  angular search (phi)    : ', nAngularSearch, &
                              100.0_defReal * nAngularSearch / q, ' %'
    print '(A,I14,F9.4,A)',   '  vertex displ. (varphi)  : ', nVertexDisplacement, &
                              100.0_defReal * nVertexDisplacement / q, ' %'
    print '(A,I14)',          '  total grid descents     : ', &
                              nQueries + nVertexDisplacement
    print '(A)', ''
  end subroutine reportPatchSearchStats

end module patchSearchStatistics_mod