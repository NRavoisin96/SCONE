module quadSurface_inter

  use numPrecision
  use universalVariables, only : INF, ZERO
  use genericProcedures,  only : solveQuadratic
  use surface_inter,      only : surface

  implicit none
  private

  !!
  !! Abstract interface to group all quadratic surfaces.
  !!
  !! At the moment it adds nothing new, except grouping the quadratic
  !! surfaces into a familly of subclasses of surface.
  !!
  type, public, abstract, extends(surface) :: quadSurface
    private
  contains
    procedure, non_overridable             :: distanceQuad
  end type quadSurface

contains

  !! Subroutine 'distanceQuad'
  !!
  !! Basic description:
  !!   Computes the distance to intersection to the quadratic surface given by the solutions to the
  !!   quadratic equation a * d^2 + 2 * k * d + c = 0.
  !!
  !! Arguments:
  !!   a [in]                         -> Coefficient of d^2 in a * d^2 + 2 * k * d + c = 0.
  !!   k [in]                         -> Coefficient of d in a * d^2 + 2 * k * d + c = 0.
  !!   c [in]                         -> Coefficient c in a * d^2 + 2 * k * d + c = 0.
  !!   surfTolCondition [in]          -> .true. if particle is within surface tolerance of the 
  !!                                     quadratic surface.
  !!
  elemental function distanceQuad(self, a, k, c, surfTolCondition) result(d)
    class(quadSurface), intent(in)           :: self
    real(defReal), intent(in)                :: a, k, c
    logical(defBool), intent(in)             :: surfTolCondition
    real(defReal)                            :: d, delta, dMax
    real(defReal), dimension(:), allocatable :: solutions
    logical(defBool)                         :: noIntersection
    
    ! Initialise d = INF and compute delta (technically, delta / 4).
    d = INF
    delta = k * k - a * c

    if (delta < ZERO) return
    solutions = solveQuadratic(a, k, c, delta)
    if (size(solutions) == 0) return

    ! Check if particle is within surface tolerance of the surface.
    if (surfTolCondition) then
      ! Update d only if k < ZERO (k >= ZERO corresponds to the particle moving on or away from the surface). 
      ! Choose maximum distance and return.
      dMax = maxval(solutions)
      if (k < ZERO .and. dMax >= self % getSurfTol()) d = dMax
      return

    end if

    ! If reached here, update d depending on the sign of c and cap distance at infinity.
    d = minval(solutions, solutions > ZERO)
    if (d <= ZERO .or. d > INF) d = INF

  end function distanceQuad
  
end module quadSurface_inter