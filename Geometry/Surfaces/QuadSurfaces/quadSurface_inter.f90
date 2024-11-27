module quadSurface_inter

  use numPrecision
  use universalVariables, only : INF, ZERO
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
    procedure, non_overridable :: distanceQuad
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
  !!   d [inout]                      -> Distance to intersection with the quadratic surface.
  !!   noIntersection [out]           -> .true. if there is no intersection with the quadratic surface.
  !!   extraCondition [in] [optional] -> Extra logical condition for which there is no intersection with
  !!                                     the quadratic surface when .true. (e.g., a = ZERO for a 
  !!                                     cylinder).
  !!
  elemental subroutine distanceQuad(self, a, k, c, surfTolCondition, d, noIntersection, extraCondition)
    class(quadSurface), intent(in)         :: self
    real(defReal), intent(in)              :: a, k, c
    logical(defBool), intent(in)           :: surfTolCondition
    real(defReal), intent(inout)           :: d
    logical(defBool), intent(out)          :: noIntersection
    logical(defBool), intent(in), optional :: extraCondition
    real(defReal)                          :: delta
    
    ! Compute delta (technically, delta / 4).
    delta = k * k - a * c

    ! If delta < ZERO, the solutions are complex. In this case, or if any other conditions are present
    ! (e.g., a = ZERO for a cylinder), there are no intersections and we can return early.
    noIntersection = delta < ZERO
    if (present(extraCondition)) noIntersection = noIntersection .or. extraCondition
    if (noIntersection) return

    ! Check if particle is within surface tolerance of the surface.
    if (surfTolCondition) then
      ! Update d only if k < ZERO (k >= ZERO corresponds to the particle moving on or away from the surface). 
      ! Choose maximum distance and return.
      if (k < ZERO) d = (-k + sqrt(delta)) / a
      return

    end if

    ! If reached here, update d depending on the sign of c and cap distance at infinity.
    d = -(k + sign(sqrt(delta), c)) / a
    if (d <= ZERO .or. d > INF) d = INF

  end subroutine distanceQuad
  
end module quadSurface_inter