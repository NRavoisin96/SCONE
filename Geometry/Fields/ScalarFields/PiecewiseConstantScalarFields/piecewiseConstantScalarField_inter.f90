module piecewiseConstantScalarField_inter

  use genericProcedures, only : fatalError, numToChar
  use numPrecision
  use scalarField_inter, only : scalarField

  implicit none
  private

  ! Public procedures.
  public :: kill

  !!
  !!
  !!
  type, public, abstract, extends(scalarField) :: piecewiseConstantScalarField
    private
    integer(shortInt)                        :: nValues = 0
    real(defReal), dimension(:), allocatable :: values
  contains
    procedure :: allocateValues
    procedure :: getValue
    procedure :: getValuesNumber
    procedure :: kill
    procedure :: setValues
  end type piecewiseConstantScalarField

contains
  !!
  !!
  !!
  subroutine allocateValues(self, nValues)
    class(piecewiseConstantScalarField), intent(inout) :: self
    integer(shortInt), intent(in)                      :: nValues
    character(*), parameter                            :: here = 'allocateValues (piecewiseConstantScalarField_inter.f90)'

    if (allocated(self % values)) call fatalError(here, 'Attempting to allocate already allocated array.')
    allocate(self % values(nValues))
    self % nValues = nValues

  end subroutine allocateValues

  !!
  !!
  !!
  function getValue(self, idx) result(value)
    class(piecewiseConstantScalarField), intent(in) :: self
    integer(shortInt), intent(in)                   :: idx
    real(defReal)                                   :: value
    character(*), parameter                         :: here = 'getValue (piecewiseConstantScalarField_inter.f90)'

    if (.not. allocated(self % values)) call fatalError(here, 'Attempting to retrieve value of unallocated field.')
    if (idx < 1 .or. size(self % values) < idx) call fatalError(here, 'Index: '//numToChar(idx)//' is out of bounds.')
    value = self % values(idx)

  end function getValue

  !!
  !!
  !!
  function getValuesNumber(self) result(nValues)
    class(piecewiseConstantScalarField), intent(in) :: self
    integer(shortInt)                               :: nValues
    character(*), parameter                         :: here = 'getValuesNumber (piecewiseConstantScalarField_inter.f90)'

    if (.not. allocated(self % values)) call fatalError(here, 'Attempting to retrieve number of values of unallocated array.')
    nValues = self % nValues

  end function getValuesNumber

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(piecewiseConstantScalarField), intent(inout) :: self

    ! Local.
    self % nValues = 0
    if (allocated(self % values)) deallocate(self % values)

  end subroutine kill

  !!
  !!
  !!
  subroutine setValues(self, values)
    class(piecewiseConstantScalarField), intent(inout) :: self
    real(defReal), dimension(:), intent(in)            :: values
    character(*), parameter                            :: here = 'setValues (piecewiseConstantScalarField_inter.f90)'

    if (size(self % values) /= size(values)) call fatalError(here, 'Size mismatch.')
    self % values = values

  end subroutine setValues

end module piecewiseConstantScalarField_inter