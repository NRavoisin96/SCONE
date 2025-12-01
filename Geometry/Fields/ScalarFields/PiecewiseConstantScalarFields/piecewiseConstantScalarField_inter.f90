module piecewiseConstantScalarField_inter

  use genericProcedures,  only : fatalError, numToChar
  use intMap_class,       only : intMap
  use numPrecision
  use scalarField_inter,  only : kill_super => kill, scalarField
  use universalVariables, only : NOT_PRESENT

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
    real(defReal)                            :: maximumValue = ZERO, minimumValue = ZERO
    real(defReal), dimension(:), allocatable :: maximumMaterialValues, minimumMaterialValues, values
    type(intMap)                             :: materialIdxsToExtremalValuesIdxsMap
  contains
    procedure                                         :: allocateValues
    procedure(createExtremalMaterialValues), deferred :: createExtremalMaterialValues
    procedure                                         :: getMaximumValue
    procedure                                         :: getMinimumValue
    procedure                                         :: getValue
    procedure                                         :: getValuesNumber
    procedure                                         :: kill
    procedure                                         :: setValues
  end type piecewiseConstantScalarField

  abstract interface
    !!
    !!
    !!
    subroutine createExtremalMaterialValues(self, maximumMaterialValues, minimumMaterialValues, map)
      import                                                :: defReal, intMap, piecewiseConstantScalarField
      class(piecewiseConstantScalarField), intent(in)       :: self
      real(defReal), dimension(:), allocatable, intent(out) :: maximumMaterialValues, minimumMaterialValues
      type(intMap), intent(out)                             :: map
    end subroutine createExtremalMaterialValues

  end interface

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
  elemental function getMaximumValue(self, defaultValue, materialIdx, mult) result(maximumValue)
    class(piecewiseConstantScalarField), intent(in) :: self
    real(defReal), intent(in)                       :: defaultValue
    integer(shortInt), intent(in), optional         :: materialIdx
    real(defReal), intent(in), optional             :: mult
    integer(shortInt)                               :: mapIdx
    real(defReal)                                   :: maximumValue

    if (.not. allocated(self % maximumMaterialValues)) then
      maximumValue = defaultValue
      return

    end if

    if (present(materialIdx)) then
      mapIdx = self % materialIdxsToExtremalValuesIdxsMap % getOrDefault(materialIdx, NOT_PRESENT)
      if (mapIdx == NOT_PRESENT) then
        maximumValue = defaultValue
        return

      end if
      maximumValue = self % maximumMaterialValues(mapIdx)

    else
      maximumValue = self % maximumValue

    end if
    if (present(mult)) maximumValue = maximumValue * mult

  end function getMaximumValue

  !!
  !!
  !!
  elemental function getMinimumValue(self, defaultValue, materialIdx, mult) result(minimumValue)
    class(piecewiseConstantScalarField), intent(in) :: self
    real(defReal), intent(in)                       :: defaultValue
    integer(shortInt), intent(in), optional         :: materialIdx
    real(defReal), intent(in), optional             :: mult
    integer(shortInt)                               :: mapIdx
    real(defReal)                                   :: minimumValue

    if (.not. allocated(self % minimumMaterialValues)) then
      minimumValue = defaultValue
      return

    end if

    if (present(materialIdx)) then
      mapIdx = self % materialIdxsToExtremalValuesIdxsMap % getOrDefault(materialIdx, NOT_PRESENT)
      if (mapIdx == NOT_PRESENT) then
        minimumValue = defaultValue
        return

      end if
      minimumValue = self % minimumMaterialValues(mapIdx)

    else
      minimumValue = self % minimumValue

    end if
    if (present(mult)) minimumValue = minimumValue * mult

  end function getMinimumValue

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

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % nValues = 0
    self % maximumValue = ZERO
    self % minimumValue = ZERO
    if (allocated(self % values)) deallocate(self % values)
    if (allocated(self % maximumMaterialValues)) deallocate(self % maximumMaterialValues)
    if (allocated(self % minimumMaterialValues)) deallocate(self % minimumMaterialValues)
    call self % materialIdxsToExtremalValuesIdxsMap % kill()

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
    call self % createExtremalMaterialValues(self % maximumMaterialValues, self % minimumMaterialValues, &
                                             self % materialIdxsToExtremalValuesIdxsMap)
    self % maximumValue = maxval(self % maximumMaterialValues)
    self % minimumValue = minval(self % minimumMaterialValues)

  end subroutine setValues

end module piecewiseConstantScalarField_inter