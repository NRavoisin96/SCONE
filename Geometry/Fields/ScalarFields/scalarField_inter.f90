module scalarField_inter

  use coordList_class,    only : coordList
  use dictionary_class,   only : dictionary
  use errors_mod,         only : fatalError
  use field_inter,        only : field
  use geometryReg_mod,    only : fieldPtrByName
  use numPrecision
  use universalVariables, only : nameHeatSource, nameTemperature

  implicit none
  private

  ! Public procedures.
  public :: getHeatSourceFieldPtr, getTemperatureFieldPtr, kill, scalarField_CptrCast

  !!
  !! Simple Real Scalar Field
  !!
  !! Access to field is via coordList to allow more fancy fields to be defined
  !! (e.g. assign value to each uniqueID etc.)
  !!
  !! Interface:
  !!   field interface
  !!   at -> Return scalar value given position coordinates
  !!
  type, public, abstract, extends(field)     :: scalarField
    real(defReal), dimension(:), allocatable :: maximumMaterialValues, minimumMaterialValues
  contains
    procedure(at), deferred        :: at
    procedure                      :: getMaximumMaterialValue
    procedure                      :: getMinimumMaterialValue
    procedure                      :: getMaximumMaterialValues
    procedure                      :: getMinimumMaterialValues
    procedure                      :: kill
    procedure                      :: setMaximumMaterialValues
    procedure                      :: setMinimumMaterialValues
    procedure(setValues), deferred :: setValues
  end type scalarField

  abstract interface

    !!
    !! Get value of the scalar field at the co-ordinate point
    !!
    !! Args:
    !!   coords [in] -> Coordinates of the position in the geometry
    !!
    !! Result:
    !!   Value of the scalar field. Real number.
    !!
    function at(self, coords) result(val)
      import :: coordList, defReal, scalarField
      class(scalarField), intent(in) :: self
      class(coordList), intent(in)   :: coords
      real(defReal)                  :: val
    end function at

    !!
    !!
    !!
    subroutine setValues(self, values)
      import                                  :: defReal, scalarField
      class(scalarField), intent(inout)       :: self
      real(defReal), dimension(:), intent(in) :: values
    end subroutine setValues

  end interface

contains
  !!
  !!
  !!
  function getHeatSourceFieldPtr() result(heatSourceFieldPtr)
    class(field), pointer       :: fieldPtr
    class(scalarField), pointer :: heatSourceFieldPtr
    character(*), parameter     :: here = 'getHeatSourceFieldPtr (scalarField_inter.f90)'

    heatSourceFieldPtr => null()
    fieldPtr => fieldPtrByName(nameHeatSource)
    if (.not. associated(fieldPtr)) return

    select type(ptr => fieldPtr)
      class is(scalarField)
        heatSourceFieldPtr => ptr

      class default
        call fatalError(here, 'Heat source field is not of type scalarField.')

    end select

  end function getHeatSourceFieldPtr

  !!
  !!
  !!
  function getMaximumMaterialValue(self, materialIdx) result(maximumValue)
    class(scalarField), intent(in) :: self
    integer(shortInt), intent(in)  :: materialIdx
    real(defReal)                  :: maximumValue
    character(*), parameter        :: here = 'getMaximumMaterialValue (scalarField_inter.f90)'

    if (.not. allocated(self % maximumMaterialValues)) &
    call fatalError(here, 'Attempting to retrieve value from unallocated array.')

    maximumValue = self % maximumMaterialValues(materialIdx)

  end function getMaximumMaterialValue

  !!
  !!
  !!
  function getMinimumMaterialValue(self, materialIdx) result(minimumValue)
    class(scalarField), intent(in) :: self
    integer(shortInt), intent(in)  :: materialIdx
    real(defReal)                  :: minimumValue
    character(*), parameter        :: here = 'getMinimumMaterialValue (scalarField_inter.f90)'

    if (.not. allocated(self % minimumMaterialValues)) &
    call fatalError(here, 'Attempting to retrieve value from unallocated array.')

    minimumValue = self % minimumMaterialValues(materialIdx)

  end function getMinimumMaterialValue

  !!
  !!
  !!
  pure function getMaximumMaterialValues(self) result(maximumMaterialValues)
    class(scalarField), intent(in)           :: self
    real(defReal), dimension(:), allocatable :: maximumMaterialValues

    if (allocated(self % maximumMaterialValues)) then
      maximumMaterialValues = self % maximumMaterialValues

    else
      allocate(maximumMaterialValues(0))

    end if

  end function getMaximumMaterialValues

  !!
  !!
  !!
  pure function getMinimumMaterialValues(self) result(minimumMaterialValues)
    class(scalarField), intent(in)           :: self
    real(defReal), dimension(:), allocatable :: minimumMaterialValues

    if (allocated(self % minimumMaterialValues)) then
      minimumMaterialValues = self % minimumMaterialValues

    else
      allocate(minimumMaterialValues(0))

    end if

  end function getMinimumMaterialValues

  !!
  !!
  !!
  function getTemperatureFieldPtr() result(temperatureFieldPtr)
    class(field), pointer       :: fieldPtr
    class(scalarField), pointer :: temperatureFieldPtr
    character(*), parameter     :: here = 'getTemperatureFieldPtr (scalarField_inter.f90)'

    temperatureFieldPtr => null()
    fieldPtr => fieldPtrByName(nameTemperature)
    if (.not. associated(fieldPtr)) return

    select type(ptr => fieldPtr)
      class is(scalarField)
        temperatureFieldPtr => ptr

      class default
        call fatalError(here, 'Temperature field is not of type scalarField.')

    end select

  end function getTemperatureFieldPtr

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(scalarField), intent(inout) :: self

    ! Local.
    if (allocated(self % maximumMaterialValues)) deallocate(self % maximumMaterialValues)
    if (allocated(self % minimumMaterialValues)) deallocate(self % minimumMaterialValues)

  end subroutine kill

  !!
  !! Cast field pointer to scalarField pointer
  !!
  !! Args:
  !!   source [in] -> source pointer of class field
  !!
  !! Result:
  !!   Null is source is not of scalarField
  !!   Pointer to source if source is scalarField class
  !!
  pure function scalarField_CptrCast(source) result(ptr)
    class(field), pointer, intent(in) :: source
    class(scalarField), pointer       :: ptr

    select type (source)
      class is (scalarField)
        ptr => source

      class default
        ptr => null()
    end select

  end function scalarField_CptrCast

  !!
  !!
  !!
  pure subroutine setMaximumMaterialValues(self, maximumMaterialValues)
    class(scalarField), intent(inout)       :: self
    real(defReal), dimension(:), intent(in) :: maximumMaterialValues

    self % maximumMaterialValues = maximumMaterialValues

  end subroutine setMaximumMaterialValues

  !!
  !!
  !!
  pure subroutine setMinimumMaterialValues(self, minimumMaterialValues)
    class(scalarField), intent(inout)       :: self
    real(defReal), dimension(:), intent(in) :: minimumMaterialValues

    self % minimumMaterialValues = minimumMaterialValues

  end subroutine setMinimumMaterialValues

end module scalarField_inter