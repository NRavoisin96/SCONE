module scalarField_inter

  use coordList_class,    only : coordList
  use dictionary_class,   only : dictionary
  use errors_mod,         only : fatalError
  use field_inter,        only : field
  use geometryReg_mod,    only : fieldPtrByName
  use numPrecision

  implicit none
  private

  ! Public procedures.
  public :: getMaximumScalarFieldValue, getMinimumScalarFieldValue, getScalarFieldValue, kill, castScalarFieldPtr

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
  type, public, abstract, extends(field) :: scalarField
    private
  contains
    procedure(at), deferred              :: at
    procedure(getMaximumValue), deferred :: getMaximumValue
    procedure(getMinimumValue), deferred :: getMinimumValue
    procedure                            :: kill
    procedure(setValues), deferred       :: setValues
  end type scalarField

  abstract interface

    !!
    !! Get value of the scalar field at the co-ordinate point
    !!
    !! Args:
    !!   defaultValue [in] -> Value to use in case the field does not cover the coordinates.
    !!   coords [in]       -> Coordinates of the position in the geometry
    !!
    !! Result:
    !!   Value of the scalar field. Real number.
    !!
    function at(self, defaultValue, coords, mult) result(val)
      import                              :: coordList, defReal, scalarField
      class(scalarField), intent(in)      :: self
      real(defReal), intent(in)           :: defaultValue
      type(coordList), intent(in)         :: coords
      real(defReal), intent(in), optional :: mult
      real(defReal)                       :: val
    end function at

    !!
    !!
    !!
    elemental function getMaximumValue(self, defaultValue, materialIdx, mult) result(maximumValue)
      import                                  :: defReal, scalarField, shortInt
      class(scalarField), intent(in)          :: self
      real(defReal), intent(in)               :: defaultValue
      integer(shortInt), intent(in), optional :: materialIdx
      real(defReal), intent(in), optional     :: mult
      real(defReal)                           :: maximumValue
    end function getMaximumValue

    !!
    !!
    !!
    elemental function getMinimumValue(self, defaultValue, materialIdx, mult) result(minimumValue)
      import                                  :: defReal, scalarField, shortInt
      class(scalarField), intent(in)          :: self
      real(defReal), intent(in)               :: defaultValue
      integer(shortInt), intent(in), optional :: materialIdx
      real(defReal), intent(in), optional     :: mult
      real(defReal)                           :: minimumValue
    end function getMinimumValue

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
  !! Cast field pointer to scalarField pointer
  !!
  !! Args:
  !!   source [in] -> source pointer of class field
  !!
  !! Result:
  !!   Null is source is not of scalarField
  !!   Pointer to source if source is scalarField class
  !!
  function castScalarFieldPtr(source, fatal) result(ptr)
    class(field), intent(in)               :: source
    logical(defBool), intent(in), optional :: fatal
    class(scalarField), pointer            :: ptr
    logical(defBool)                       :: throwError
    character(*), parameter                :: HERE = 'castScalarFieldPtr (scalarField_inter.f90)'

    select type (temp => source)
      class is (scalarField)
        ptr => temp

      class default
        ptr => null()

    end select

    ! By default throw error if pointer is unassociated.
    throwError = .true.
    if (present(fatal)) throwError = fatal
    if (throwError .and. .not. associated(ptr)) call fatalError(HERE, "Field is not of class 'scalarField.")

  end function castScalarFieldPtr

  !!
  !!
  !!
  function getMaximumScalarFieldValue(fieldName, defaultValue, materialIdx, mult) result(maximumValue)
    character(nameLen), intent(in)          :: fieldName
    real(defReal), intent(in)               :: defaultValue
    integer(shortInt), intent(in), optional :: materialIdx
    real(defReal), intent(in), optional     :: mult
    class(field), pointer                   :: fieldPtr
    class(scalarField), pointer             :: scalarFieldPtr
    real(defReal)                           :: maximumValue

    fieldPtr => fieldPtrByName(fieldName)
    if (.not. associated(fieldPtr)) then
      maximumValue = defaultValue
      return

    end if

    scalarFieldPtr => castScalarFieldPtr(fieldPtr)
    maximumValue = scalarFieldPtr % getMaximumValue(defaultValue, materialIdx, mult)

  end function getMaximumScalarFieldValue

  !!
  !!
  !!
  function getMinimumScalarFieldValue(fieldName, defaultValue, materialIdx, mult) result(minimumValue)
    character(nameLen), intent(in)          :: fieldName
    real(defReal), intent(in)               :: defaultValue
    integer(shortInt), intent(in), optional :: materialIdx
    real(defReal), intent(in), optional     :: mult
    class(field), pointer                   :: fieldPtr
    class(scalarField), pointer             :: scalarFieldPtr
    real(defReal)                           :: minimumValue

    fieldPtr => fieldPtrByName(fieldName)
    if (.not. associated(fieldPtr)) then
      minimumValue = defaultValue
      return

    end if

    scalarFieldPtr => castScalarFieldPtr(fieldPtr)
    minimumValue = scalarFieldPtr % getMinimumValue(defaultValue, materialIdx, mult)

  end function getMinimumScalarFieldValue

  !!
  !!
  !!
  function getScalarFieldValue(fieldName, defaultValue, coords, mult) result(value)
    character(nameLen), intent(in)      :: fieldName
    real(defReal), intent(in)           :: defaultValue
    type(coordList), intent(in)         :: coords
    real(defReal), intent(in), optional :: mult
    class(field), pointer               :: fieldPtr
    class(scalarField), pointer         :: scalarFieldPtr
    real(defReal)                       :: value

    fieldPtr => fieldPtrByName(fieldName)
    if (.not. associated(fieldPtr)) then
      value = defaultValue
      return

    end if
    
    scalarFieldPtr => castScalarFieldPtr(fieldPtr)
    value = scalarFieldPtr % at(defaultValue, coords, mult)

  end function getScalarFieldValue

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(scalarField), intent(inout) :: self

    ! Local.

  end subroutine kill

end module scalarField_inter