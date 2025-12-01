module uniformScalarField_class

  use coordList_class,   only : coordList
  use dictionary_class,  only : dictionary
  use errors_mod,        only : fatalError
  use field_inter,       only : field
  use numPrecision
  use scalarField_inter, only : kill_super => kill, scalarField

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: uniformScalarField_TptrCast

  !!
  !! Uniform Scalar Field
  !!
  !! Always returns the same value
  !!
  !! Sample Dictionary Input:
  !!   field { type uniformScalarField; value 3.0;}
  !!
  !! Public Members:
  !!   val -> Value of the field
  !!
  !! Interface:
  !!   scalarField interface
  !!
  type, public, extends(scalarField) :: uniformScalarField
    real(defReal) :: val = ZERO
  contains
    ! Superclass interface
    procedure :: at
    procedure :: getMaximumValue
    procedure :: getMinimumValue
    procedure :: init
    procedure :: kill
    procedure :: setValues
  end type uniformScalarField

contains
  !!
  !! Get value of the scalar field at the co-ordinate point
  !!
  !! See scalarField_inter for details
  !!
  function at(self, defaultValue, coords, mult) result(val)
    class(uniformScalarField), intent(in) :: self
    real(defReal), intent(in)             :: defaultValue
    type(coordList), intent(in)           :: coords
    real(defReal), intent(in), optional   :: mult
    real(defReal)                         :: val

    val = self % val
    if (present(mult)) val = val * mult

  end function at

  !!
  !! Initialise from dictionary
  !!
  !! See field_inter for details
  !!
  subroutine init(self, dict)
    class(uniformScalarField), intent(inout) :: self
    class(dictionary), intent(in)            :: dict

    ! Load value
    call dict % get(self % val, 'value')

  end subroutine init

  !!
  !!
  !!
  elemental function getMaximumValue(self, defaultValue, materialIdx, mult) result(maximumValue)
    class(uniformScalarField), intent(in)   :: self
    real(defReal), intent(in)               :: defaultValue
    integer(shortInt), intent(in), optional :: materialIdx
    real(defReal), intent(in), optional     :: mult
    real(defReal)                           :: maximumValue

    maximumValue = self % val
    if (present(mult)) maximumValue = maximumValue * mult

  end function getMaximumValue

  !!
  !!
  !!
  elemental function getMinimumValue(self, defaultValue, materialIdx, mult) result(minimumValue)
    class(uniformScalarField), intent(in)   :: self
    real(defReal), intent(in)               :: defaultValue
    integer(shortInt), intent(in), optional :: materialIdx
    real(defReal), intent(in), optional     :: mult
    real(defReal)                           :: minimumValue

    minimumValue = self % val
    if (present(mult)) minimumValue = minimumValue * mult

  end function getMinimumValue

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(uniformScalarField), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % val = ZERO

  end subroutine kill

  !!
  !!
  !!
  subroutine setValues(self, values)
    class(uniformScalarField), intent(inout) :: self
    real(defReal), dimension(:), intent(in)  :: values
    character(*), parameter                  :: here = 'setValues (uniformScalarField_class.f90)'

    if (1 < size(values)) call fatalError(here, 'Attempting to set more than one value for a uniform scalar field.')
    self % val = values(1)

  end subroutine setValues

  !!
  !! Cast field pointer to uniformScalarField pointer
  !!
  !! Args:
  !!   source [in] -> source pointer of class field
  !!
  !! Result:
  !!   Null is source is not of uniformScalarField
  !!   Pointer to source if source is uniformScalarField type
  !!
  pure function uniformScalarField_TptrCast(source) result(ptr)
    class(field), pointer, intent(in) :: source
    type(uniformScalarField), pointer :: ptr

    select type (source)
      type is (uniformScalarField)
        ptr => source

      class default
        ptr => null()
    end select

  end function uniformScalarField_TptrCast

end module uniformScalarField_class