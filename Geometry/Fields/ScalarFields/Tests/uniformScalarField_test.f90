module uniformScalarField_test

  use coordList_class,          only : coordList
  use dictionary_class,         only : dictionary
  use field_inter,              only : field
  use funit
  use numPrecision
  use scalarField_inter,        only : scalarField, castScalarFieldPtr
  use uniformScalarField_class, only : uniformScalarField, uniformScalarField_TptrCast

  implicit none

contains

  !!
  !! Test Uniform Scalar Field
  !!
@Test
  subroutine test_uniformScalarField()
    type(uniformScalarField), target  :: fieldT
    class(field), pointer             :: ref
    class(scalarField), pointer       :: ptr
    type(uniformScalarField), pointer :: ptr2
    type(dictionary)                  :: dict
    type(coordList)                   :: coords
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Test invalid pointers
    ref => null()

    ptr2 => uniformScalarField_TptrCast(ref)
    @assertFalse(associated(ptr2))

    ! Test valid pointers
    ref => fieldT

    ptr => castScalarFieldPtr(ref, .false.)
    ptr2 => uniformScalarField_TptrCast(ref)

    @assertTrue(associated(ptr, fieldT))
    @assertTrue(associated(ptr2, fieldT))

    ! Initialise field
    call dict % init(2)
    call dict % store('type', 'uniformVectorField')
    call dict % store('value', 9.6_defReal)

    call fieldT % init(dict)

    ! Check value
    @assertEqual(9.6_defReal, fieldT % at(ZERO, coords), TOL)

    ! Kill
    call fieldT % kill()

  end subroutine test_uniformScalarField

end module uniformScalarField_test
