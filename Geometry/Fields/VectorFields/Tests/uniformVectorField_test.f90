module uniformVectorField_test

  use dictionary_class,          only : dictionary
  use field_inter,               only : field
  use funit
  use numPrecision
  use testTransportObject_class, only : testTransportObject
  use uniformVectorField_class,  only : uniformVectorField, uniformVectorField_TptrCast
  use vectorField_inter,         only : vectorField, vectorField_CptrCast

  implicit none

  ! Variables.
  type(testTransportObject)        :: testObject
  type(uniformVectorField), target :: fieldT

contains
@Before
  !!
  !!
  !!
  subroutine setUp()
    type(dictionary) :: dict

    ! Initialise field.
    call dict % init(2)
    call dict % store('type', 'uniformVectorField')
    call dict % store('value', [9.6_defReal, -8.0_defReal, 9.7_defReal])
    call fieldT % init(dict)

  end subroutine setUp

@After
  !!
  !!
  !!
  subroutine tearDown()

    ! Clean up.
    call fieldT % kill()

  end subroutine tearDown

  !!
  !! Test Uniform Scalar Field
  !!
@Test
  subroutine test_uniformVectorField()
    class(field), pointer             :: ref
    class(vectorField), pointer       :: vectorFieldPtr
    type(dictionary)                  :: dict
    type(uniformVectorField), pointer :: uniformVectorFieldPtr
    real(defReal), parameter          :: TOL = 1.0E-7_defReal

    ! Test invalid pointers
    ref => null()
    vectorFieldPtr => vectorField_CptrCast(ref)
    uniformVectorFieldPtr => uniformVectorField_TptrCast(ref)
    @assertFalse(associated(vectorFieldPtr))
    @assertFalse(associated(uniformVectorFieldPtr))

    ! Test valid pointers
    ref => fieldT
    vectorFieldPtr => vectorField_CptrCast(ref)
    uniformVectorFieldPtr => uniformVectorField_TptrCast(ref)
    @assertTrue(associated(vectorFieldPtr, fieldT))
    @assertTrue(associated(uniformVectorFieldPtr, fieldT))

    ! Check value
    @assertEqual([9.6_defReal, -8.0_defReal, 9.7_defReal], fieldT % at(testObject), TOL)

  end subroutine test_uniformVectorField

end module uniformVectorField_test