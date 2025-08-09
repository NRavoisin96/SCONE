module OpenFOAMScalarField_class

  use dictionary_class,                               only : dictionary
  use genericProcedures,                              only : fatalError, openToRead
  use numPrecision
  use unstructuredPiecewiseConstantScalarField_inter, only : init_super => init, unstructuredPiecewiseConstantScalarField

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(unstructuredPiecewiseConstantScalarField) :: OpenFOAMScalarField
  contains
    procedure :: init
  end type OpenFOAMScalarField

contains
  !!
  !!
  !!
  subroutine init(self, dict)
    class(OpenFOAMScalarField), intent(inout) :: self
    class(dictionary), intent(in)             :: dict
    character(10)                             :: fieldType
    character(pathLen)                        :: path
    character(256)                            :: buffer
    character(:), allocatable                 :: dataString, trimmedPath
    integer(shortInt)                         :: idx, nValues
    logical(defBool)                          :: exists
    real(defReal)                             :: uniformValue
    real(defReal), dimension(:), allocatable  :: values
    integer(shortInt), parameter              :: unit = 10
    character(*), parameter                   :: here = 'init (OpenFOAMScalarField_class.f90)'

    ! Initialise superclass.
    call init_super(self, dict)

    ! Open the source file specified in the dictionary.
    call dict % get(path, 'path')
    trimmedPath = trim(path)
    inquire(file = trimmedPath, exist = exists)
    if (.not. exists) call fatalError(here, trimmedPath//' is not a valid file.')
    call openToRead(unit, trimmedPath)

    ! Parse the data contained in the file.
    read(unit, '(a)') buffer
    do while (index(buffer, 'internalField') == 0)
      read(unit, '(a)') buffer

    end do

    ! We are on the correct line. Read the next word to determine if field is uniform or nonuniform.
    read(buffer(index(buffer, 'internalField') + 13:), *) fieldType
    
    ! Handle logic based on whether the field is uniform or not.
    select case(trim(fieldType))
      case('uniform')
        read(buffer(index(buffer, 'uniform') + 7:), *) uniformValue
        allocate(values(self % getValuesNumber()))
        values = uniformValue

      case('nonuniform')
        idx = index(buffer, '(')
        if (0 < idx) then
          ! Single-line format.
          read(buffer(index(buffer, 'List<scalar>') + 12:idx - 1), *) nValues
          dataString = trim(buffer(idx + 1:index(buffer, ')') - 1))

        else
          ! Multi-line format. Move onto the next line to get the number of values.
          read(unit, '(a)') buffer
          read(buffer, *) nValues

          ! Move onto the next line, which should contain the opening parenthesis.
          read(unit, '(a)') buffer

          ! Initialise dataString and parse.
          dataString = ''
          do
            read(unit, '(a)') buffer
            if (0 < index(buffer, ')')) exit
            dataString = dataString//' '//trim(buffer)

          end do

        end if

        ! Check that the number of values matches the number of values in the field.
        if (0 < len_trim(dataString)) then
          allocate(values(nValues))
          read(dataString, *) values

        else
          call fatalError(here, 'Unable to parse internal field data.')

        end if

      case default
        call fatalError(here, 'Unknown internal scalar field type: '//trim(fieldType)//'.')

    end select

    ! Close the file and set the values in the field.
    close(unit)
    call self % setValues(values)

  end subroutine init

end module OpenFOAMScalarField_class