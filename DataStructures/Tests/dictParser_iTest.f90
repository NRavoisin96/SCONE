module dictParser_iTest
  
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : fileToDict
  use funit
  use numPrecision
  
  implicit none

contains

  !!
  !! Test Reading a Dictionary from a File
  !!
@Test
  subroutine testFromFile()
    type(dictionary)                             :: dict
    integer(shortInt)                            :: tempInt
    integer(longInt)                             :: tempLongInt
    real(defReal)                                :: tempReal
    class(dictionary), pointer                   :: dictPtr
    integer(shortInt), dimension(:), allocatable :: tempIntArray
    integer(longInt), dimension(:), allocatable  :: tempLongIntArray
    real(defReal), dimension(:), allocatable     :: tempRealArray

    call fileToDict(dict, './IntegrationTestFiles/testDictionary')

    ! Verify integer values
    call dict % get(tempInt, 'myInt')
    call dict % get(tempIntArray, 'intArray')

    @assertEqual(7, tempInt)
    @assertEqual([1, 2, 4, 5], tempIntArray)

    ! Verify longInt values.
    call dict % get(tempLongInt, 'myLongInt')
    call dict % get(tempLongIntArray, 'longIntArray')
    
    @assertEqual(777777777777_longInt, tempLongInt)
    @assertEqual([51234567890_longInt, -678998846213_longInt, 44498876563976_longInt], tempLongIntArray)

    ! Verify real values
    call dict % get(tempReal, 'myReal')
    call dict % get(tempRealArray, 'realArray')

    @assertEqual(1.3_defReal, tempReal)
    @assertEqual([1.0_defReal, 2.2_defReal, 3.5_defReal], tempRealArray)

    ! Verify nested dictionary
    dictPtr => dict % getDictPtr('subDict')
    call dictPtr % get(tempInt, 'myInt')
    call dictPtr % get(tempLongInt, 'myLongInt')
    call dictPtr % get(tempReal, 'myReal')

    @assertEqual(3, tempInt)
    @assertEqual(4_longInt, tempLongInt)
    @assertEqual(3.2_defReal, tempReal)

  end subroutine testFromFile


end module dictParser_iTest
