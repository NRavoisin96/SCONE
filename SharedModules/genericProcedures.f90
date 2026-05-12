module genericProcedures

  use endfConstants
  use errors_mod,        only : fatalError
  use iso_fortran_env,   only : iostat_end ! Intrinsic fortran Modules
  use numPrecision
  use universalVariables

  implicit none

  interface areWithinTolerance
    module procedure areWithinTolerance_defReal
  end interface areWithinTolerance

  interface binarySearch
    module procedure binaryFloorIdxClosed_defReal
  end interface binarySearch

  interface ENDFInterpolate
    module procedure ENDFInterpolate_defReal_defReal
  end interface ENDFInterpolate

  interface getOrDefault
    module procedure getOrDefault_defBool
    module procedure getOrDefault_defReal
    module procedure getOrDefault_longInt
    module procedure getOrDefault_shortInt
  end interface getOrDefault

  interface hasDuplicates
    module procedure hasDuplicates_defReal
    module procedure hasDuplicates_longInt
    module procedure hasDuplicates_shortInt
  end interface hasDuplicates

  interface isIn
    module procedure isIn_defReal
    module procedure isIn_longInt
    module procedure isIn_shortInt
  end interface isIn

  interface isSortedAscending
    module procedure isSortedAscending_defReal
    module procedure isSortedAscending_longInt
    module procedure isSortedAscending_shortInt
  end interface isSortedAscending

  interface isSortedDescending
    module procedure isSortedDescending_defReal
    module procedure isSortedDescending_longInt
    module procedure isSortedDescending_shortInt
  end interface isSortedDescending

  interface linearLinearInterpolate
    module procedure linearLinearInterpolate_defReal_defReal
  end interface linearLinearInterpolate

  interface linearFind
    module procedure linearFind_char
    module procedure linearFind_defReal
    module procedure linearFind_longInt
    module procedure linearFind_shortInt
  end interface linearFind

  interface linearSearchCeil
    module procedure linearSearchCeilingIdxOpen_shortInt
  end interface linearSearchCeil

  interface linearSearchFloor
    module procedure linearSearchFloorIdxClosed_defReal
    module procedure linearSearchFloorIdxClosed_shortInt
  end interface linearSearchFloor

  interface LomutoPartition
    module procedure LomutoPartition_char
    module procedure LomutoPartition_defReal
    module procedure LomutoPartition_defReal_defReal
    module procedure LomutoPartition_longInt
    module procedure LomutoPartition_shortInt
  end interface LomutoPartition

  interface mergeSort
    module procedure mergeSort_defReal
    module procedure mergeSort_defReal_defReal
    module procedure mergeSort_shortInt
  end interface mergeSort

  interface numToChar
    module procedure numToChar_defReal
    module procedure numToChar_defRealArray
    module procedure numToChar_longInt
    module procedure numToChar_longIntArray
    module procedure numToChar_shortInt
    module procedure numToChar_shortIntArray
  end interface numToChar

  interface quickSort
    module procedure quickSort_char
    module procedure quickSort_defReal
    module procedure quickSort_defReal_defReal
    module procedure quickSort_longInt
    module procedure quickSort_shortInt
  end interface quickSort

  interface readArray
    module procedure readArray_defReal
    module procedure readArray_shortInt
  end interface readArray

  interface removeDuplicates
    module procedure removeDuplicates_char
    module procedure removeDuplicates_defReal
    module procedure removeDuplicates_shortInt
  end interface removeDuplicates

  interface swap
    module procedure swap_char_nameLen
    module procedure swap_defBool
    module procedure swap_defReal
    module procedure swap_defBool_defBool
    module procedure swap_defReal_defReal
    module procedure swap_shortInt
  end interface swap

  interface swapInArray
    module procedure swapInArray_char
    module procedure swapInArray_defReal
    module procedure swapInArray_longInt
    module procedure swapInArray_shortInt
  end interface swapInArray

contains
  !!
  !! Returns .true. if two floating point numbers are within tolerance of one another. 
  !!
  !! Due to floating point artihmetic and rounding-off errors being slightly different 
  !! across different architectures (eg, Intel vs ARM), it is necessary to use some small 
  !! tolerances to assert equality between two floating point numbers.
  !!
  elemental function areWithinTolerance_defReal(a, b, tolerance, relativeTolerance) result(areThey)
    real(defReal), intent(in)           :: a, b
    real(defReal), intent(in), optional :: tolerance, relativeTolerance
    logical(defBool)                    :: areThey
    real(defReal)                       :: absDiff, relativeToleranceToUse, toleranceToUse

    ! Initialise areThey = .true. and check for perfect (to the bit) equality, since we can 
    ! return early in this case.
    areThey = .true.
    if(a == b) return

    ! Determine tolerance to use.
    if(present(tolerance)) then
      toleranceToUse = tolerance

    else
      toleranceToUse = floatTol

    end if

    ! Compute the absolute value of the difference between the two floating point numbers.
    ! Note that if a and b are both very large and of opposite signs this can cause overflow.
    absDiff = abs(a - b)

    ! Check if absDiff is less than toleranceToUse and return if yes.
    if(absDiff < toleranceToUse) return

    ! Determine relative tolerance to use.
    if(present(relativeTolerance)) then
      relativeToleranceToUse = relativeTolerance

    else
      relativeToleranceToUse = FP_REL_TOL

    end if

    ! Check if a and b are within some small relative tolerance of each other and return if
    ! yes. Note that if a and b are both very small numbers, then multiplying by a small
    ! tolerance can cause underflow. This is why we check absolute tolerance first.
    if(absDiff < max(abs(a), abs(b)) * relativeToleranceToUse) return

    ! If reached here, a and b are not within absolute or relative tolerance of each other.
    ! update equal = .false.
    areThey = .false.

  end function areWithinTolerance_defReal

  !!
  !! Binary search for the largest smaller-or-equal element in the array
  !!
  !! Finds a location of a value in a sorted array by a binary search. Returns index of the
  !! "floor" of the bin in which value lies
  !!
  !! Args:
  !!   array [in] -> Sorted array of reals. Must be in increasing order (a_i <= a_j for j > i)
  !!   value [in] -> Value, which location is to be found.
  !!
  !! Returns:
  !!   Index of the "floor" in the array for the value.
  !!   -> if value == a(N) returns N-1
  !!   -> if value == a(1) returns 1
  !!
  !! Errors:
  !!   idx == valueOutideArray if value < a(1) .or. value > a(N)
  !!   idx == tooManyIter if search fails to terminate
  !!
  pure function binaryFloorIdxClosed_defReal(array, value) result(idx)
    real(defReal), dimension(:), intent(in) :: array
    real(defReal), intent(in)               :: value
    integer(shortInt)                       :: bottom, top, idx

    ! Initialise bottom = 1 and top = size(array).
    bottom = 1
    top = size(array)

    if(top == 1 .or. value < array(bottom) .or. array(top) < value) then
      ! If array contains only one element, or if value is outside array bounds, set idx = valueOutideArray and return.
      idx = valueOutsideArray
      return

    end if

    do while(1 < top - bottom)
      ! Calculate mid point. Note: add bottom and compute top - bottom to avoid potential integer overflows.
      idx = bottom + (top - bottom) / 2

      ! Binary Step.
      if(array(idx) <= value) then
        bottom = idx

      else
        top = idx

      end if

    end do
    idx = bottom

  end function binaryFloorIdxClosed_defReal

  !!
  !! Compares strings for equality. Ignores leading blanks.
  !!
  elemental function charCmp(char1, char2) result(areEqual)
    character(*), intent(in) :: char1, char2
    logical(defBool)         :: areEqual

    areEqual = adjustl(char1) == adjustl(char2)

  end function charCmp

  !!
  !! Convert character to an Integer
  !!
  !! Looks at first 20 places in the character and if they contain a valid integer
  !! it performs conversion.
  !!
  !! E.G
  !! '2' -> 2 OK!
  !! '-003' -> -3 OK!
  !! '  2E+03' -> FAIL!
  !! ' 7 is Swell' -> FAIL!
  !! '7                   A' -> 7 OK!
  !!
  !! Args:
  !!   str [in]      -> character to convert
  !!   error [inout] -> Optional. Set to .TRUE. if conversion has failed
  !! Result:
  !!   An integer contained within first 20 fields of str.
  !! Error:
  !!   Fortran Intrinsic error upon failed conversion Unless error is present.
  !!   If error argument is present it is set to .TRUE. if conversion failed by any reason and
  !!   the value of i becomes undefined
  !!
  function charToInt(str, error) result(i)
    character(*), intent(in)                 :: str
    logical(defBool),intent(inout), optional :: error
    integer(shortInt)                        :: i
    integer(shortInt)                        :: state

    if(present(error)) then
      read(str, '(I20)', IOSTAT = state) i
      error = state /= 0

    else
      read(str, '(I20)') i

    end if

  end function charToInt

  !!
  !! Cross product for 3D vectors
  !!
  pure function crossProduct(a, b) result(c)
    real(defReal), dimension(3), intent(in) :: a,b
    real(defReal), dimension(3)             :: c

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)

  end function crossProduct

  !!
  !! Perform one of ENDF defined interpolation types
  !!
  function ENDFInterpolate_defReal_defReal(xMin, xMax, yMin, yMax, x, endfNum) result(y)
    real(defReal), intent(in)     :: xMin, xMax, yMin, yMax, x
    integer(shortInt), intent(in) :: endfNum
    real(defReal)                 :: y
    character(*), parameter       :: HERE = 'ENDFInterpolate_defReal_defReal (genericProcedures.f90)'

    select case (endfNum) ! Naming Convention for ENDF interpolation (inY-inX) i.e. log-lin => logarithmic in y; linear in x
      case (histogramInterpolation)
        y = yMin

      case (linLinInterpolation)
        y = linearLinearInterpolate(xMin, xMax, yMin, yMax, x)

      case (linLogInterpolation)
        y = linearLinearInterpolate(log(xMin), log(xMax), yMin, yMax, log(x))

      case (logLinInterpolation)
        y = exp(linearLinearInterpolate(xMin, xMax, log(yMin), log(yMax), x))

      case (loglogInterpolation)
        y = exp(linearLinearInterpolate(log(xMin), log(xMax), log(yMin), log(yMax), log(x)))

      case (chargedParticleInterpolation)
        ! Not implemented
        call fatalError(HERE, 'ENDF interpolation law for charged Particles is not implemented.')

      case default
        call fatalError(HERE, 'Unknown ENDF interpolation number.')

    end select

  end function ENDFInterpolate_defReal_defReal

  !! Function 'getOrDefault_defBool'
  !!
  !! Description:
  !!   Returns optional defBool if present. Else, returns a default defBool.
  !!
  !! Arguments:
  !!   defaultVaue [in]   -> Default defBool in case optionalValue is absent.
  !!   optionalValue [in] -> Optional defBool.
  !!
  !! Result:
  !!   value -> defBool that is equal to optionalValue if present, else is equal to defaultValue.
  !!
  elemental function getOrDefault_defBool(defaultValue, optionalValue) result(value)
    logical(defBool), intent(in)           :: defaultValue
    logical(defBool), intent(in), optional :: optionalValue
    logical(defBool)                       :: value

    if(present(optionalValue)) then
      value = optionalValue

    else
      value = defaultValue

    end if

  end function getOrDefault_defBool

  !! Function 'getOrDefault_defReal'
  !!
  !! Description:
  !!   Returns optional defReal if present. Else, returns a default defReal.
  !!
  !! Arguments:
  !!   defaultVaue [in]   -> Default defReal in case optionalValue is absent.
  !!   optionalValue [in] -> Optional defReal.
  !!
  !! Result:
  !!   value -> defReal that is equal to optionalValue if present, else is equal to defaultValue.
  !!
  elemental function getOrDefault_defReal(defaultValue, optionalValue) result(value)
    real(defReal), intent(in)           :: defaultValue
    real(defReal), intent(in), optional :: optionalValue
    real(defReal)                       :: value

    if(present(optionalValue)) then
      value = optionalValue

    else
      value = defaultValue

    end if

  end function getOrDefault_defReal

  !! Function 'getOrDefault_longInt'
  !!
  !! Description:
  !!   Returns optional longInt if present. Else, returns a default longInt.
  !!
  !! Arguments:
  !!   defaultVaue [in]   -> Default longInt in case optionalValue is absent.
  !!   optionalValue [in] -> Optional longInt.
  !!
  !! Result:
  !!   value -> longInt that is equal to optionalValue if present, else is equal to defaultValue.
  !!
  elemental function getOrDefault_longInt(defaultValue, optionalValue) result(value)
    integer(longInt), intent(in)           :: defaultValue
    integer(longInt), intent(in), optional :: optionalValue
    integer(longInt)                       :: value

    if(present(optionalValue)) then
      value = optionalValue

    else
      value = defaultValue

    end if

  end function getOrDefault_longInt

  !! Function 'getOrDefault_shortInt'
  !!
  !! Description:
  !!   Returns optional shortInt if present. Else, returns a default shortInt.
  !!
  !! Arguments:
  !!   defaultVaue [in]   -> Default shortInt in case optionalValue is absent.
  !!   optionalValue [in] -> Optional shortInt.
  !!
  !! Result:
  !!   value -> shortInt that is equal to optionalValue if present, else is equal to defaultValue.
  !!
  elemental function getOrDefault_shortInt(defaultValue, optionalValue) result(value)
    integer(shortInt), intent(in)           :: defaultValue
    integer(shortInt), intent(in), optional :: optionalValue
    integer(shortInt)                       :: value

    if(present(optionalValue)) then
      value = optionalValue

    else
      value = defaultValue

    end if

  end function getOrDefault_shortInt

  !!
  !! Returns true if array contains duplicates
  !!
  pure function hasDuplicates_defReal(array) result(doesIt)
    real(defReal), dimension(:), intent(in) :: array
    integer(shortInt)                       :: i
    logical(defBool)                        :: doesIt
    real(defReal), dimension(size(array))   :: temp

    ! Copy and sort array.
    temp = array
    call quickSort(temp)

    ! Loop through the array and return as soon as a duplicate is found.
    doesIt = .true.
    do i = 2, size(array)
      if(temp(i) == temp(i - 1)) return

    end do
    doesIt = .false.

  end function hasDuplicates_defReal

  !!
  !! Returns true if array contains duplicates
  !!
  pure function hasDuplicates_longInt(array) result(doesIt)
    integer(longInt), dimension(:), intent(in) :: array
    integer(shortInt)                          :: i
    integer(longInt), dimension(size(array))   :: temp
    logical(defBool)                           :: doesIt

    ! Copy and sort array.
    temp = array
    call quickSort(temp)

    ! Loop through the array and return as soon as a duplicate is found.
    doesIt = .true.
    do i = 2, size(array)
      if(temp(i) == temp(i - 1)) return

    end do
    doesIt = .false.

  end function hasDuplicates_longInt

  !!
  !! Returns true if array contains duplicates
  !!
  pure function hasDuplicates_shortInt(array) result(doesIt)
    integer(shortInt), dimension(:), intent(in) :: array
    integer(shortInt)                           :: i
    integer(shortInt), dimension(size(array))   :: temp
    logical(defBool)                            :: doesIt

    ! Copy and sort array.
    temp = array
    call quickSort(temp)

    ! Loop through the array and return as soon as a duplicate is found.
    doesIt = .true.
    do i = 2, size(array)
      if(temp(i) == temp(i - 1)) return

    end do
    doesIt = .false.

  end function hasDuplicates_shortInt

  !!
  !! Return true if key is in the array
  !!
  !! Args:
  !!   array [in] -> Array of data
  !!   key [in]   -> Required key
  !!
  !! Result:
  !!   True if key is in array. False otherwise.
  !!
  pure function isIn_defReal(array, key) result(isIt)
    real(defReal), dimension(:), intent(in) :: array
    real(defReal), intent(in)               :: key
    logical(defBool)                        :: isIt

    isIt = targetNotFound /= linearFind(array, key)

  end function isIn_defReal

  !!
  !! Return true if key is in the array
  !!
  !! Args:
  !!   array [in] -> Array of data
  !!   key [in]   -> Required key
  !!
  !! Result:
  !!   True if key is in array. False otherwise.
  !!
  pure function isIn_longInt(array, key) result(isIt)
    integer(longInt), dimension(:), intent(in) :: array
    integer(longInt), intent(in)               :: key
    logical(defBool)                           :: isIt

    isIt = targetNotFound /= linearFind(array, key)

  end function isIn_longInt

  !!
  !! Return true if key is in the array
  !!
  !! Args:
  !!   array [in] -> Array of data
  !!   key [in]   -> Required key
  !!
  !! Result:
  !!   True if key is in array. False otherwise.
  !!
  pure function isIn_shortInt(array, key) result(isIt)
    integer(shortInt), dimension(:), intent(in) :: array
    integer(shortInt), intent(in)               :: key
    logical(defBool)                            :: isIt

    isIt = targetNotFound /= linearFind(array, key)

  end function isIn_shortInt

  !!
  !! Checks if the provided float is an integer.
  !!
  elemental function isInteger(float) result (isIt)
    real(defReal), intent(in) :: float
    logical(defBool)          :: isIt

    ! Check if distance to nearest integer falls below small tolerance.
    isIt = abs(float - anint(float)) <= floatTol

  end function isInteger

  !!
  !! Function that check if the array is sorted in ascending order (a(i) >= a(i-1) for all i).
  !!
  pure function isSortedAscending_defReal(array) result (isIt)
    real(defReal), dimension(:), intent(in) :: array
    integer(shortInt)                       :: i
    logical(defBool)                        :: isIt

    isIt = .false.
    do i = 2, size(array)
      if(array(i) < array(i - 1)) return

    end do
    isIt = .true.

  end function isSortedAscending_defReal

  !!
  !! Function that check if the array is sorted in ascending order (a(i) >= a(i-1) for all i).
  !!
  pure function isSortedAscending_longInt(array) result (isIt)
    integer(longInt), dimension(:), intent(in) :: array
    integer(shortInt)                          :: i
    logical(defBool)                           :: isIt

    isIt = .false.
    do i = 2, size(array)
      if(array(i) < array(i - 1)) return

    end do
    isIt = .true.

  end function isSortedAscending_longInt

  !!
  !! Function that check if the array is sorted in ascending order (a(i) >= a(i-1) for all i).
  !!
  pure function isSortedAscending_shortInt(array) result (isIt)
    integer(shortInt), dimension(:), intent(in) :: array
    integer(shortInt)                           :: i
    logical(defBool)                            :: isIt

    isIt = .false.
    do i = 2, size(array)
      if(array(i) < array(i - 1)) return

    end do
    isIt = .true.

  end function isSortedAscending_shortInt

  !!
  !! Function that check if the array is sorted in descending order (a(i) <= a(i-1) for all i).
  !!
  function isSortedDescending_defReal(array) result (isIt)
    real(defReal),dimension(:),intent(in) :: array
    integer(shortInt)                     :: i
    logical(defBool)                      :: isIt

    isIt = .false.
    do i = 2, size(array)
      if(array(i - 1) < array(i)) return

    end do
    isIt = .true.

  end function isSortedDescending_defReal

  !!
  !! Function that check if the array is sorted in descending order (a(i) <= a(i-1) for all i).
  !!
  function isSortedDescending_longInt(array) result (isIt)
    integer(longInt),dimension(:),intent(in) :: array
    integer(shortInt)                        :: i
    logical(defBool)                         :: isIt

    isIt = .false.
    do i = 2, size(array)
      if(array(i - 1) < array(i)) return

    end do
    isIt = .true.

  end function isSortedDescending_longInt

  !!
  !! Function that check if the array is sorted in descending order (a(i) <= a(i-1) for all i).
  !!
  function isSortedDescending_shortInt(array) result (isIt)
    integer(shortInt),dimension(:),intent(in) :: array
    integer(shortInt)                         :: i
    logical(defBool)                          :: isIt

    isIt = .false.
    do i = 2, size(array)
      if(array(i - 1) < array(i)) return

    end do
    isIt = .true.

  end function isSortedDescending_shortInt

  !!
  !! Searches linearly for the occurance of target in charArray. Removes left blanks.
  !! Following Errors can occur:
  !! targetNotFound -> target is not present in the array
  !!
  pure function linearFind_char(array, target) result(idx)
    character(*), dimension(:), intent(in) :: array
    character(*), intent(in)               :: target
    integer(shortInt)                      :: idx

    do idx = 1, size(array)
      if(array(idx) == target) return

    end do
    idx = targetNotFound

  end function linearFind_char

  !!
  !! Searches linearly for the occurance of target in defRealArray.
  !!
  pure function linearFind_defReal(array, target, tolerance) result (idx)
    real(defReal), dimension(:), intent(in) :: array
    real(defReal), intent(in)               :: target
    real(defReal), intent(in), optional     :: tolerance
    integer(shortInt)                       :: idx
    real(defReal)                           :: toleranceToUse

    if(present(tolerance)) then
      toleranceToUse = tolerance

    else
      toleranceToUse = ZERO

    end if

    do idx= 1, size(array)
      if(abs(array(idx) - target) <= toleranceToUse) return

    end do
    idx = targetNotFound

  end function linearFind_defReal

  !!
  !! Searches linearly for the occurance of target in longIntArray.
  !!
  pure function linearFind_longInt(array, target) result (idx)
    integer(longInt), dimension(:), intent(in) :: array
    integer(longInt), intent(in)               :: target
    integer(shortInt)                          :: idx

    do idx = 1, size(array)
      if(array(idx) == target) return

    end do
    idx = targetNotFound

  end function linearFind_longInt

  !!
  !! Searches linearly for the occurance of target in shortIntArray.
  !!
  pure function linearFind_shortInt(array, target) result (idx)
    integer(shortInt), dimension(:), intent(in) :: array
    integer(shortInt), intent(in)               :: target
    integer(shortInt)                           :: idx

    do idx = 1, size(array)
      if(array(idx) == target) return

    end do
    idx = targetNotFound

  end function linearFind_shortInt

  !!
  !! Perform linear interpolation between defReals
  !!
  elemental function linearLinearInterpolate_defReal_defReal(xMin, xMax, yMin, yMax, x) result(y)
    real(defReal), intent(in) :: xMin, xMax, yMin, yMax, x
    real(defReal)             :: interpolationFactor, y

    interpolationFactor = (x - xMin) / (xMax - xMin)
    y = yMax * interpolationFactor + (ONE - interpolationFactor) * yMin

  end function linearLinearInterpolate_defReal_defReal

  !!
  !! Linear search for the smallest larger-or-equal element in an array of integers
  !!
  !! Finds a location of a value in a sorted array. Returns index of the
  !! "ceiling" of the bin in which the value lies
  !!
  !! Search is "open" below the array i.e. for value <= a(1) returns 1.
  !!
  !! Args:
  !!   array [in] -> Sorted array of shortInts. Must be in increasing order (a_i <= a_j for j > i).
  !!   value [in] -> Value, which location is to be found.
  !!
  !! Result:
  !!   Index of the "ceiling" in the array.
  !!   -> if value <= array(1) returns 1.
  !!   -> if array(size(array)) <= value return valueOutsideArray.
  !!
  pure function linearSearchCeilingIdxOpen_shortInt(array, value) result(idx)
    integer(shortInt), dimension(:), intent(in) :: array
    integer(shortInt), intent(in)               :: value
    integer(shortInt)                           :: idx

    do idx = 1, size(array)
      if(value <= array(idx)) return

    end do
    
    ! Value is larger than the upper bound of the array.
    idx = valueOutsideArray

  end function linearSearchCeilingIdxOpen_shortInt

  !!
  !! Linear search for the largest smaller-or-equal element in the array
  !!
  !! Finds a location of a value in a sorted array. Returns index of the
  !! "floor" of the bin in which value lies
  !!
  !! Args:
  !!   array [in] -> Sorted array of reals. Must be in increasing order (a_i <= a_j for j > i).
  !!   value [in] -> Value, which location is to be found.
  !!
  !! Result:
  !!   Index of the "floor" in the array for the value.
  !!   -> if value == array(1) returns 1.
  !!   -> if value == array(size(array)) returns size(array) - 1.
  !!   -> if value < array(1) .or. array(size(array)) < value returns valueOutsideArray.
  !!
  pure function linearSearchFloorIdxClosed_defReal(array, value) result (idx)
    real(defReal), dimension(:), intent(in) :: array
    real(defReal), intent(in)               :: value
    integer(shortInt)                       :: arraySize, idx

    arraySize = size(array)
    if(value < array(1) .or. array(arraySize) < value) then
      idx = valueOutsideArray
      return

    end if

    do idx = arraySize - 1, 1, -1
      if(array(idx) <= value) return

    end do

  end function linearSearchFloorIdxClosed_defReal

  !!
  !! Linear search for the largest smaller-or-equal element in the array of integers
  !!
  !! Finds a location of a value in a sorted array. Returns index of the
  !! "floor" of the bin in which value lies
  !!
  !! Args:
  !!   array [in] -> Sorted array of shortInts. Must be in increasing order (a_i <= a_j for j > i).
  !!   value [in] -> Value, which location is to be found.
  !!
  !! Result:
  !!   Index of the "floor" in the array for the value.
  !!   -> if value == array(1) returns 1.
  !!   -> if value == array(N) returns N - 1.
  !!
  !! Errors:
  !!   fatalError if value < array(1) .or. array(size(array)) <= value.
  !!
  function linearSearchFloorIdxClosed_shortInt(array, value) result(idx)
    integer(shortInt), dimension(:), intent(in) :: array
    integer(shortInt), intent(in)               :: value
    integer(shortInt)                           :: arraySize, idx
    character(*) ,parameter                     :: HERE = 'linearSearchFloorIdxClosed_shortInt (genericProcedures.f90)'

    ! Check if the value is outside array bounds.
    arraySize = size(array)
    if(value < array(1) .or. array(arraySize) <= value) call fatalError(HERE, 'Value is outside array bounds.')

    do idx = arraySize, 1, -1
      if(array(idx) <= value) return

    end do

  end function linearSearchFloorIdxClosed_shortInt

  !!
  !! Lomuto partitioning for character string.
  !!
  pure subroutine LomutoPartition_char(low, high, array, pivotIdx)
    integer(shortInt), intent(in)                   :: low, high
    character(len = *), dimension(:), intent(inout) :: array
    integer(shortInt), intent(out)                  :: pivotIdx
    character(len = size(array))                    :: pivot
    integer(shortInt)                               :: i, i_plus_one, j, middle

    ! Compute middle index then swap array(middle) with array(high).
    middle = low + (high - low) / 2
    call swapInArray_char(middle, high, array)

    ! Initialise pivot then sort.
    pivot = array(high)
    i = low - 1
    do j = low, high - 1
      if(array(j) <= pivot) then
        ! Increment i then swap array(i) with array(j).
        i = i + 1
        call swapInArray_char(i, j, array)

      end if

    end do

    ! Now set the pivot into its sorted position and return the partition index.
    i_plus_one = i + 1
    call swapInArray(i_plus_one, high, array)
    pivotIdx = i_plus_one

  end subroutine LomutoPartition_char

  !!
  !! Lomuto partitioning for defReal.
  !!
  pure subroutine LomutoPartition_defReal(low, high, array, pivotIdx)
    integer(shortInt), intent(in)              :: low, high
    real(defReal), dimension(:), intent(inout) :: array
    integer(shortInt), intent(out)             :: pivotIdx
    integer(shortInt)                          :: i, i_plus_one, j, middle
    real(defReal)                              :: pivot

    ! Compute middle index then swap array(middle) with array(high).
    middle = low + (high - low) / 2
    call swapInArray_defReal(middle, high, array)

    ! Initialise pivot then sort.
    pivot = array(high)
    i = low - 1
    do j = low, high - 1
      if(array(j) <= pivot) then
        ! Increment i then swap array(i) with array(j).
        i = i + 1
        call swapInArray_defReal(i, j, array)

      end if

    end do

    ! Now set the pivot into its sorted position and return the pivot index.
    i_plus_one = i + 1
    call swapInArray_defReal(i_plus_one, high, array)
    pivotIdx = i_plus_one

  end subroutine LomutoPartition_defReal

  !!
  !! Lomuto partitioning for defReal.
  !!
  pure subroutine LomutoPartition_defReal_defReal(low, high, array1, array2, pivotIdx)
    integer(shortInt), intent(in)              :: low, high
    real(defReal), dimension(:), intent(inout) :: array1, array2
    integer(shortInt), intent(out)             :: pivotIdx
    integer(shortInt)                          :: i, i_plus_one, j, middle
    real(defReal)                              :: pivot

    ! Compute middle index then swap elements in array1 and array2.
    middle = low + (high - low) / 2
    call swapInArray_defReal(middle, high, array1)
    call swapInArray_defReal(middle, high, array2)

    ! Initialise pivot then sort.
    pivot = array1(high)
    i = low - 1
    do j = low, high - 1
      if(array1(j) <= pivot) then
        ! Increment i then swap elements in array1 and array2.
        i = i + 1
        call swapInArray_defReal(i, j, array1)
        call swapInArray_defReal(i, j, array2)

      end if

    end do

    ! Now set the pivot into its sorted position and return the pivot index.
    i_plus_one = i + 1
    call swapInArray_defReal(i_plus_one, high, array1)
    call swapInArray_defReal(i_plus_one, high, array2)
    pivotIdx = i_plus_one

  end subroutine LomutoPartition_defReal_defReal

  !!
  !! Lomuto partitioning for shortInt.
  !!
  pure subroutine LomutoPartition_longInt(low, high, array, pivotIdx)
    integer(shortInt), intent(in)                 :: low, high
    integer(longInt), dimension(:), intent(inout) :: array
    integer(shortInt), intent(out)                :: pivotIdx
    integer(longInt)                              :: pivot
    integer(shortInt)                             :: i, i_plus_one, j, middle

    ! Compute middle index then swap array(middle) with array(high).
    middle = low + (high - low) / 2
    call swapInArray_longInt(middle, high, array)

    ! Initialise pivot then sort.
    pivot = array(high)
    i = low - 1
    do j = low, high - 1
      if(array(j) <= pivot) then
        ! Increment i then swap array(i) with array(j).
        i = i + 1
        call swapInArray_longInt(i, j, array)

      end if

    end do

    ! Now set the pivot into its sorted position and return the pivot index.
    i_plus_one = i + 1
    call swapInArray_longInt(i_plus_one, high, array)
    pivotIdx = i_plus_one

  end subroutine LomutoPartition_longInt

  !!
  !! Lomuto partitioning for shortInt.
  !!
  pure subroutine LomutoPartition_shortInt(low, high, array, pivotIdx)
    integer(shortInt), intent(in)                  :: low, high
    integer(shortInt), dimension(:), intent(inout) :: array
    integer(shortInt), intent(out)                 :: pivotIdx
    integer(shortInt)                              :: i, i_plus_one, j, middle, pivot

    ! Compute middle index then swap array(middle) with array(high).
    middle = low + (high - low) / 2
    call swapInArray_shortInt(middle, high, array)

    ! Initialise pivot then sort.
    pivot = array(high)
    i = low - 1
    do j = low, high - 1
      if(array(j) <= pivot) then
        ! Increment i then swap array(i) with array(j).
        i = i + 1
        call swapInArray_shortInt(i, j, array)

      end if

    end do

    ! Now set the pivot into its sorted position and return the pivot index.
    i_plus_one = i + 1
    call swapInArray_shortInt(i_plus_one, high, array)
    pivotIdx = i_plus_one

  end subroutine LomutoPartition_shortInt

  !! Subroutine 'mergeSort_defReal'
  !!
  !! Description:
  !!   Merge sort algorithm for defReal arrays. Has a time complexity of O(N * log(N)) and a space complexity of O(N), where
  !!   N = size(array). Merge sort is a stable sorting algorithm, meaning that the relative order of identical elements in
  !!   the input is preserved.
  !!
  !! Arguments:
  !!   array [inout] -> defReal array to be sorted.
  !!
  pure subroutine mergeSort_defReal(array)
    real(defReal), dimension(:), intent(inout) :: array
    integer(shortInt)                          :: arraySize, high, i, j, k, l, low, middle, two_widths, width
    real(defReal), dimension(:), allocatable   :: temp

    ! Compute arraySize and immediately return if arraySize < 2.
    arraySize = size(array)
    if(arraySize < 2) return

    ! Allocate temp, initialise width = 1 then sort array.
    allocate(temp(arraySize))
    width = 1
    do while(width < arraySize)
      ! Initialise i = 1 then pre-compute two_widths.
      i = 1
      two_widths = 2 * width
      do while(i < arraySize)
        ! Initialise low, middle, and high.
        low = i
        middle = min(i + width - 1, arraySize)
        high = min(i + two_widths - 1, arraySize)

        ! Initialise j, k, and l then build temp array.
        j = low
        k = middle + 1
        l = low
        do while(j <= middle .and. k <= high)
          ! Compare values and pull from left or right half accordingly.
          if(array(j) <= array(k)) then
            temp(l) = array(j)
            j = j + 1

          else
            temp(l) = array(k)
            k = k + 1

          end if
          ! Update current index in temp array.
          l = l + 1

        end do

        ! Add any leftovers in the left half.
        do while(j <= middle)
          temp(l) = array(j)
          j = j + 1
          l = l + 1

        end do

        ! Add any leftovers in the right half.
        do while(k <= high)
          temp(l) = array(k)
          k = k + 1
          l = l + 1

        end do
        ! Update i.
        i = i + two_widths

      end do
      ! Update array and width.
      array = temp
      width = two_widths

    end do

  end subroutine mergeSort_defReal

  !! Subroutine 'mergeSort_defReal_defReal'
  !!
  !! Description:
  !!   Merge sort algorithm for pairs of real arrays. array2 is sorted according to array1. Has a time complexity of 
  !!   O(N * log(N)) and a space complexity of O(N), where N = size(array1). Merge sort is a stable sorting algorithm, meaning 
  !!   that the relative order of identical elements in the input is preserved.
  !!
  !! Arguments:
  !!   array1 [inout] -> defReal array to be sorted.
  !!   array2 [inout] -> defReal array sorted according to array1.
  !!
  !! Errors:
  !!   - fatalError if size(array1) /= size(array2).
  !!
  subroutine mergeSort_defReal_defReal(array1, array2)
    real(defReal), dimension(:), intent(inout) :: array1, array2
    integer(shortInt)                          :: arraySize, high, i, j, k, l, low, middle, two_widths, width
    real(defReal), dimension(:), allocatable   :: temp1, temp2
    character(*), parameter                    :: HERE = 'mergeSort_defReal_defReal (genericProcedures.f90)'

    ! Compute arraySize and immediately return if arraySize < 2.
    arraySize = size(array1)
    if(arraySize /= size(array2)) call fatalError(HERE, 'Arrays have different sizes.')
    if(arraySize < 2) return

    ! Allocate temp, initialise width = 1 then sort array.
    allocate(temp1(arraySize), temp2(arraySize))
    width = 1
    do while(width < arraySize)
      ! Initialise i = 1 then pre-compute two_widths.
      i = 1
      two_widths = 2 * width
      do while(i < arraySize)
        ! Initialise low, middle, and high.
        low = i
        middle = min(i + width - 1, arraySize)
        high = min(i + two_widths - 1, arraySize)

        ! Initialise j, k, and l then build temp array.
        j = low
        k = middle + 1
        l = low
        do while(j <= middle .and. k <= high)
          ! Compare values and pull from left or right half accordingly.
          if(array1(j) <= array1(k)) then
            temp1(l) = array1(j)
            temp2(l) = array2(j)
            j = j + 1

          else
            temp1(l) = array1(k)
            temp2(l) = array2(k)
            k = k + 1

          end if
          ! Update current index in temp array.
          l = l + 1

        end do

        ! Add any leftovers in the left half.
        do while(j <= middle)
          temp1(l) = array1(j)
          temp2(l) = array2(j)
          j = j + 1
          l = l + 1

        end do

        ! Add any leftovers in the right half.
        do while(k <= high)
          temp1(l) = array1(k)
          temp2(l) = array2(k)
          k = k + 1
          l = l + 1

        end do
        ! Update i.
        i = i + two_widths

      end do
      ! Update array and width.
      array1 = temp1
      array2 = temp2
      width = two_widths

    end do

  end subroutine mergeSort_defReal_defReal

  !! Subroutine 'mergeSort_shortInt'
  !!
  !! Description:
  !!   Merge sort algorithm for shortInt arrays. Has a time complexity of O(N * log(N)) and a space complexity of O(N), where
  !!   N = size(array). Merge sort is a stable sorting algorithm, meaning that the relative order of identical elements in
  !!   the input is preserved.
  !!
  !! Arguments:
  !!   array [inout] -> shortInt array to be sorted.
  !!
  pure subroutine mergeSort_shortInt(array)
    integer(shortInt), dimension(:), intent(inout) :: array
    integer(shortInt)                              :: arraySize, high, i, j, k, l, low, middle, two_widths, width
    integer(shortInt), dimension(:), allocatable   :: temp

    ! Compute arraySize and immediately return if arraySize < 2.
    arraySize = size(array)
    if(arraySize < 2) return

    ! Allocate temp, initialise width = 1 then sort array.
    allocate(temp(arraySize))
    width = 1
    do while(width < arraySize)
      ! Initialise i = 1 then pre-compute two_widths.
      i = 1
      two_widths = 2 * width
      do while(i < arraySize)
        ! Initialise low, middle, and high.
        low = i
        middle = min(i + width - 1, arraySize)
        high = min(i + two_widths - 1, arraySize)

        ! Initialise j, k, and l then build temp array.
        j = low
        k = middle + 1
        l = low
        do while(j <= middle .and. k <= high)
          ! Compare values and pull from left or right half accordingly.
          if(array(j) <= array(k)) then
            temp(l) = array(j)
            j = j + 1

          else
            temp(l) = array(k)
            k = k + 1

          end if
          ! Update current index in temp array.
          l = l + 1

        end do

        ! Add any leftovers in the left half.
        do while(j <= middle)
          temp(l) = array(j)
          j = j + 1
          l = l + 1

        end do

        ! Add any leftovers in the right half.
        do while(k <= high)
          temp(l) = array(k)
          k = k + 1
          l = l + 1

        end do
        ! Update i.
        i = i + two_widths

      end do
      ! Update array and width.
      array = temp
      width = two_widths

    end do

  end subroutine mergeSort_shortInt

  !!
  !! Converts defReal to character.
  !!
  pure function numToChar_defReal(x) result(c)
    real(defReal), intent(in) :: x
    character(:), allocatable :: c
    character(30)             :: tempChar ! Note: slightly larger than theoretical maximum (about 25) just to be extra safe.

    write(tempChar, '(g0)') x
    c = trim(adjustl(tempChar))

  end function numToChar_defReal

  !!
  !! Converts defReal array to character
  !!
  pure function numToChar_defRealArray(array) result(c)
    real(defReal), dimension(:),intent(in)  :: array
    character(:), allocatable               :: buffer, c
    character(30)                           :: tempChar ! Note: slightly larger than theoretical maximum (about 25) just to be extra safe.
    integer(shortInt)                       :: arraySize, i, length, position

    ! Compute array size and handle degenerate cases.
    arraySize = size(array)
    if(arraySize == 0) then
      c = ''
      return

    end if

    ! Allocate memory. Note: add 1 to len(tempChar) to account for blanks between consecutive elements.
    allocate(character((len(tempChar) + 1) * arraySize) :: buffer)
    
    ! Fill each element into the buffer.
    position = 1
    do i = 1, arraySize
      tempChar = numToChar_defReal(array(i))
      length = len_trim(tempChar)
      buffer(position:position + length - 1) = tempChar(1:length)
      position = position + length
      if(i < arraySize) then
        ! Add a space unless we are at the last element.
        buffer(position:position) = ' '
        position = position + 1

      end if

    end do
    c = buffer(1:position - 1)

  end function numToChar_defRealArray

  !!
  !! Converts longInt to character
  !!
  pure function numToChar_longInt(x) result(c)
    integer(longInt),intent(in) :: x
    character(:), allocatable   :: c
    character(20)               :: tempChar 

    write(tempChar, '(I0)') x
    c = trim(adjustl(tempChar))

  end function numToChar_longInt

  !!
  !! Converts longInt array to character
  !!
  pure function numToChar_longIntArray(array) result(c)
    integer(longInt), dimension(:),intent(in) :: array
    character(:), allocatable                 :: buffer, c
    character(20)                             :: tempChar
    integer(shortInt)                         :: arraySize, i, length, position

    ! Compute array size and handle degenerate cases.
    arraySize = size(array)
    if(arraySize == 0) then
      c = ''
      return

    end if

    ! Allocate memory. Note: add 1 to len(tempChar) to account for blanks between consecutive elements.
    allocate(character((len(tempChar) + 1) * arraySize) :: buffer)
    
    ! Fill each element into the buffer.
    position = 1
    do i = 1, arraySize
      tempChar = numToChar_longInt(array(i))
      length = len_trim(tempChar)
      buffer(position:position + length - 1) = tempChar(1:length)
      position = position + length
      if(i < arraySize) then
        ! Add a space unless we are at the last element.
        buffer(position:position) = ' '
        position = position + 1

      end if

    end do
    c = buffer(1:position - 1)

  end function numToChar_longIntArray

  !!
  !! Converts shortInt to character
  !!
  pure function numToChar_shortInt(x) result(c)
    integer(shortInt),intent(in) :: x
    character(:), allocatable    :: c
    character(11)                :: tempChar ! Note: 11 is the worst-case scenario for a shortInt.

    write(tempChar, '(I0)') x
    c = trim(adjustl(tempChar))

  end function numToChar_shortInt

  !!
  !! Converts shortInt array to character
  !!
  pure function numToChar_shortIntArray(array) result(c)
    integer(shortInt), dimension(:), intent(in) :: array
    character(11)                               :: tempChar ! Note: 11 is the worst-case scenario for a shortInt.
    character(:), allocatable                   :: buffer, c
    integer(shortInt)                           :: arraySize, i, length, position

    ! Compute array size and handle degenerate cases.
    arraySize = size(array)
    if(arraySize == 0) then
      c = ''
      return

    end if

    ! Allocate memory. Note: add 1 to len(tempChar) to account for blanks between consecutive elements.
    allocate(character((len(tempChar) + 1) * arraySize) :: buffer)
    
    ! Fill each element into the buffer.
    position = 1
    do i = 1, arraySize
      tempChar = numToChar_shortInt(array(i))
      length = len_trim(tempChar)
      buffer(position:position + length - 1) = tempChar(1:length)
      position = position + length
      if(i < arraySize) then
        ! Add a space unless we are at the last element.
        buffer(position:position) = ' '
        position = position + 1

      end if

    end do
    c = buffer(1:position - 1)

  end function numToChar_shortIntArray

  !!
  !! Open "file" for reading under with "unitNum" reference
  !!
  subroutine openToRead(unitNum, file)
    integer(shortInt), intent(in) :: unitNum
    character(*), intent(in)      :: file
    character(:), allocatable     :: errorMsg
    integer(shortInt)             :: errorStat
    character(*), parameter       :: HERE = 'openToRead (genericProcedures.f90)'

    open(unit = unitNum, file = file, status = "old", action = "read", iostat = errorStat, iomsg  = errorMsg)
    if(0 < errorStat) call fatalError(HERE, errorMsg)

  end subroutine openToRead

  !!
  !! Convert Particle Type to string
  !!
  !! Args:
  !!   type [in] -> particle type
  !!
  !! Result:
  !!   Allocatable String that describes what particle is this
  !!
  !! Errors:
  !!   For unknown type prints "Unknown <int>" where int is number in type
  !!
  pure function printParticleType(type) result(str)
    integer(shortInt), intent(in) :: type
    character(:), allocatable     :: str

    select case(type)
      case(P_NEUTRON_CE)
        str = "CE Neutron"

      case(P_NEUTRON_MG)
        str = "MG Neutron"

      case default
        str = "Unknown "//numToChar(type)

    end select

  end function printParticleType

  !!
  !! Quicksort for characters.
  !!
  pure subroutine quickSort_char(array)
    character(*), dimension(:), intent(inout) :: array
    integer(shortInt)                         :: arraySize, high, low, pivotIdx, pivotIdx_minus_one, pivotIdx_plus_one, stackTop
    integer(shortInt), dimension(128)         :: stack

    ! Handle degenerate cases.
    arraySize = len(array)
    if(arraySize < 2) return

    ! Initialise stack.
    stackTop = 1
    stack(stackTop) = 1
    stackTop = stackTop + 1
    stack(stackTop) = arraySize

    ! Loop until stack is empty.
    do while(0 < stackTop)
      ! Pop from the stack.
      high = stack(stackTop)
      stackTop = stackTop - 1
      low = stack(stackTop)
      stackTop = stackTop - 1

      ! Partition the array.
      call LomutoPartition_char(low, high, array, pivotIdx)

      ! If there are elements on the left side of the pivot, push the left sub-array bounds into the stack.
      pivotIdx_minus_one = pivotIdx - 1
      pivotIdx_plus_one = pivotIdx + 1

      if(pivotIdx_minus_one - low < high - pivotIdx_plus_one) then
        if(pivotIdx_plus_one < high) then
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_plus_one
          stackTop = stackTop + 1
          stack(stackTop) = high

        end if

        if(low < pivotIdx_minus_one) then
          stackTop = stackTop + 1
          stack(stackTop) = low
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_minus_one

        end if

      else
        if(low < pivotIdx_minus_one) then
          stackTop = stackTop + 1
          stack(stackTop) = low
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_minus_one

        end if

        if(pivotIdx_plus_one < high) then
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_plus_one
          stackTop = stackTop + 1
          stack(stackTop) = high

        end if

      end if

    end do

  end subroutine quickSort_char

  !!
  !! Quicksort for real array
  !!
  pure subroutine quickSort_defReal(array)
    real(defReal), dimension(:), intent(inout) :: array
    integer(shortInt)                          :: arraySize, high, low, pivotIdx, pivotIdx_minus_one, pivotIdx_plus_one, stackTop
    integer(shortInt), dimension(128)          :: stack

    ! Handle degenerate cases.
    arraySize = size(array)
    if(arraySize < 2) return

    ! Initialise stack.
    stack(1) = 1
    stack(2) = arraySize
    stackTop = 2

    ! Loop until stack is empty.
    do while(0 < stackTop)
      ! Pop from the stack.
      high = stack(stackTop)
      stackTop = stackTop - 1
      low = stack(stackTop)
      stackTop = stackTop - 1

      ! Partition the array.
      call LomutoPartition_defReal(low, high, array, pivotIdx)

      ! If there are elements on the left side of the pivot, push the left sub-array bounds into the stack.
      pivotIdx_minus_one = pivotIdx - 1
      pivotIdx_plus_one = pivotIdx + 1

      if(pivotIdx_minus_one - low < high - pivotIdx_plus_one) then
        if(pivotIdx_plus_one < high) then
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_plus_one
          stackTop = stackTop + 1
          stack(stackTop) = high

        end if

        if(low < pivotIdx_minus_one) then
          stackTop = stackTop + 1
          stack(stackTop) = low
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_minus_one

        end if

      else
        if(low < pivotIdx_minus_one) then
          stackTop = stackTop + 1
          stack(stackTop) = low
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_minus_one

        end if

        if(pivotIdx_plus_one < high) then
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_plus_one
          stackTop = stackTop + 1
          stack(stackTop) = high

        end if

      end if

    end do

  end subroutine quickSort_defReal

  !!
  !! Quicksort for array1 and array2 of reals by array1.
  !!
  subroutine quickSort_defReal_defReal(array1, array2)
    real(defReal), dimension(:), intent(inout) :: array1, array2
    integer(shortInt)                          :: array1Size, high, low, pivotIdx, pivotIdx_minus_one, pivotIdx_plus_one, &
                                                  stackTop
    integer(shortInt), dimension(128)          :: stack
    character(*), parameter                    :: HERE = 'quickSort_defReal_defReal (genericProcedures.f90)'

    ! Handle degenerate cases.
    array1Size = size(array1)
    if(array1Size /= size(array2)) call fatalError(HERE, 'Arrays have different sizes.')
    if(array1Size < 2) return

    ! Initialise stack.
    stack(1) = 1
    stack(2) = array1Size
    stackTop = 2

    ! Loop until stack is empty.
    do while(0 < stackTop)
      ! Pop from the stack.
      high = stack(stackTop)
      stackTop = stackTop - 1
      low = stack(stackTop)
      stackTop = stackTop - 1

      ! Partition the array.
      call LomutoPartition_defReal_defReal(low, high, array1, array2, pivotIdx)

      ! If there are elements on the left side of the pivot, push the left sub-array bounds into the stack.
      pivotIdx_minus_one = pivotIdx - 1
      pivotIdx_plus_one = pivotIdx + 1

      if(pivotIdx_minus_one - low < high - pivotIdx_plus_one) then
        if(pivotIdx_plus_one < high) then
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_plus_one
          stackTop = stackTop + 1
          stack(stackTop) = high

        end if

        if(low < pivotIdx_minus_one) then
          stackTop = stackTop + 1
          stack(stackTop) = low
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_minus_one

        end if

      else
        if(low < pivotIdx_minus_one) then
          stackTop = stackTop + 1
          stack(stackTop) = low
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_minus_one

        end if

        if(pivotIdx_plus_one < high) then
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_plus_one
          stackTop = stackTop + 1
          stack(stackTop) = high

        end if

      end if

    end do

  end subroutine quickSort_defReal_defReal

  !!
  !! Quicksort for longInt array.
  !!
  pure subroutine quickSort_longInt(array)
    integer(longInt), dimension(:), intent(inout) :: array
    integer(shortInt)                             :: arraySize, high, low, pivotIdx, pivotIdx_minus_one, pivotIdx_plus_one, &
                                                     stackTop
    integer(shortInt), dimension(128)             :: stack

    ! Handle degenerate cases.
    arraySize = size(array)
    if(arraySize < 2) return

    ! Initialise stack.
    stackTop = 1
    stack(stackTop) = 1
    stackTop = stackTop + 1
    stack(stackTop) = arraySize

    ! Loop until stack is empty.
    do while(0 < stackTop)
      ! Pop from the stack.
      high = stack(stackTop)
      stackTop = stackTop - 1
      low = stack(stackTop)
      stackTop = stackTop - 1

      ! Partition the array.
      call LomutoPartition_longInt(low, high, array, pivotIdx)

      ! If there are elements on the left side of the pivot, push the left sub-array bounds into the stack.
      pivotIdx_minus_one = pivotIdx - 1
      pivotIdx_plus_one = pivotIdx + 1

      if(pivotIdx_minus_one - low < high - pivotIdx_plus_one) then
        if(pivotIdx_plus_one < high) then
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_plus_one
          stackTop = stackTop + 1
          stack(stackTop) = high

        end if

        if(low < pivotIdx_minus_one) then
          stackTop = stackTop + 1
          stack(stackTop) = low
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_minus_one

        end if

      else
        if(low < pivotIdx_minus_one) then
          stackTop = stackTop + 1
          stack(stackTop) = low
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_minus_one

        end if

        if(pivotIdx_plus_one < high) then
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_plus_one
          stackTop = stackTop + 1
          stack(stackTop) = high

        end if

      end if

    end do

  end subroutine quickSort_longInt

  !!
  !! Quicksort for integer array
  !!
  pure subroutine quickSort_shortInt(array)
    integer(shortInt), dimension(:), intent(inout) :: array
    integer(shortInt)                              :: arraySize, high, low, pivotIdx, pivotIdx_minus_one, &
                                                      pivotIdx_plus_one, stackTop
    integer(shortInt), dimension(128)              :: stack

    ! Handle degenerate cases.
    arraySize = size(array)
    if(arraySize < 2) return

    ! Initialise stack.
    stackTop = 1
    stack(stackTop) = 1
    stackTop = stackTop + 1
    stack(stackTop) = arraySize

    ! Loop until stack is empty.
    do while(0 < stackTop)
      ! Pop from the stack.
      high = stack(stackTop)
      stackTop = stackTop - 1
      low = stack(stackTop)
      stackTop = stackTop - 1

      ! Partition the array.
      call LomutoPartition_shortInt(low, high, array, pivotIdx)

      ! If there are elements on the left side of the pivot, push the left sub-array bounds into the stack.
      pivotIdx_minus_one = pivotIdx - 1
      pivotIdx_plus_one = pivotIdx + 1

      if(pivotIdx_minus_one - low < high - pivotIdx_plus_one) then
        if(pivotIdx_plus_one < high) then
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_plus_one
          stackTop = stackTop + 1
          stack(stackTop) = high

        end if

        if(low < pivotIdx_minus_one) then
          stackTop = stackTop + 1
          stack(stackTop) = low
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_minus_one

        end if

      else
        if(low < pivotIdx_minus_one) then
          stackTop = stackTop + 1
          stack(stackTop) = low
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_minus_one

        end if

        if(pivotIdx_plus_one < high) then
          stackTop = stackTop + 1
          stack(stackTop) = pivotIdx_plus_one
          stackTop = stackTop + 1
          stack(stackTop) = high

        end if

      end if

    end do

  end subroutine quickSort_shortInt

  !!
  !! Read a line from the source file in ASCII or binary format
  !! EOF is a logical output that is set to true if end of file is reached and false otherwise
  !!
  subroutine readArray_defReal(unit, readBinary, output, EOF)
    integer(shortInt), intent(in)              :: unit
    logical(defBool), intent(in)               :: readBinary
    real(defReal), dimension(:), intent(inout) :: output
    logical(defBool), intent(out)              :: EOF
    integer(shortInt)                          :: errorCode
    character(*), parameter                    :: HERE = 'readArray_defReal (genericProcedures.f90)'

    ! Read directly into the output array.
    if(readBinary) then
      read(unit, iostat = errorCode) output

    else
      read(unit, *, iostat = errorCode) output

    end if
    
    select case(errorCode)
      case(0)
        EOF = .false.

      case(iostat_end)
        EOF = .true.

      case default
        call fatalError(HERE, 'Error reading source file.')

    end select

  end subroutine readArray_defReal

  !!
  !! Read a line from the source file in ASCII or binary format
  !! EOF is a logical output that is set to true if end of file is reached and false otherwise
  !!
  subroutine readArray_shortInt(unit, readBinary, output, EOF)
    integer(shortInt), intent(in)                  :: unit
    logical(defBool), intent(in)                   :: readBinary
    integer(shortInt), dimension(:), intent(inout) :: output
    logical(defBool), intent(out)                  :: EOF
    integer(shortInt)                              :: errorCode
    character(*), parameter                        :: HERE = 'readArray_shortInt (genericProcedures.f90)'

    ! Read directly into the output array.
    if(readBinary) then
      read(unit, iostat = errorCode) output 

    else
      read(unit, *, iostat = errorCode) output

    end if

    select case(errorCode)
      case(0)
        EOF = .false.

      case(iostat_end)
        EOF = .true.

      case default
        call fatalError(HERE, 'Error reading source file.')

    end select

  end subroutine readArray_shortInt

  !!
  !! Removes duplicates from a nameLen character array.
  !! Does not preserve the order.
  !!
  pure function removeDuplicates_char(array, isArraySorted) result(out)
    character(nameLen), dimension(:), intent(in)  :: array
    logical(defBool), intent(in), optional        :: isArraySorted
    character(nameLen), dimension(size(array))    :: temp
    character(nameLen), dimension(:), allocatable :: out
    integer(shortInt)                             :: arraySize, i, j
    logical(defBool)                              :: sorted

    ! Handle degenerate cases.
    arraySize = size(array)
    if(arraySize < 2) then
      allocate(out(arraySize))
      if(arraySize == 1) out(1) = array(1)
      return

    end if

    ! Initialise sorted.
    if(present(isArraySorted)) then
      sorted = isArraySorted

    else
      sorted = .false.

    end if

    ! Copy and sort the array if necessary.
    temp = array
    if(.not. sorted) call quickSort(temp)

    ! Find unique elements from sorted array.
    j = 1
    do i = 2, arraySize
      if(temp(i) /= temp(j)) then
        j = j + 1
        temp(j) = temp(i)

      end if

    end do

    ! Allocate output.
    out = temp(1:j)

  end function removeDuplicates_char

  !!
  !! Removes duplicates from a defReal array.
  !! Does not preserve the order
  !!
  pure function removeDuplicates_defReal(array, isArraySorted) result(out)
    real(defReal), dimension(:), intent(in)  :: array
    logical(defBool), intent(in), optional   :: isArraySorted
    integer(shortInt)                        :: arraySize, i, j
    logical(defBool)                         :: sorted
    real(defReal), dimension(size(array))    :: temp
    real(defReal), dimension(:), allocatable :: out

    ! Handle degenerate cases.
    arraySize = size(array)
    if(arraySize < 2) then
      allocate(out(arraySize))
      if(arraySize == 1) out(1) = array(1)
      return

    end if

    ! Initialise sorted.
    if(present(isArraySorted)) then
      sorted = isArraySorted

    else
      sorted = .false.

    end if

    ! Copy and sort the array if necessary.
    temp = array
    if(.not. sorted) call quickSort(temp)

    ! Find unique elements from sorted array.
    j = 1
    do i = 2, arraySize
      if(temp(i) /= temp(j)) then
        j = j + 1
        temp(j) = temp(i)

      end if

    end do

    ! Allocate output.
    out = temp(1:j)

  end function removeDuplicates_defReal

  !!
  !! Removes duplicates from a shortInt array.
  !! Does not preserve the order.
  !!
  pure function removeDuplicates_shortInt(array, isArraySorted) result(out)
    integer(shortInt), dimension(:), intent(in)  :: array
    logical(defBool), intent(in), optional       :: isArraySorted
    integer(shortInt)                            :: arraySize, i, j
    integer(shortInt), dimension(size(array))    :: temp
    integer(shortInt), dimension(:), allocatable :: out
    logical(defBool)                             :: sorted

    ! Handle degenerate cases.
    arraySize = size(array)
    if(arraySize < 2) then
      allocate(out(arraySize))
      if(arraySize == 1) out(1) = array(1)
      return

    end if

    ! Initialise sorted.
    if(present(isArraySorted)) then
      sorted = isArraySorted

    else
      sorted = .false.

    end if

    ! Copy and sort the array if necessary.
    temp = array
    if(.not. sorted) call quickSort(temp)

    ! Find unique elements from sorted array.
    j = 1
    do i = 2, arraySize
      if(temp(i) /= temp(j)) then
        j = j + 1
        temp(j) = temp(i)

      end if

    end do

    ! Allocate output.
    out = temp(1:j)

  end function removeDuplicates_shortInt

  !!
  !! Replaces all symbols "old" with "new" in string
  !!
  pure subroutine replaceChar(old, new, string)
    character(1), intent(in)    :: old, new
    character(*), intent(inout) :: string
    integer(shortInt)           :: i

    do i = 1, len(string)
      if(string(i:i) == old) string(i:i) = new

    end do

  end subroutine replaceChar

  !!
  !! Subroutine takes a normilised direction vector dir and rotates it by cosine of a polar angle
  !! mu and azimuthal angle phi (in radians).
  !! Procedure will produce incorrect results WITHOUT error message if dir is not normalised
  !!
  pure function rotateVector(dir, mu, phi) result(newDir)
    real(defReal), dimension(3), intent(in)    :: dir
    real(defReal), intent(in)                  :: mu, phi
    real(defReal), dimension(3)                :: newDir
    real(defReal)                              :: cosPhi, inverse_sineIncidentAngle, sineIncidentAngle, sinePhi, &
                                                  sinePolarScatteringAngle

    ! Precalculate cosine and sine of polar angle
    sinePhi = sin(phi)
    cosPhi = cos(phi)

    ! Perform standard rotation. Note that indexes are parameterised
    sinePolarScatteringAngle = sqrt(max(ZERO, ONE - mu * mu))
    sineIncidentAngle = sqrt(max(ZERO, ONE - dir(3) * dir(3)))

    if(.not. areWithinTolerance(ZERO, sineIncidentAngle)) then
      inverse_sineIncidentAngle = ONE / sineIncidentAngle
      newDir(1) = mu * dir(1) + &
                  sinePolarScatteringAngle * (dir(1) * dir(3) * cosPhi - dir(2) * sinePhi) * inverse_sineIncidentAngle
      newDir(2) = mu * dir(2) + &
                  sinePolarScatteringAngle * (dir(2) * dir(3) * cosPhi + dir(1) * sinePhi) * inverse_sineIncidentAngle
      newDir(3) = mu * dir(3) - sinePolarScatteringAngle * sineIncidentAngle * cosPhi

    else
      ! Direction is parallel to the z-axis (u = v = 0).
      newDir(1) = sinePolarScatteringAngle * cosPhi
      newDir(2) = sinePolarScatteringAngle * sinePhi
      newDir(3) = mu * sign(ONE, dir(3))

    end if

  end function rotateVector

  !!
  !! Generate Euler rotation matrix using ZXZ convention
  !!
  !! Args:
  !!   matrix [out] -> Space for the matrix dimension(3,3)
  !!   phi [in] -> Initial rotation over Z axis [deg]. In 0 to 360.
  !!   theta [in] -> 2nd rotation over tranfromed X' axis [deg]. In 0 to 180.
  !!   psi [in] -> Final rotation over transformed Z' axis [deg]. In 0 to 360.
  !!
  !! Errors:
  !!   fatalError if any angle is beyond its range
  !!
  pure subroutine rotationMatrix(matrix, phi, theta, psi)
    real(defReal), dimension(3,3), intent(out) :: matrix
    real(defReal), intent(in)                  :: phi, theta, psi
    real(defReal)                              :: conv, cos_phi, cos_psi, cos_theta, phi_rad, psi_rad, sin_phi, sin_psi, &
                                                  sin_theta, theta_rad

    ! Evaluate trigonometric functions.
    conv = TWO_PI / 360.0_defReal
    phi_rad = phi * conv
    psi_rad = psi * conv
    theta_rad = theta * conv

    cos_phi = cos(phi_rad)
    cos_psi = cos(psi_rad)
    cos_theta = cos(theta_rad)
    sin_phi = sin(phi_rad)
    sin_psi = sin(psi_rad)
    sin_theta = sin(theta_rad)

    ! Assign matrix elements.
    matrix(1, 1) = cos_psi * cos_phi - cos_theta * sin_phi * sin_psi
    matrix(2, 1) = -sin_psi * cos_phi - cos_theta * sin_phi * cos_psi
    matrix(3, 1) = sin_theta * sin_phi

    matrix(1, 2) = cos_psi * sin_phi + cos_theta * cos_phi * sin_psi
    matrix(2, 2) = -sin_psi * sin_phi + cos_theta * cos_phi * cos_psi
    matrix(3, 2) = -sin_theta * cos_phi

    matrix(1, 3) = sin_psi * sin_theta
    matrix(2, 3) = cos_psi * sin_theta
    matrix(3, 3) = cos_theta

  end subroutine rotationMatrix

  !!
  !! Swap character of length nameLen
  !!
  elemental subroutine swap_char_nameLen(c1, c2)
    character(nameLen), intent(inout) :: c1, c2
    character(nameLen)                :: temp

    temp = c1
    c1 = c2
    c2 = temp

  end subroutine swap_char_nameLen

  !!
  !! Swap two bools
  !!
  elemental subroutine swap_defBool(l1, l2)
    logical(defBool), intent(inout) :: l1, l2
    logical(defBool)                :: temp

    temp = l1
    l1 = l2
    l2 = temp

  end subroutine swap_defBool

  !!
  !! Swap two reals
  !!
  elemental subroutine swap_defReal(r1,r2)
    real(defReal), intent(inout) :: r1, r2
    real(defReal)                :: temp

    temp = r1
    r1 = r2
    r2 = temp

  end subroutine swap_defReal

  !!
  !! Swap two integers
  !!
  elemental subroutine swap_shortInt(i1, i2)
    integer(shortInt), intent(inout) :: i1, i2
    integer(shortInt)                :: temp

    temp = i1
    i1 = i2
    i2 = temp

  end subroutine swap_shortInt

  !!
  !! Swap two pair of bools
  !!
  elemental subroutine swap_defBool_defBool(r1_1, r1_2, r2_1, r2_2)
    logical(defBool), intent(inout) :: r1_1
    logical(defBool), intent(inout) :: r1_2
    logical(defBool), intent(inout) :: r2_1
    logical(defBool), intent(inout) :: r2_2
    logical(defBool)                :: temp1, temp2

    ! Load first pair into temps
    temp1 = r1_1
    temp2 = r1_2

    ! Assign values of 2nd pair to 1st pair
    r1_1 = r2_1
    r1_2 = r2_2

    ! Assign values of 1st pair to 2nd pair
    r2_1 = temp1
    r2_2 = temp2

  end subroutine swap_defBool_defBool

  !!
  !! Swap two pair of reals
  !!
  elemental subroutine swap_defReal_defReal(r1_1, r1_2, r2_1, r2_2)
    real(defReal), intent(inout) :: r1_1
    real(defReal), intent(inout) :: r1_2
    real(defReal), intent(inout) :: r2_1
    real(defReal), intent(inout) :: r2_2
    real(defReal)                :: temp1, temp2

    ! Load first pair into temps
    temp1 = r1_1
    temp2 = r1_2

    ! Assign values of 2nd pair to 1st pair
    r1_1 = r2_1
    r1_2 = r2_2

    ! Assign values of 1st pair to 2nd pair
    r2_1 = temp1
    r2_2 = temp2

  end subroutine swap_defReal_defReal

  !! Subroutine 'swapInArray_char'
  !!
  !! Description:
  !!   Swaps two elements within a character array.
  !!
  !! Arguments:
  !!   idx1 [in]     -> Index of the first element.
  !!   idx2 [in]     -> Index of the second element.
  !!   array [inout] -> character array.
  !!
  pure subroutine swapInArray_char(idx1, idx2, array)
    integer(shortInt), intent(in)             :: idx1, idx2
    character(*), dimension(:), intent(inout) :: array
    character(:), allocatable                 :: temp

    temp = array(idx1)
    array(idx1) = array(idx2)
    array(idx2) = temp

  end subroutine swapInArray_char

  !! Subroutine 'swapInArray_defReal'
  !!
  !! Description:
  !!   Swaps two elements within a defReal array.
  !!
  !! Arguments:
  !!   idx1 [in]     -> Index of the first element.
  !!   idx2 [in]     -> Index of the second element.
  !!   array [inout] -> defReal array.
  !!
  pure subroutine swapInArray_defReal(idx1, idx2, array)
    integer(shortInt), intent(in)              :: idx1, idx2
    real(defReal), dimension(:), intent(inout) :: array
    real(defReal)                              :: temp

    temp = array(idx1)
    array(idx1) = array(idx2)
    array(idx2) = temp

  end subroutine swapInArray_defReal

  !! Subroutine 'swapInArray_longInt'
  !!
  !! Description:
  !!   Swaps two elements within a longInt array.
  !!
  !! Arguments:
  !!   idx1 [in]     -> Index of the first element.
  !!   idx2 [in]     -> Index of the second element.
  !!   array [inout] -> longInt array.
  !!
  pure subroutine swapInArray_longInt(idx1, idx2, array)
    integer(shortInt), intent(in)                 :: idx1, idx2
    integer(longInt), dimension(:), intent(inout) :: array
    integer(longInt)                              :: temp

    temp = array(idx1)
    array(idx1) = array(idx2)
    array(idx2) = temp

  end subroutine swapInArray_longInt

  !! Subroutine 'swapInArray_shortInt'
  !!
  !! Description:
  !!   Swaps two elements within a shortInt array.
  !!
  !! Arguments:
  !!   idx1 [in]     -> Index of the first element.
  !!   idx2 [in]     -> Index of the second element.
  !!   array [inout] -> shortInt array.
  !!
  pure subroutine swapInArray_shortInt(idx1, idx2, array)
    integer(shortInt), intent(in)                  :: idx1, idx2
    integer(shortInt), dimension(:), intent(inout) :: array
    integer(shortInt)                              :: temp

    temp = array(idx1)
    array(idx1) = array(idx2)
    array(idx2) = temp

  end subroutine swapInArray_shortInt

end module genericProcedures