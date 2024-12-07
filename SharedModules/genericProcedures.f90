module genericProcedures
  ! Intrinsic fortran Modules
  use iso_fortran_env, only : compiler_version

  use numPrecision
  use openmp_func, only : ompGetMaxThreads
  use errors_mod,  only : fatalError
  use endfConstants
  use universalVariables

  implicit none

  interface append
    module procedure append_defReal
    module procedure append_shortInt
    module procedure append_shortIntArray
  end interface

  interface countCharacters
    module procedure countCharacters_shortInt
    module procedure countCharacters_longInt
    module procedure countCharacters_defReal
  end interface
  
  interface swap
    module procedure swap_shortInt
    module procedure swap_shortIntArray
    module procedure swap_defReal
    module procedure swap_char_nameLen
    module procedure swap_defReal_defReal
  end interface

  interface quickSort
    module procedure quickSort_shortInt
    module procedure quickSort_defReal
    module procedure quickSort_defReal_defReal
  end interface

  interface hasDuplicates
    module procedure hasDuplicates_char
    module procedure hasDuplicates_shortInt
    module procedure hasDuplicates_defReal
  end interface

  interface hasDuplicatesSorted
    module procedure hasDuplicatesSorted_defReal
  end interface

  interface removeDuplicates
    module procedure removeDuplicates_Char
    module procedure removeDuplicates_shortInt
    module procedure removeDuplicates_Real
  end interface removeDuplicates

  interface removeDuplicatesSorted
    module procedure removeDuplicatesSorted_Real
  end interface

  interface linFind
    module procedure linFind_Char
    module procedure linFind_defReal
    module procedure linFind_defReal_withTOL
    module procedure linFind_shortInt
  end interface

  interface findCommon
    module procedure findCommon_shortInt
  end interface

  interface findDifferent
    module procedure findDifferent_shortInt
  end interface

  interface findDuplicates
    module procedure findDuplicates_Char
  end interface

  interface findDuplicatesSorted
    module procedure findDuplicatesSorted_Real
  end interface

  interface binarySearch
    module procedure binaryFloorIdxClosed_Real
  end interface

  interface endfInterpolate
    module procedure RealReal_endf_interpolate
  end interface

  interface interpolate
    module procedure RealReal_linlin_elemental_interpolate
  end interface

  interface isSorted
    module procedure isSorted_defReal
    module procedure isSorted_shortInt
  end interface

  interface isDescending
    module procedure isDescending_defReal
    module procedure isDescending_shortInt
  end interface

  interface numToChar
    module procedure numToChar_shortInt
    module procedure numToChar_longInt
    module procedure numToChar_defReal
    module procedure numToChar_shortIntArray
    module procedure numToChar_defRealArray
  end interface

  interface concatenate
    module procedure concatenateArrays_Real
  end interface

  interface areEqual
    module procedure areEqual_defReal
    module procedure areEqual_defRealArray
  end interface

contains
  !!
  !! Computes solutions of quadratic equation a * x² + 2 * b * x + c = 0
  !!
  pure function solveQuadratic(a, b, c, delta) result(solutions)
    real(defReal), intent(in)                :: a, b, c, delta
    real(defReal), dimension(:), allocatable :: solutions
    real(defReal)                            :: inverseA
    
    ! Handle special cases first. If a = ZERO, the quadratic equation becomes linear.
    if (areEqual(a, ZERO)) then
      ! If b = ZERO, this is a degenerate case and there are no solutions. Allocate the 
      ! solutions array to 0-size return.
      if (areEqual(b, ZERO)) then
        allocate(solutions(0))
        return

      end if

      ! If reached here, there is only one solution given by d = -c / (2 * b).
      allocate(solutions(1))
      solutions(1) = -HALF * c / b
      return

    end if

    ! If reached here, there are two solutions. Pre-compute 1 / a and compute solutions.
    allocate(solutions(2))
    inverseA = ONE / a
    solutions(1) = -(b + sqrt(delta)) * inverseA
    solutions(2) = -(b - sqrt(delta)) * inverseA

  end function solveQuadratic

  !! Subroutine 'append_shortInt'
  !!
  !! Basic description:
  !!   Adds a shortInt value to a shortInt array. Allocates the array if it is empty.
  !!
  !! Result:
  !!   array -> Array of shortInt integer entries.
  !!
  !! Error:
  !!   None.
  !!
  pure subroutine append_shortInt(array, value)
    integer(shortInt), dimension(:), allocatable, intent(inout) :: array
    integer(shortInt), intent(in)                               :: value
    
    if (.not. allocated(array)) then
      allocate(array(1))
      array = value
      return

    end if
    
    array = [array, value]

  end subroutine append_shortInt
  
  !! Subroutine 'append_shortIntArray'
  !!
  !! Basic description:
  !!   Adds a shortInt array to a shortInt array. Allocates the array if it is empty.
  !!
  !! Result:
  !!   array -> Array of shortInt integer entries.
  !!
  !! Error:
  !!   None.
  !!
  pure subroutine append_shortIntArray(array, values)
    integer(shortInt), dimension(:), allocatable, intent(inout) :: array
    integer(shortInt), dimension(:), intent(in)                 :: values
    
    if (.not. allocated(array)) then
      allocate(array(size(values)))
      array = values
      return

    end if
    
    array = [array, values]

  end subroutine append_shortIntArray

  !! Subroutine 'append_defReal'
  !!
  !! Basic description:
  !!   Adds a defReal entry to a defReal array. Allocates the array if it is empty.
  !!
  !! Result:
  !!   array -> Array of defReal real entries.
  !!
  !! Error:
  !!   None.
  !!
  pure subroutine append_defReal(array, value)
    integer(defReal), dimension(:), allocatable, intent(inout) :: array
    integer(defReal), intent(in)                               :: value
    
    if (.not. allocated(array)) then
      allocate(array(1))
      array = value
      return

    end if
    
    array = [array, value]

  end subroutine append_defReal

  pure function computePyramidCentre(array) result(centre)
    real(defReal), dimension(2, 3), intent(in) :: array
    real(defReal), dimension(3)                :: centre

    centre = FOURTH * sum(array, 1)

  end function computePyramidCentre

  pure function computePyramidVolume(array) result(volume)
    real(defReal), dimension(2, 3), intent(in) :: array
    real(defReal)                              :: volume

    volume = THIRD * abs(dot_product(array(1, :), array(2, :)))

  end function computePyramidVolume

  !! Function 'computeTetrahedronCentre'
  !!
  !! Basic description:
  !!   Computes the centre of the tetrahedron given a 3x4 array containing the
  !!   3-D coordinates of its vertices.
  !!
  !! Arguments:
  !!   array [in] -> 3x4 array containing the 3-D coordinates of the tetrahedron's
  !!                 vertices.
  !!
  !! Result:
  !!   centre     -> 3-D coordinates of the tetrahedron's centre.
  !!
  pure function computeTetrahedronCentre(array) result(centre)
    real(defReal), dimension(4, 3), intent(in) :: array
    real(defReal), dimension(3)                :: centre

    centre = FOURTH * sum(array, 1)

  end function computeTetrahedronCentre

  !! Function 'computeTetrahedronVolume'
  !!
  !! Basic description:
  !!   Computes the volume of the tetrahedron given a 3x4 array containing the
  !!   3-D coordinates of its vertices.
  !!
  !! Arguments:
  !!   array [in] -> 3x4 array containing the 3-D coordinates of the tetrahedron's
  !!                 vertices.
  !!
  !! Result:
  !!   volume     -> Volume of the tetrahedron.
  !!
  pure function computeTetrahedronVolume(array) result(volume)
    real(defReal), dimension(4, 3), intent(in) :: array
    real(defReal)                              :: volume
    real(defReal), dimension(3)                :: A

    A = array(1, :)
    volume = SIXTH * abs(dot_product(crossProduct(array(2, :) - A, array(3, :) - A), array(4, :) - A))

  end function computeTetrahedronVolume

  !! Function 'computeTriangleAra'
  !!
  !! Basic description:
  !!   Computes the area of a triangle from its normal vector.
  !!
  !! Arguments:
  !!   normal [in] -> 3-D components of the triangle's normal vector.
  !!
  pure function computeTriangleArea(normal) result(area)
    real(defReal), dimension(3), intent(in) :: normal
    real(defReal)                           :: area

    area = HALF * norm2(normal)

  end function computeTriangleArea

  !! Function 'computeTriangleCentre'
  !!
  !! Basic description:
  !!   Computes the centre of a triangle given a 3x3 array containing the
  !!   3-D coordinates of its vertices.
  !!
  !! Arguments:
  !!   array [in] -> 3x3 array containing the 3-D coordinates of the triangle's
  !!                 vertices.
  !!
  !! Result:
  !!   centre     -> 3-D coordinates of the triangle's centre.
  !!
  pure function computeTriangleCentre(array) result(centre)
    real(defReal), dimension(3, 3), intent(in) :: array
    real(defReal), dimension(3)                :: centre

    centre = THIRD * sum(array, 1)

  end function computeTriangleCentre

  !! Function 'computeTriangleNormal'
  !!
  !! Basic description:
  !!   Computes the normal vector of a triangle given a 3x3 array containing the
  !!   3-D coordinates of its vertices.
  !!
  !! Arguments:
  !!   array [in] -> 3x3 array containing the 3-D coordinates of the triangle's
  !!                 vertices.
  !!
  !! Result:
  !!   normal     -> 3-D components of the triangle's normal vector.
  !!
  pure function computeTriangleNormal(array) result(normal)
    real(defReal), dimension(3, 3), intent(in) :: array
    real(defReal), dimension(3)                :: normal

    normal = crossProduct(array(2, :) - array(1, :), array(3, :) - array(1, :))

  end function computeTriangleNormal

  !! Function 'countCharacters_shortInt'
  !!
  !! Basic description:
  !!   Counts the number of characters in a shortInt integer.
  !!
  !! Arguments:
  !!   n [in] -> A shortInt integer.
  !!
  !! Result:
  !!   nChars -> Number of characters in the shortInt integer.
  !!
  elemental function countCharacters_shortInt(n) result(nChars)
    integer(shortInt), intent(in) :: n
    integer(shortInt)             :: nChars, absValue

    ! Initialise nChars to 0 or 1 depending on whether n is positive or negative.
    nChars = max(0, sign(1, -n))

    ! Compute the absolute value of n, then increment nChars by 1 for each power of 10 in absValue.
    absValue = abs(n)
    do while (absValue > 0)
      absValue = absValue / 10
      nChars = nChars + 1

    end do

    ! Catch special case where n = 0 (since nChars = 0 in this case).
    nChars = max(nChars, 1)

  end function countCharacters_shortInt

  !! Function 'countCharacters_longInt'
  !!
  !! Basic description:
  !!   Counts the number of characters in a longInt integer.
  !!
  !! Arguments:
  !!   n [in] -> A longInt integer.
  !!
  !! Result:
  !!   nChars -> Number of characters in the longInt integer.
  !!
  elemental function countCharacters_longInt(n) result(nChars)
    integer(longInt), intent(in) :: n
    integer(shortInt)            :: nChars
    integer(longInt)             :: absValue

    ! Initialise nChars to 0 or 1 depending on whether n is positive or negative.
    nChars = 0
    if (n < 0) nChars = 1

    ! Compute the absolute value of n, then increment nChars by 1 for each power of 10 in absValue.
    absValue = abs(n)
    do while (absValue > 0)
      absValue = absValue / 10
      nChars = nChars + 1

    end do

    ! Catch special case where n = 0 (since nChars = 0 in this case).
    nChars = max(nChars, 1)

  end function countCharacters_longInt

  !! Function 'countCharacters_defReal'
  !!
  !! Basic description:
  !!   Counts the number of characters in a defReal real.
  !!
  !! Arguments:
  !!   r [in] -> A defReal real.
  !!
  !! Result:
  !!   nChars -> Number of characters in the defReal real.
  !!
  elemental function countCharacters_defReal(r) result(nChars)
    real(defReal), intent(in) :: r
    integer(shortInt)         :: nChars
    real(defReal)             :: absValue
    real(defReal), parameter  :: lowerBound = 0.1_defReal, upperBound = 1.0E17_defReal

    ! Initialise nChars to 18 or 19 depending on whether r is positive or negative.
    ! Note: here we need to be careful, as to Fortran ZERO = -ZERO, but sign(a, ZERO) /= sign(a, -ZERO)!
    nChars = 18
    if (sign(ONE, r) < ZERO) nChars = nChars + 1

    ! If r = ZERO or r = -ZERO (this is equivalent to Fortran) return early.
    if (r == ZERO) return

    ! Compute the absolute value of r.
    absValue = abs(r)
    
    ! If lowerBound <= absValue < ONE Fortran outputs one more character. Increase nChars 
    ! by 1 and return.
    if (lowerBound <= absValue .and. absValue < ONE) then
      nChars = nChars + 1
      return

    end if

    ! If absValue < lowerBound or upperBound <= absValue then Fortran outputs the number in
    ! scientific notation. Add 5 characters in this case.
    if (absValue < lowerBound .or. upperBound <= absValue) nChars = nChars + 5

  end function countCharacters_defReal

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
  pure function binaryFloorIdxClosed_Real(array, value) result(idx)
    real(defReal), dimension(:), intent(in) :: array
    real(defReal), intent(in)               :: value
    integer(shortInt)                       :: idx, bottom, top, i

    ! Find Top and Bottom Index Array
    bottom = 1
    top = size(array)

    ! Check if the element is in array bounds
    if ( value < array(bottom) .or. value >array(top)) then
      idx = valueOutsideArray
      return
    end if

    do i = 1,70
      !Calculate mid point
      idx = (top + bottom)/2

      ! Termination condition
      if (bottom == idx) return

      ! Binary Step
      if (array(idx) <= value ) then
        bottom = idx
      else
        top = idx
      end if
    end do

    ! Failed to end in 70 steps
    idx = tooManyIter

  end function binaryFloorIdxClosed_Real

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
  !!   -> if value == a(N) returns N-1
  !!   -> if value == a(1) returns 1
  !!
  !! Errors:
  !!   idx == valueOutideArray if value < a(1) .or. value > a(N)
  !!
  pure function linearFloorIdxClosed_Real(array,value) result (idx)
    real(defReal),dimension(:),intent(in) :: array
    real(defReal),intent(in)              :: value
    integer(shortInt)                     :: idx

    if (value > array(size(array)) .or. value < array(1)) then
      idx = valueOutsideArray
      return
    end if

    do idx = size(array)-1,1,-1
      if( array(idx) <= value) return
    end do

  end function linearFloorIdxClosed_Real

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
  !!   -> if value == a(N) returns N-1
  !!   -> if value == a(1) returns 1
  !!
  !! Errors:
  !!   idx == valueOutideArray if value < a(1) .or. value > a(N)
  !!
  function linearFloorIdxClosed_shortInt(array,value) result(idx)
    integer(shortInt),dimension(:),intent(in) :: Array
    integer(shortInt),intent(in)              :: Value
    integer(shortInt)                         :: idx
    character(100),parameter                  :: Here='linearFloorIdxClosed_shortInt (genericProcedures.f90)'

    ! Check if the value is above the bounds of an array
    if ( Value >= array(size(array))) call fatalError(Here,'Value is above upper bound of the array')

    do idx=size(array),1,-1
      if ( array(idx) <= value ) return
    end do

    call fatalError(Here,'Value is below lower bound of the array')

  end function linearFloorIdxClosed_shortInt

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
  !!   If value <= a(a) returns 1.
  !!
  !! Errors:
  !!   idx == valueOutideArray if value > a(N)
  !!
  pure function linearCeilingIdxOpen_shortInt(array,value) result(idx)
    integer(shortInt),dimension(:),intent(in) :: Array
    integer(shortInt),intent(in)              :: Value
    integer(shortInt)                         :: idx
    character(100),parameter                  :: Here='linearCeilingIdxOpen_shortInt (genericProcedures.f90)'

    do idx=1,size(array)
      if ( array(idx) >= value ) return
    end do

    ! Value is larger than the upper bound of the array
    idx = valueOutsideArray

  end function linearCeilingIdxOpen_shortInt

  !!
  !! Subroutine that checks whether there was an error during search and returns approperiate
  !! message.
  !!
  !! TODO: IT IS POINTLESS AND WILL BE REMOVED.
  !! ERROR HANDLING FOR SEARCHES SHOULD BE DONE BY CLIENT
  !!
  !! Although it can be replaced by ERROR_CODE -> ERROR_MESSAGE translation function to
  !! give nice strings instead of raw integers in the error messages...
  !!
  subroutine searchError(idx,Here)
    integer(shortInt),intent(in)  :: idx
    character(*),intent(in)       :: Here

    if (idx < 0) then
      select case (idx)
        case (valueOutsideArray)
          call fatalError(Here,'The requested value was outide the array bounds')

        case (tooManyIter)
          call fatalError(Here,'Search did not terminate in hardcoded number of iterations')

        case (targetNotFound)
          call fatalError(Here,'Search failed to find target in the array')

        case default
          call fatalError(Here,'Search returned unknown error flag (negative index)')

      end select
    end if
  end subroutine searchError

  !!
  !! Open "File" for reading under with "unitNum" reference
  !!
  subroutine openToRead(unitNum,File)
    integer(shortInt), intent(in)  :: unitNum
    character(*), intent(in)       :: File
    integer(shortInt)              :: errorStat
    character(99)                  :: errorMsg

        open ( unit   = unitNum,   &
               file   = File,      &
               status = "old",     &
               action = "read",    &
               iostat = errorStat, &
               iomsg  = errorMsg)

        if (errorStat > 0) call fatalError('openToRead subroutine (genericProcedures.f90)', &
                                           errorMsg )
  end subroutine openToRead


  !!
  !! Removes duplicates from an unsorted intArray
  !! Does not preserve the order
  !!
  pure function removeDuplicates_shortInt(intArray) result(out)
    integer(shortInt), dimension(:), intent(in) :: intArray
    integer(shortInt), dimension(:),allocatable :: out
    integer(shortInt),dimension(size(intArray)) :: array
    integer(shortInt)                           :: i, j, N, test, idx

    N = size(intArray)

    ! If provided array is of size 0, return size 0 array
    if (N == 0 ) then
      allocate(out(0))
      return
    end if

    array(1) = intArray(1)    ! Copy first element
    j = 2                     ! Set new empty index
    do i=2,N
      test = intArray(i)
      idx = linFind(intArray(1:i-1),test)  ! Check if element is present before i in the array
      if  (idx == targetNotFound ) then    ! If it isn't copy it to array without duplicates
        array(j) = test
        j = j + 1

      end if
    end do

    out = array(1:j-1)  ! Allocate output

  end function removeDuplicates_shortInt

  !!
  !! Removes duplicates from an unsorted realArray
  !! Does not preserve the order
  !!
  pure function removeDuplicates_Real(realArray) result(out)
    real(defReal), dimension(:), intent(in)   :: realArray
    real(defReal), dimension(:),allocatable   :: out
    real(defReal), dimension(size(realArray)) :: array
    real(defReal)                             :: test
    integer(shortInt)                         :: i, j, N

    N = size(realArray)

    ! If provided array is of size 0, return size 0 array
    if (N == 0 ) then
      allocate(out(0))
      return
    end if

    array(1) = realArray(1)    ! Copy first element
    j = 2                     ! Set new empty index
    do i=2,N
      test = realArray(i)
       if (any(array(1:j-1) == test)) cycle   ! Skip cycle if a duplicate is found
        array(j) = test
        j = j + 1
    end do

    out = array(1:j-1)  ! Allocate output

  end function removeDuplicates_Real

  !!
  !! Function that removes duplicates from input character array. It returns array of equal or
  !! smaller size. Unfortunatly Fortran requires output character to have specified length. Length
  !! of 100 is hardcoded at the moment. Function returns fatal error if input characters are of
  !! length greater then 100
  !!
  function removeDuplicates_Char(charArray) result(out)
    character(nameLen),dimension(:),intent(in)   :: charArray
    character(nameLen),dimension(:),allocatable  :: out
    logical(defBool),dimension(:),allocatable    :: unique
    integer(shortInt)                            :: i,j

    if (len(charArray) > len(out)) call fatalError('removeDuplicates_Char (genericProcedures.f90)',&
                                                   'Maximum length of input character is 100 ')
    if (size(charArray) == 1) then
      out = charArray
    else
      allocate(unique(size(charArray)))
      unique = .true.
      ! For every element search if it matches any previous element. Change uniqe to false if it is
      ! repeted.
        do i = 1,size(charArray)
          search: &
          do j = 1, i-1
            if( trim(charArray(i)) == trim(charArray(j)) ) then
              unique(i) = .false.
              exit search
            end if
          end do search
        end do
        ! Select elements from charArray for which unique == .true.
        out = pack(charArray, unique)
    end if
  end function removeDuplicates_Char

  !!
  !! Removes duplicates from a sorted realArray
  !! Preserves the order
  !!
  pure function removeDuplicatesSorted_Real(realArray) result(out)
    real(defReal), dimension(:), intent(in)   :: realArray
    real(defReal), dimension(:),allocatable   :: out
    real(defReal), dimension(size(realArray)) :: array
    real(defReal)                             :: test
    integer(shortInt)                         :: i, j, N

    N = size(realArray)

    ! If provided array is of size 0, return size 0 array
    if (N == 0 ) then
      allocate(out(0))
      return
    end if

    array(1) = realArray(1)    ! Copy first element
    j = 2                     ! Set new empty index
    do i=2,N
      test = realArray(i)
       if (array(j-1) == test) cycle   ! Skip cycle if a duplicate is found
        array(j) = test
        j = j + 1
    end do

    out = array(1:j-1)  ! Allocate output

  end function removeDuplicatesSorted_Real

  !!
  !! Find duplicates from a sorted realArray
  !! Returns an array with the index of the repeated element
  !!
  pure function findDuplicatesSorted_Real(realArray) result(out)
    real(defReal), dimension(:), intent(in)       :: realArray
    integer(shortInt), dimension(:),allocatable   :: out
    integer(shortInt), dimension(size(realArray)) :: array
    real(defReal)                                 :: test
    integer(shortInt)                             :: i, j, N

    N = size(realArray)

    ! If provided array is of size 0, return size 0 array
    if (N == 0 ) then
      allocate(out(0))
      return
    end if

    j = 1                    ! Set new empty index
    do i=2,N
      test = realArray(i)
       if (realArray(i-1) /= test) cycle   ! Skip cycle if a duplicate is not found
        array(j) = i
        j = j + 1
    end do

    out = array(1:j-1)  ! Allocate output

  end function findDuplicatesSorted_Real

  !!
  !! Function that finds duplicates in array of characters. Returns array that contains repeted
  !! element. Unfortunatly Fortran requires output character to have specified length. Length
  !! of 100 is hardcoded at the moment. Function returns fatal error if input characters are of
  !! length greater then 100
  !!
  function findDuplicates_Char(charArray) result(out)
    character(nameLen),dimension(:),intent(in)  :: charArray
    character(nameLen),dimension(:),allocatable :: out
    logical(defBool),dimension(:),allocatable   :: unique
    integer(shortInt)                           :: i,j

    if (len(charArray) > len(out)) call fatalError('removeDuplicates_Char (genericProcedures.f90)',&
                                                   'Maximum length of input character is 100 ')
    if (size(charArray) == 1) then
      out = charArray
    else
      allocate(unique(size(charArray)))
      unique = .true.
      ! For every element search if it matches any previous element. Change uniqe to false if it is
      ! repeted.
        do i = 1,size(charArray)
          search: &
          do j = 1, i-1
            if( trim(charArray(i)) == trim(charArray(j)) ) then
              unique(i) = .false.
              exit search
            end if
          end do search
        end do
        ! Select elements from charArray for which unique == .true.
        out = pack(charArray, .not.unique)
        out = removeDuplicates(out)
    end if
  end function findDuplicates_Char

  !! Function 'findCommon_shortInt'
  !!
  !! Basic description:
  !!   Find common entries between two shortInt arrays.
  !!
  !! Result:
  !!   commons -> shortInt array containing the common entries between the two supplied arrays.
  !!
  !! Error:
  !!   None.
  !!
  pure function findCommon_shortInt(array1, array2) result(commons)
    integer(shortInt), dimension(:), intent(in)  :: array1, array2
    integer(shortInt), dimension(:), allocatable :: commons
    integer(shortInt), dimension(size(array1))   :: array1Copy
    integer(shortInt), dimension(size(array2))   :: array2Copy
    integer(shortInt)                            :: i, j
    
    ! Initialise common array to 0-size, copy the two input arrays and sort them.
    allocate(commons(0))
    array1Copy = array1
    array2Copy = array2
    call quickSort(array1Copy)
    call quickSort(array2Copy)

    ! Initialise i = 1 and j = 1 then search for all entries in each array.
    i = 1
    j = 1
    do while (i <= size(array1Copy) .and. j <= size(array2Copy))
      if (array1Copy(i) == array2Copy(j)) then
        call append(commons, array1Copy(i))
        i = i + 1
        j = j + 1

      else if (array1Copy(i) < array2Copy(j)) then
        i = i + 1

      else
        j = j + 1

      end if

    end do

  end function findCommon_shortInt

  !! Function 'findDifferent_shortInt'
  !!
  !! Basic description:
  !!   Returns the element of a size-2 shortInt array not equal to value. Useful for mesh tracking
  !!   where we need to use lots of connectivity information.
  !!
  !! Arguments:
  !!   array [in] -> A shortInt array of length 2.
  !!   value [in] -> Value to exclude.
  !!
  !! Result:
  !!   different -> The element in the array not equal to value.
  !!
  pure function findDifferent_shortInt(array, value) result(different)
    integer(shortInt), dimension(2), intent(in) :: array
    integer(shortInt), intent(in)               :: value
    integer(shortInt)                           :: different
    
    if (array(1) == value) then
      different = array(2)

    else
      different = array(1)

    end if
    
  end function findDifferent_shortInt

  !!
  !! Searches linearly for the occurance of target in charArray. Removes left blanks.
  !! Following Errors can occur:
  !! targetNotFound -> target is not present in the array
  !!
  pure function linFind_Char(charArray,target) result(idx)
    character(*),dimension(:),intent(in) :: charArray
    character(*),intent(in)              :: target
    integer(shortInt)                    :: idx

    do idx=1,size(charArray)
     ! if( trim(charArray(idx)) == trim(target) ) return
      if( adjustl(charArray(idx)) == adjustl(target)) return
    end do
    idx = targetNotFound

  end function linFind_Char

  !!
  !! Searches linearly for the occurance of target in defRealArray. Following Errors can occur
  !! valueOutsideArray -> target is not present in the array
  !!
  pure function linFind_defReal(defRealArray,target) result (idx)
    real(defReal), dimension(:), intent(in)  :: defRealArray
    real(defReal), intent(in)                :: target
    integer(shortInt)                        :: idx

    do idx=1,size(defRealArray)
      if (defRealArray(idx) == target) return
    end do
    idx = targetNotFound
  end function linFind_defReal

  !!
  !! Searches linearly for the occurance of target in defRealArray.
  !!
  !! Uses relative tolerance of TOL
  !!
  !! Args:
  !!   defRealArray [in] -> array to search
  !!   target [in]       -> traget to search for
  !!   tol [in]          -> relative tolerance (+ve)
  !!
  pure function linFind_defReal_withTOL(defRealArray, target, tol) result (idx)
    real(defReal), dimension(:), intent(in)  :: defRealArray
    real(defReal), intent(in)                :: target
    real(defReal), intent(in)                :: tol
    integer(shortInt)                        :: idx
    real(defReal)                            :: inv_T

    inv_T = ONE/target

    do idx=1,size(defRealArray)
      if (abs(defRealArray(idx) * inv_T - ONE) < tol) return
    end do
    idx = targetNotFound
  end function linFind_defReal_withTOL

  !!
  !! Searches linearly for the occurance of target in shortIntArray. Following Errors can occur
  !! valueOutsideArray -> target is not present in the array
  !!
  pure function linFind_shortInt(shortIntArray,target) result (idx)
    integer(shortInt), dimension(:), intent(in)  :: shortIntArray
    integer(shortInt), intent(in)                :: target
    integer(shortInt)                            :: idx

    do idx=1,size(shortIntArray)
      if (shortIntArray(idx) == target) return
    end do
    idx = targetNotFound
  end function linFind_shortInt

  !!
  !! Given 2 arrays, it outputs a third array which is a concatenation of them
  !!
  pure function concatenateArrays_Real(array1,array2) result(out)
    real(defReal),dimension(:),intent(in)              :: array1
    real(defReal),dimension(:),intent(in)              :: array2
    real(defReal),dimension(size(array1)+size(array2)) :: out
    integer(shortInt)                                  :: N, M

    N = size(array1)
    M = size(array2)
    out(1:N) = array1
    out(N+1:N+M) = array2

  end function concatenateArrays_Real

  !!
  !! Returns .true. if two floating point numbers are equal. 
  !!
  !! Due to floating point artihmetic and rounding-off errors being slightly different 
  !! across different architectures (eg, Intel vs ARM), it is necessary to use some small 
  !! tolerances to assert equality between two floating point numbers. These tolerances 
  !! are specified in numPrecision.f90.
  !!
  elemental function areEqual_defReal(a, b) result(equal)
    real(defReal), intent(in) :: a, b
    logical(defBool)          :: equal
    real(defReal)             :: absDiff

    ! Initialise equal = .true. and check for perfect (to the bit) equality, since we can 
    ! return early in this case.
    equal = .true.
    if (a == b) return

    ! Compute the absolute value of the difference between the two floating point numbers.
    ! Note that if a and b are both very large and of opposite signs this can cause overflow.
    absDiff = abs(a - b)

    ! Check if absDiff is less than some absolute very small tolerance first and return if yes.
    if (absDiff < floatTol) return

    ! Check if a and b are within some small relative tolerance of each other and return if
    ! yes. Note that if a and b are both very small numbers, then multiplying by a small
    ! tolerance can cause underflow. This is why we check absolute tolerance first.
    if (absDiff < max(abs(a), abs(b)) * FP_REL_TOL) return

    ! If reached here, a and b are not within absolute or relative tolerance of each other.
    ! update equal = .false.
    equal = .false.

  end function areEqual_defReal

  !!
  !! Returns .true. if all floating point numbers in an array are equal to a given value. 
  !!
  !! Due to floating point artihmetic and rounding-off errors being slightly different 
  !! across different architectures (eg, Intel vs ARM), it is necessary to use some small 
  !! tolerances to assert equality between two floating point numbers. These tolerances 
  !! are specified in numPrecision.f90.
  !!
  pure function areEqual_defRealArray(array, b) result(equal)
    real(defReal), dimension(:), intent(in)     :: array
    real(defReal), intent(in)                   :: b
    logical(defBool)                            :: equal
    integer(shortInt)                           :: i
    real(defReal)                               :: a, absDiff

    ! Initialise equal = .true. and loop over all element in the array.
    equal = .true.
    do i = 1, size(array)
      ! Retrieve current element of the array and check for perfect (to the bit) equality. 
      ! Cycle to the next element if yes.
      a = array(i)
      if (a == b) cycle

      ! Compute the absolute value of the difference between the two floating point numbers.
      ! Note that if a and b are both very large and of opposite signs this can cause overflow.
      absDiff = abs(a - b)

      ! Check if absDiff is less than some absolute very small tolerance first and cycle if yes.
      if (absDiff < floatTol) cycle

      ! Check if a and b are within some small relative tolerance of each other and cycle if
      ! yes. Note that if a and b are both very small numbers, then multiplying by a small
      ! tolerance can cause underflow. This is why we check absolute tolerance first.
      if (absDiff < max(abs(a), abs(b)) * FP_REL_TOL) cycle

      ! If reached here, a and b are not within absolute or relative tolerance of each other.
      ! update equal = .false. and return.
      equal = .false.
      return

    end do

  end function areEqual_defRealArray

  !!
  !! Concatenate strings from an array into a single long character (tape). Asjusts left and trims
  !! elements of char Array. Adds a blank at the end of a line
  !!
  pure function arrayConcat(charArray) result(out)
    character(*),dimension(:),intent(in)       :: charArray
    character(:),allocatable                   :: out
    integer(shortInt)                          :: elementLen, trimLen , i

    ! Find total length of elements of charArray after adjusting left and trimming
    trimLen=0
    do i=1,size(charArray)
      elementLen =  len( trim( adjustl( charArray(i))))
      trimLen = trimLen + elementLen
    end do

    ! Allocate output array
    i = trimLen+size(charArray)
    allocate(character(i) :: out)

    ! Make shure output tape is empty
    out = ''

    ! Write elements of the array to output tape
    do i=1,size(charArray)
      out = out // trim( adjustl(charArray(i))) // ' '
    end do

  end function arrayConcat

  !!
  !! Function that searches counts all occurences of a "symbol" in a "string"
  !!
  pure function countSymbol(string,symbol) result(num)
    character(*),intent(in)  :: string
    character(1),intent(in)  :: symbol
    integer(shortInt)        :: num
    integer(shortInt)        :: start, end, pos

    start = 1
    end   = len(string)
    num   = 0

    do
      pos = index(string(start:end),symbol)
      if (pos == 0) return ! No more occurences of symbol in "string"
      num = num + 1
      start = start + pos ! pos is returned relative to start
    end do

  end function countSymbol

  !!
  !! Goes through the string and adds +1 to balance for each leftS and -1 for each rightS. It
  !! terminates and returns -1 when balance becomes -ve.
  !!
  pure function symbolBalance(str,leftS,rightS) result (balance)
    character(*),intent(in)       :: str
    character(1),intent(in)       :: leftS
    character(1),intent(in)       :: rightS
    integer(shortInt)             :: balance
    integer(shortInt)             :: i

    balance = 0
    do i=1,len(str)

      if(str(i:i) == leftS) then
        balance = balance + 1

      elseif(str(i:i) == rightS) then
        balance = balance -1

      end if

      if (balance < 0) return

    end do

  end function symbolBalance

  !!
  !! Finds index of next character in signs in the string
  !!
  pure function indexOfNext(signs,string) result (idx)
    character(1),dimension(:),intent(in)     :: signs
    character(*), intent(in)                 :: string
    integer(shortInt)                        :: idx
    integer(shortInt),dimension(size(signs)) :: temp_idx

    temp_idx = index(string,signs)
    idx = minval(temp_idx,temp_idx > 0)
    if (idx == huge(temp_idx)) idx = 0

  end function indexOfNext

  !!
  !! Replaces all repeated blanks in string with a single blank
  !!
  pure subroutine compressBlanks(string)
    character(*), intent(inout)     :: string
    character(len(string))          :: stringCopy
    integer(shortInt)               :: i, j
    logical(defBool)                :: lastBlank


    lastBlank = .false.
    stringCopy = ''
    j = 1

    do i=1,len(string)
      if (lastBlank) then
        if (string(i:i) /= " ") then
          lastBlank = .false.
          stringCopy(j:j) = string(i:i)
          j=j+1
        end if

      else
        stringCopy(j:j) = string(i:i)
        j=j+1
        if (string(i:i) == " ") then
          lastBlank = .true.
        endif

      end if
    end do
    string = stringCopy

  end subroutine compressBlanks

  !!
  !! Replaces all symbols "oldS" with "newS" in string
  !!
  pure subroutine replaceChar(string,oldS,newS)
    character(*), intent(inout) :: string
    character(1), intent(in)    :: oldS
    character(1), intent(in)    :: newS
    integer(shortInt)           :: i

    do i=1,len(string)
      if(string(i:i) == oldS) string(i:i) = newS
    end do

  end subroutine replaceChar

  !!
  !! Compares strings for equality. Ignores leading blanks.
  !!
  elemental function charCmp(char1, char2) result(areEqual)
    character(*), intent(in)  :: char1
    character(*), intent(in)  :: char2
    logical(defBool)          :: areEqual

    areEqual = (trim(adjustl(char1)) == trim(adjustl(char2)))

  end function charCmp

  !!
  !! Perform linear interpolation between defReals
  !!
  elemental function RealReal_linlin_elemental_interpolate(xMin,xMax,yMin,yMax,x) result(y)
    real(defReal), intent(in) :: xMin, xMax, yMin, yMax, x
    real(defReal)             :: y
    real(defReal)             :: interFactor

    interFactor = (x-xMin)/(xMax-xMin)
    y = yMax * interFactor + (1-interFactor)*yMin

  end function RealReal_linlin_elemental_interpolate

  !!
  !! Perform one of ENDF defined interpolation types
  !!
  function RealReal_endf_interpolate(xMin,xMax,yMin,yMax,x,endfNum) result(y)
    real(defReal), intent(in)     :: xMin, xMax, yMin, yMax, x
    integer(shortInt), intent(in) :: endfNum
    real(defReal)                 :: y
    character(100),parameter      :: Here='RealReal_endf_interpolate (genericProcedures.f90)'

    select case (endfNum) ! Naming Convention for ENDF interpolation (inY-inX) i.e. log-lin => logarithmic in y; linear in x
      case (histogramInterpolation)
        y = yMin

      case (linLinInterpolation)
        y = interpolate(xMin,xMax,yMin,yMax,x)

      case (linLogInterpolation)
        y = interpolate(log(xMin),log(xMax),yMin,yMax,log(x))

      case (logLinInterpolation)
        y = interpolate(xMin,xMax,log(yMin),log(yMax),x)
        y = exp(y)

      case (loglogInterpolation)
        y = interpolate(log(xMin),log(xMax),log(yMin),log(yMax),log(x))
        y = exp(y)

      case (chargedParticleInterpolation)
        ! Not implemented
        call fatalError(Here, 'ENDF interpolation law for charged Particles is not implemented')

      case default
        call fatalError(Here, 'Unknown ENDF interpolation number')
    end select

  end function RealReal_endf_interpolate

  !!
  !! Checks if the provided float is an integer. It may not be rebust and requires further
  !! testing. It uses the fact that ceiling and floor of an integer are the same.
  !!
  elemental function isInteger(float) result (isIt)
    real(defReal), intent(in) :: float
    logical(defBool)          :: isIt

    ! It is Fortran 2008 Standard Compliant. Really!
    isIt = (floor(float,longInt) == ceiling(float,longInt))

  end function isInteger

  !!
  !! Function that check if the array is sorted in ascending order (a(i) >= a(i-1) for all i).
  !!
  pure function isSorted_defReal(array) result (isIt)
    real(defReal),dimension(:),intent(in) :: array
    logical(defBool)                      :: isIt
    integer(shortInt)                     :: i

    do i=2,size(array)
      if (array(i) < array(i-1)) then
        isIt = .false.
        return
      end if
    end do

    isIt = .true.

  end function isSorted_defReal

  !!
  !! Function that check if the array is sorted in ascending order (a(i) >= a(i-1) for all i).
  !!
  pure function isSorted_shortInt(array) result (isIt)
    integer(shortInt),dimension(:),intent(in) :: array
    logical(defBool)                          :: isIt
    integer(shortInt)                         :: i

    do i=2,size(array)
      if (array(i) < array(i-1)) then
        isIt = .false.
        return
      end if
    end do

    isIt = .true.

  end function isSorted_shortInt

  !!
  !! Function that check if the array is sorted in descending order (a(i) <= a(i-1) for all i).
  !!
  pure function isDescending_defReal(array) result (isIt)
    real(defReal),dimension(:),intent(in) :: array
    logical(defBool)                      :: isIt
    integer(shortInt)                     :: i

    do i=2,size(array)
      if (array(i) > array(i-1)) then
        isIt = .false.
        return
      end if
    end do

    isIt = .true.

  end function isDescending_defReal

  !!
  !! Function that check if the array is sorted in descending order (a(i) <= a(i-1) for all i).
  !!
  pure function isDescending_shortInt(array) result (isIt)
    integer(shortInt),dimension(:),intent(in) :: array
    logical(defBool)                          :: isIt
    integer(shortInt)                         :: i

    do i=2,size(array)
      if (array(i) > array(i-1)) then
        isIt = .false.
        return
      end if
    end do

    isIt = .true.

  end function isDescending_shortInt

  !!
  !! Convert shortInt to character
  !!
  function numToChar_shortInt(x) result(c)
    integer(shortInt),intent(in)  :: x
    character(countCharacters(x)) :: c

    write(c, '(I0)') x

  end function numToChar_shortInt

  !!
  !! Convert shortInt array to character
  !!
  function numToChar_shortIntArray(x) result(c)
    integer(shortInt), dimension(:), intent(in) :: x
    character(:), allocatable                   :: c, tempChar
    integer(shortInt)                           :: nIntegers, n, i

    c = ''
    nIntegers = size(x)
    if (nIntegers == 0) return

    n = x(1)
    allocate(character(countCharacters(n)) :: tempChar)
    write(tempChar,'(I0)') n
    c = tempChar
    deallocate(tempChar)

    do i = 2, nIntegers
      n = x(i)
      allocate(character(countCharacters(n)) :: tempChar)
      write(tempChar,'(I0)') n
      c = c//' '//tempChar
      deallocate(tempChar)

    end do

  end function numToChar_shortIntArray

  !!
  !! Convert longInt to character
  !!
  function numToChar_longInt(x) result(c)
    integer(longInt), intent(in)  :: x
    character(countCharacters(x)) :: c

    write(c, '(I0)') x

  end function numToChar_longInt

  !!
  !! Convert defReal to character
  !!
  function numToChar_defReal(x) result(c)
    real(defReal),intent(in)      :: x
    character(countCharacters(x)) :: c

    write(c, *) x

  end function numToChar_defReal

  !!
  !! Convert defReal array to character
  !!
  function numToChar_defRealArray(x) result(c)
    real(defReal), dimension(:),intent(in) :: x
    character(:), allocatable              :: c, tempChar
    integer(shortInt)                      :: nReals, i
    real(defReal)                          :: r

    c = ''
    nReals = size(x)
    if (nReals == 0) return

    r = x(1)
    allocate(character(countCharacters(r)) :: tempChar)
    write(tempChar, *) r
    c = tempChar
    deallocate(tempChar)

    do i = 2, nReals
      r = x(i)
      allocate(character(countCharacters(r)) :: tempChar)
      write(tempChar, *) r
      c = c//' '//tempChar
      deallocate(tempChar)

    end do

  end function numToChar_defRealArray

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
      read(str,'(I20)', IOSTAT=state) i
      error = state /= 0
    else
      read(str,'(I20)') i
    end if

  end function charToInt

  !!
  !! Subroutine takes a normilised direction vector dir and rotates it by cosine of a polar angle
  !! mu and azimuthal angle phi (in radians).
  !! Procedure will produce incorrect results WITHOUT error message if dir is not normalised
  !!
  pure function rotateVector(dir, mu, phi) result(newDir)
    real(defReal), dimension(3), intent(in)    :: dir
    real(defReal), intent(in)                  :: mu
    real(defReal), intent(in)                  :: phi
    real(defReal), dimension(3)                :: newDir
    real(defReal)                              :: u,v,w
    real(defReal)                              :: sinPol, cosPol, A, B

    ! Precalculate cosine and sine of polar angle
    sinPol = sin(phi)
    cosPol = cos(phi)

    ! Load Cartesian components of direction.
    u = dir(1)
    v = dir(2)
    w = dir(3)

    ! Perform standard roatation. Note that indexes are parametrised
    A = sqrt(max(ZERO, ONE - mu*mu))
    B = sqrt(max(ZERO, ONE - w*w  ))


    if ( B > 1E-8) then
      newDir(1) = mu * u + A * (u*w*cosPol - v * sinPol) / B
      newDir(2) = mu * v + A * (v*w*cosPol + u * sinPol) / B
      newDir(3) = mu * w - A * B * cosPol

    else
      B = sqrt(max(ZERO, ONE - v*v))
      newDir(1) = mu * u + A *(u*v*cosPol + w * sinPol) / B
      newDir(2) = mu * v - A * B * cosPol
      newDir(3) = mu * w + A* (v*w*cosPol - u * sinPol) /B

    end if

    newDir = newDir / norm2(newDir)

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
  subroutine rotationMatrix(matrix, phi, theta, psi)
    real(defReal), dimension(3,3), intent(out) :: matrix
    real(defReal), intent(in)                  :: phi
    real(defReal), intent(in)                  :: theta
    real(defReal), intent(in)                  :: psi
    real(defReal) :: sin_phi, cos_phi, sin_t, cos_t, sin_psi, cos_psi, conv
    character(100), parameter :: Here = 'rotationMatrix (genericProcedures.f90)'

    ! Check input data
    if (phi < ZERO .or. phi >= 360.0_defReal) then
      call fatalError(Here, 'Angle phi must be in <0;360). Is: '//numToChar(phi))

    else if (theta < ZERO .or. theta > 180.0_defReal) then
      call fatalError(Here, 'Angle theta must be in <0;180>. Is: '//numToChar(theta))

    else if (psi < ZERO .or. psi >= 360.0_defReal) then
      call fatalError(Here, 'Angle psi must be in <0;180>. Is: '//numToChar(theta))

    end if

    ! Evaluate trigonometric functions
    conv = TWO_PI / 360.0_defReal
    sin_phi = sin(phi * conv)
    cos_phi = cos(phi * conv)
    sin_t   = sin(theta * conv)
    cos_t   = cos(theta * conv)
    sin_psi = sin(psi * conv)
    cos_psi = cos(psi * conv)

    ! Assign matrix elemets
    matrix(1,1) = cos_psi * cos_phi - cos_t * sin_phi * sin_psi
    matrix(1,2) = cos_psi * sin_phi + cos_t * cos_phi * sin_psi
    matrix(1,3) = sin_psi * sin_phi

    matrix(2,1) = -sin_psi * cos_phi - cos_t * sin_phi * cos_psi
    matrix(2,2) = -sin_psi * sin_phi + cos_t * cos_phi * cos_psi
    matrix(2,3) = cos_psi * sin_t

    matrix(3,1) = sin_t * sin_phi
    matrix(3,2) = -sin_t * cos_phi
    matrix(3,3) = cos_t

  end subroutine rotationMatrix

  !!
  !! Cross product for 3D vectors
  !!
  pure function crossProduct(a, b) result(c)
    real(defReal), dimension(3), intent(in) :: a, b
    real(defReal), dimension(3)             :: c

    c = [a(2)*b(3) - a(3)*b(2), &
         a(3)*b(1) - a(1)*b(3), &
         a(1)*b(2) - a(2)*b(1)]

  end function crossProduct

  !!
  !! Return true if key is in the array
  !!
  !! Args:
  !!   array [in] -> Array of data
  !!   key [in]   -> Required key
  !!
  !! Result:
  !!   True if key is in array. False otherwise
  !!
  pure function isIn(array, key)
    integer(shortInt), dimension(:), intent(in) :: array
    integer(shortInt), intent(in)               :: key
    logical(defBool)                            :: isIn

    isIn =  targetNotFound /= linFind(array, key)

  end function isIn

  !!
  !! Returns .true. if string array contains duplicates.
  !!
  pure function hasDuplicates_char(array) result(doesIt)
    character(nameLen), dimension(:), intent(in) :: array
    integer(shortInt)                            :: nStrings, i, j
    logical(defBool)                             :: doesIt
    character(:), allocatable                    :: testString
    
    ! Initialise doesIt = .true. and compute number of strings in the array.
    doesIt = .true.
    nStrings = size(array)

    ! Loop from first till second last element in the array.
    do i = 1, nStrings - 1
      ! Create adjusted and trimmed current string and check for duplicates
      ! in remaining strings of the array.
      testString = adjustl(trim(array(i)))
      do j = i + 1, nStrings
        if (testString == adjustl(trim(array(j)))) return

      end do

    end do
    
    ! If reached here, update doesIt = .false.
    doesIt = .false.

  end function hasDuplicates_char

  !!
  !! Returns true if array contains duplicates
  !!
  pure function hasDuplicates_shortInt(array) result(doesIt)
    integer(shortInt), dimension(:), intent(in) :: array
    logical(defBool)                            :: doesIt
    integer(shortInt), dimension(size(array))   :: arrayCopy
    integer(shortInt)                           :: i

    ! Copy and sort array
    arrayCopy = array
    call quickSort(arrayCopy)

    ! Search through the array looking for duplicates
    doesIt = .false.
    do i=2,size(array)
      doesIt = doesIt .or. arrayCopy(i) == arrayCopy(i-1)

    end do

  end function hasDuplicates_shortInt

  !!
  !! Returns tue if array contains duplicates
  !!
  pure function hasDuplicates_defReal(array) result(doesIt)
    real(defReal), dimension(:), intent(in) :: array
    logical(defBool)                        :: doesIt
    real(defReal), dimension(size(array))   :: arrayCopy
    integer(shortInt)                       :: i

    ! Copy and sort array
    arrayCopy = array
    call quickSort(arrayCopy)

    ! Search through the array looking for duplicates
    doesIt = .false.
    do i=2,size(array)
      doesIt = doesIt .or. arrayCopy(i) == arrayCopy(i-1)

    end do

  end function hasDuplicates_defReal

  !!
  !! Returns tue if array contains duplicates
  !! Only for already sorted arrays
  !!
  pure function hasDuplicatesSorted_defReal(array) result(doesIt)
    real(defReal), dimension(:), intent(in) :: array
    logical(defBool)                        :: doesIt
    integer(shortInt)                       :: i

    ! Search through the array looking for duplicates
    doesIt = .false.
    do i=2,size(array)
      doesIt = doesIt .or. array(i) == array(i-1)
    end do

  end function hasDuplicatesSorted_defReal

  !!
  !! Quicksort for integer array
  !!
  recursive pure subroutine quickSort_shortInt(array)
    integer(shortInt), dimension(:), intent(inout) :: array
    integer(shortInt)                              :: pivot
    integer(shortInt)                              :: i, maxSmall

    if (size(array) > 1 ) then
      ! Set a pivot to the rightmost element
      pivot = size(array)

      ! Move all elements <= pivot to the LHS of the pivot
      ! Find position of the pivot in the array at the end (maxSmall)
      maxSmall = 0
      do i=1,size(array)

        if( array(i) <= array(pivot)) then
          maxSmall = maxSmall + 1
          call swap(array(i),array(maxSmall))
        end if
      end do

      ! Recursivly sort the sub arrays divided by the pivot
      call quickSort(array(1:maxSmall-1))
      call quickSort(array(maxSmall+1:size(array)))
    end if

  end subroutine quickSort_shortInt

  !!
  !! Quicksort for real array
  !!
  recursive pure subroutine quickSort_defReal(array)
    real(defReal), dimension(:), intent(inout) :: array
    integer(shortInt)                          :: i, maxSmall, pivot

    if (size(array) > 1 ) then
      ! Set a pivot to the rightmost element
      pivot = size(array)

      ! Move all elements <= pivot to the LHS of the pivot
      ! Find position of the pivot in the array at the end (maxSmall)
      maxSmall = 0
      do i=1,size(array)

        if( array(i) <= array(pivot)) then
          maxSmall = maxSmall + 1
          call swap(array(i),array(maxSmall))
        end if
      end do

      ! Recursivly sort the sub arrays divided by the pivot
      call quickSort(array(1:maxSmall-1))
      call quickSort(array(maxSmall+1:size(array)))
    end if

  end subroutine quickSort_defReal

  !!
  !! Quicksort for array1 and array2 of reals by array1
  !! If size(array2) > size(array1) ignores extra elements (does not sort them)
  !! If size(array1) > size(array2) behaviour is undefined.
  !!
  recursive subroutine quickSort_defReal_defReal(array1, array2)
    real(defReal), dimension(:), intent(inout) :: array1
    real(defReal), dimension(:), intent(inout) :: array2
    integer(shortInt)                          :: i, maxSmall, pivot
    character(100),parameter :: Here = 'quickSort_defReal_defReal (genericProcdures.f90)'

    if(size(array1) /= size(array2)) then
      call fatalError(Here,'Arrays have diffrent size!')
    end if

    if (size(array1) > 1 ) then
      ! Set a pivot to the rightmost element
      pivot = size(array1)

      ! Move all elements <= pivot to the LHS of the pivot
      ! Find position of the pivot in the array1 at the end (maxSmall)
      maxSmall = 0
      do i=1,size(array1)

        if( array1(i) <= array1(pivot)) then
          maxSmall = maxSmall + 1
          call swap(array1(i), array2(i), array1(maxSmall), array2(maxSmall))
        end if
      end do

      ! Recursivly sort the sub arrays divided by the pivot
      call quickSort(array1(1:maxSmall-1), array2(1:maxSmall-1))
      call quickSort(array1(maxSmall+1:size(array1)), array2(maxSmall+1:size(array1)))
    end if

  end subroutine quickSort_defReal_defReal

  !!
  !! Swap two integers
  !! Use XOR to avoid extra space (which is completly unnecessary)
  !!
  elemental subroutine swap_shortInt(i1,i2)
    integer(shortInt), intent(inout) :: i1
    integer(shortInt), intent(inout) :: i2
    integer(shortInt)                :: temp

    temp = i1
    i1 = i2
    i2 = temp

  end subroutine swap_shortInt

  !! Subroutine 'swap_shortIntArray'
  !!
  !! Basic description:
  !!   Swaps two entries within a shortInt array.
  !!
  !! Error:
  !!   None.
  !!
  pure subroutine swap_shortIntArray(array, idx1, idx2)
    integer(shortInt), dimension(:), intent(inout) :: array
    integer(shortInt), intent(in)                  :: idx1, idx2
    integer(shortInt)                              :: temp
    
    temp = array(idx1)
    array(idx1) = array(idx2)
    array(idx2) = temp
  
  end subroutine swap_shortIntArray

  !!
  !! Swap two reals
  !!
  elemental subroutine swap_defReal(r1,r2)
    real(defReal), intent(inout) :: r1
    real(defReal), intent(inout) :: r2
    real(defReal)                :: temp

    temp = r1
    r1 = r2
    r2 = temp

  end subroutine swap_defReal

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

  !!
  !! Swap character of length nameLen
  !!
  elemental subroutine swap_char_nameLen(c1,c2)
    character(nameLen), intent(inout) :: c1
    character(nameLen), intent(inout) :: c2
    character(nameLen)                :: temp

    temp = c1
    c1 = c2
    c2 = temp

  end subroutine swap_char_nameLen

  pure function charCshift(string,N) result(shifted)
    character(*), intent(in)      :: string
    integer(shortInt), intent(in) :: N
    character(len(string))        :: shifted
    integer(shortInt)             :: i, i_swap

    ! Loop over character and copy character with an offset
    do i=1,len(string)
      i_swap = modulo(i+N-1,len(string))+1
      shifted(i_swap:i_swap) = string(i:i)
    end do

  end function charCShift

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
  function printParticleType(type) result(str)
    integer(shortInt), intent(in) :: type
    character(:), allocatable     :: str

    select case(type)
      case(P_NEUTRON_CE)
        str = "CE Neutron"

      case(P_NEUTRON_MG)
        str = "MG Neutron"

      case default
        str = "Unknown "// numToChar(type)
    end select
  end function printParticleType


  !!
  !! Prints Scone ACII Header
  !!
  subroutine printStart()
    print *, repeat(" ><((((*> ",10)
    print *, ''
    print * ,"        _____ __________  _   ________  "
    print * ,"       / ___// ____/ __ \/ | / / ____/  "
    print * ,"       \__ \/ /   / / / /  |/ / __/     "
    print * ,"      ___/ / /___/ /_/ / /|  / /___     "
    print * ,"     /____/\____/\____/_/ |_/_____/     "
    print * , ''
    print * , ''
    print * , "Compiler Info :   ", compiler_version()
#ifdef _OPENMP
    print '(A, I4)', " OpenMP Threads: ", ompGetMaxThreads()
#endif
    print *
    print *, repeat(" <*((((>< ",10)

    ! TODO: Add extra info like date & time

  end subroutine printStart

  !!
  !! Prints line of fishes swiming right with an offset
  !!
  subroutine printFishLineR(offset)
    integer(shortInt),intent(in) :: offset
    integer(shortInt)            :: offset_L
    character(100), parameter    :: line = repeat(" ><((((*> ",10)
    character(100),dimension(10), parameter :: lines = [ &
    " ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*> " ,&
    "  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>" ,&
    ">  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*" ,&
    "*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((" ,&
    "(*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><(((" ,&
    "((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((" ,&
    "(((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><(" ,&
    "((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><" ,&
    "<((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  >" ,&
    "><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  " ]

    offset_L = modulo(offset,10)

    print *, lines(offset_L+1)

  end subroutine  printFishLineR

end module genericProcedures