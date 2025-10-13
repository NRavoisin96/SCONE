module cartesianGenericProcedures

  use vertexShelf_class,            only : vertexShelf
  use edgeShelf_class,              only : edgeShelf
  use faceShelf_class,              only : faceShelf
  use elementShelf_class,           only : elementShelf
  use genericProcedures,            only : findCommon
  use numPrecision   

  implicit none

  !!
  !! DERIVED TYPE: num_freq
  !! A custom data structure to hold a number and its corresponding frequency.
  !! Used for: sortByHighestFrequency
  type :: num_freq
      integer :: number
      integer :: frequency
  end type num_freq

contains

  !!
  !!
  !!
  pure function findMinFaceAngle(edges, faces) result(maxCosValue)
    class(edgeShelf), intent(in)                        :: edges
    class(faceShelf), intent(in)                        :: faces
    integer(shortInt)                                   :: i, j, k, l ,m, sign1, sign2
    real(defReal)                                       :: currCosValue, maxCosValue
    integer(shortInt), dimension(:), allocatable        :: currFaceEdgeIdxs
    integer(shortInt), dimension(2)                     :: currEdge1VertexIdxs, currEdge2VertexIdxs

    maxCosValue = -1.0d0

    do i = 1, faces % getSize()
      currFaceEdgeIdxs = faces % getFaceEdgeIdxs(i)

      do j = 1, size(currFaceEdgeIdxs) - 1
        currEdge1VertexIdxs = edges % getEdgeVertexIdxs(currFaceEdgeIdxs(j))

        do k = j + 1, size(currFaceEdgeIdxs)
          currEdge2VertexIdxs = edges % getEdgeVertexIdxs(currFaceEdgeIdxs(k))

          do l = 1, 2

            do m = 1, 2
              if (currEdge1VertexIdxs(l) == currEdge2VertexIdxs(m)) then

                ! correct the direction of unit vector of each edge
                if (l == 1) then
                  sign1 = 1
                else 
                  sign1 = -1
                end if

                if (m == 1) then
                  sign2 = 1
                else 
                  sign2 = -1
                end if

                ! calculate cosine value
                currCosValue = dot_product(edges % getEdgeUnitvector(currFaceEdgeIdxs(j)), &
                                           edges % getEdgeUnitvector(currFaceEdgeIdxs(k)))*sign1*sign2

                ! update maxCosValue
                if (maxCosValue < currCosValue) maxCosValue = currCosValue

              end if

            end do

          end do

        end do

      end do

    end do

  end function findMinFaceAngle

  !!
  !!
  !!
  pure function findMinDihedralAngle(edges, faces, elements) result(maxCosValue)
    class(edgeShelf), intent(in)                        :: edges
    class(faceShelf), intent(in)                        :: faces
    class(elementShelf), intent(in)                     :: elements
    integer(shortInt)                                   :: i, j, k
    real(defReal)                                       :: currCosValue, maxCosValue
    integer(shortInt), dimension(:), allocatable        :: currElementFaceIdxs, currElementEdgeIdxs
    integer(shortInt), dimension(2)                     :: candidateFaceIdxs, candidateElementIdxs, signArray

    maxCosValue = -1.0d0

    do i = 1, elements % getSize()
      currElementFaceIdxs = abs(elements % getElementFaceIdxs(i))
      currElementEdgeIdxs = elements % getElementEdgeIdxs(i)

      do j = 1, size(currElementEdgeIdxs)
        candidateFaceIdxs = findCommon(currElementFaceIdxs, edges % getEdgeFaceIdxs(currElementEdgeIdxs(j)))

        do k = 1, 2
          candidateElementIdxs = faces % getFaceElementIdxs(candidateFaceIdxs(k))

          if (i < candidateElementIdxs(1) .OR. i < candidateElementIdxs(2)) then
            signArray(k) = 1
          else
            signArray(k) = -1
          end if

        end do

        currCosValue = dot_product(faces % getFaceNormal(candidateFaceIdxs(1)), &
                                   faces % getFaceNormal(candidateFaceIdxs(2)))*signArray(1)*signArray(2)*(-1)

        ! update maxCosValue
        if (maxCosValue < currCosValue) maxCosValue = currCosValue
        
      end do

    end do

  end function findMinDihedralAngle

  !!
  !!
  !!


  !! (needs to be changed) (can be accelerated: use subroutine to calculate baseIntegerCoord of each vertex and)
  !! (write vertices subroutine to set it and vertices function to get it. Then, use ishft and iand for each layer)
  !! (, and we do not have to calculate indices fresh for each layer)
  !!
  !! AABBIndices = [xmin, ymin, zmin, xmax, ymax, zmax]
  pure function constructAABB(vertices, currVertexIdxs, gridBounds_min, gridSpacing) result(AABBIndices)
    class(vertexShelf), intent(in)                         :: vertices
    integer(shortInt), dimension(:), intent(in)            :: currVertexIdxs
    real(defReal), dimension(3), intent(in)                :: gridBounds_min
    real(defReal), intent(in)                              :: gridSpacing
    integer(shortInt), dimension(6)                        :: AABBIndices
    real(defreal), dimension(3)                            :: xyz_min, xyz_max, currVertexCoords
    integer(shortInt)                                      :: i, j

    ! initialise xyz_min and xyz_max using the first vertex
    currVertexCoords = vertices % getVertexCoordinates(currVertexIdxs(1))
    xyz_max = currVertexCoords
    xyz_min = currVertexCoords

    ! find xyz_min and xyz_max 
    do i = 2, size(currVertexIdxs)
        currVertexCoords = vertices % getVertexCoordinates(currVertexIdxs(i))

        do j = 1, 3
            if (xyz_min(j) > currVertexCoords(j)) then
                 xyz_min(j) = currVertexCoords(j)
            elseif (xyz_max(j) < currVertexCoords(j)) then
                xyz_max(j) = currVertexCoords(j)
            end if
        end do

    end do

    ! find AABBIndices
    do i = 1, 3
        AABBIndices(i) = ceiling((xyz_min(i) - gridBounds_min(i))/(gridSpacing)) 
        AABBIndices(3+i) = ceiling((xyz_max(i) - gridBounds_min(i))/(gridSpacing)) 
    end do

  end function constructAABB

  !! insertion sort O(N^2). There exists cheaper sorting algorithm (quickSort O(N logN)) 
  !! but for N < 8 (which is mostly the case in FEM), insertion sort is better because quicksort
  !! needs extra procedures such as selecting pivots, ....
  !!
  !! sorts arraysReal so that its values are increasing with index. Array Int are 
  !! sorted using the exactly the same swaps made during arrayReal sorting process.
  ! (needs to be changed) (move to genericProcedures)
  pure subroutine sortPairs(arrayReal, arrayInt)
    real(defReal), dimension(:), intent(inout)            :: arrayReal
    integer(shortInt), dimension(:), intent(inout)        :: arrayInt
    integer(shortInt)                                     :: i, j
    real(defReal)                                         :: key_real
    integer(shortInt)                                     :: key_Int

    do i = 2, size(arrayReal)
       key_real = arrayReal(i);  key_Int = arrayInt(i)
       j = i - 1
       do while (j >= 1 .and. arrayReal(j) > key_real)
          arrayReal(j+1) = arrayReal(j)
          arrayInt(j+1) = arrayInt(j)
          j      = j - 1
       end do
       arrayReal(j+1) = key_real
       arrayInt(j+1) = key_Int
    end do

  end subroutine sortPairs

  !!
  !!
  !! (needs to be changed) (the binarySearch subroutine in generic procedure does not allow)
  !! ("value" to be outside the bounds of "array". So rewritten. Can be moved to genericProcedure Later)
  !! (needs to be changed) (for boundary edges, value > array(size(array)) or value < array(1) cases might not)
  !! (work using this subroutine)
  function binarySearchAngle(array, value) result(idx)
    real(defReal), dimension(:), intent(in)             :: array
    real(defReal), intent(in)                           :: value
    integer(shortInt)                                   :: idx, bottom, top, i

    ! in case of value being outside the ranges of "array", manually assign idx = size(array).
    if (value > array(size(array))) then
      idx = size(array)
      return
    elseif (value < array(1)) then
      idx = size(array)
      return
    end if

    ! Find Top and Bottom Index Array
    bottom = 1
    top = size(array)

    do i = 1,100
      !Calculate mid point
      idx = (top + bottom)*0.5

      ! Termination condition
      if (bottom == idx) return

      ! Binary Step
      if (array(idx) <= value) then
        bottom = idx
      else
        top = idx
      end if
    end do

  end function binarySearchAngle

  !!
  !!
  !!
  ! (needs to be changed) (move to genericProcedures.f90)
  elemental function testIntervalIntersection(min1, max1, min2, max2) result(intersect)
    real(defReal), intent(in)                           :: min1, max1, min2, max2
    logical                                             :: intersect

    intersect = ((min1 <= max2) .AND. (min2 <= max1))

  end function testIntervalIntersection

  !!
  !!
  !!
  ! (needs to be changed) (written for temp operation; can be optimised)
  elemental function calculateAvgEdgeLength(edges) result(avgEdgeLength)
    class(edgeShelf), intent(in)                        :: edges
    integer(shortInt)                                   :: i
    real(defReal)                                       :: sumEdgeLength, avgEdgeLength

    sumEdgeLength = 0
    do i = 1, edges % getSize()
      sumEdgeLength = sumEdgeLength + edges % getEdgeLength(i)
    end do

    avgEdgeLength = sumEdgeLength / (edges % getSize())

  end function calculateAvgEdgeLength

  !!
  !!
  !!
  pure function getGlobalIdx(baseIntegerCoord, shift) result(globalIdx)
    integer(shortInt), dimension(3), intent(in)     :: baseIntegerCoord, shift
    integer(shortInt), dimension(3)                 :: globalIdx
    integer(shortInt)                               :: i

    globalIdx = ishft(baseIntegerCoord, -shift) + 1
    
  end function getGlobalIdx

  !!
  !!
  !!
  pure function getLocalIdx(baseIntegerCoord, shift, mask) result(localIdx)
    integer(shortInt), dimension(3), intent(in)     :: baseIntegerCoord, shift, mask
    integer(shortInt), dimension(3)                 :: localIdx
    integer(shortInt)                               :: i

    localIdx = iand( ishft(baseIntegerCoord, -shift),  mask ) + 1
    
  end function getLocalIdx

  !!
  !!
  !! (needs to be changed) (remove shift)
  pure function getLocalIdxFinest(baseIntegerCoord, shift, mask) result(localIdx)
    integer(shortInt), dimension(3), intent(in)     :: baseIntegerCoord, shift, mask
    integer(shortInt), dimension(3)                 :: localIdx
    integer(shortInt)                               :: i

    localIdx = iand( baseIntegerCoord,  mask ) + 1
    
  end function getLocalIdxFinest

  !!
  !!
  !! (needs to be changed) (can go to genericProcedures)
  function getUniqueSortedArr(arr) result(res)
    integer(shortInt), intent(in)   :: arr(:)
    integer(shortInt), allocatable  :: res(:)
    integer(shortInt), allocatable  :: temp(:)
    integer(shortInt)               :: n, i, m

    n = size(arr)
    if (n == 0) then
        allocate(res(0))
        return
    end if

    allocate(temp(n))
    temp = arr

    ! sort temp in-place
    call quicksort(temp, 1, n)

    ! compact unique values
    m = 0
    do i = 1, n
        if (i == 1 .or. temp(i) /= temp(i-1)) then
          m = m + 1
          temp(m) = temp(i)
        end if
    end do

    allocate(res(m))
    res = temp(1:m)

  end function getUniqueSortedArr

  !!
  !!
  !! (needs to be changed) (there is the same in genericProcedures already?)
  recursive subroutine quicksort(a, left, right)
    integer(shortInt), intent(inout) :: a(:)
    integer(shortInt), intent(in)    :: left, right
    integer(shortInt)                :: i, j, pivot, tmp

    if (left >= right) return

    pivot = a((left + right) / 2)
    i = left
    j = right

    do
       do while (a(i) < pivot)
          i = i + 1
       end do
       do while (a(j) > pivot)
          j = j - 1
       end do

       if (i <= j) then
          tmp = a(i)
          a(i) = a(j)
          a(j) = tmp
          i = i + 1
          j = j - 1
       end if

       if (i > j) exit
    end do

    if (left < j) call quicksort(a, left, j)
    if (i < right) call quicksort(a, i, right)
  end subroutine quicksort

  !!
  !!
  !!
  function intBinarySearch(arr, value) result(idx)
      integer(shortInt), dimension(:), intent(in) :: arr
      integer(shortInt), intent(in)               :: value
      integer(shortInt)                           :: idx
      integer(shortInt)                           :: left, right, mid

      left  = 1
      right = size(arr)
      idx   = 0

      do while (left <= right)
          mid = (left + right) / 2
          if (arr(mid) == value) then
              idx = mid
              exit
          else if (arr(mid) < value) then
              left = mid + 1
          else
              right = mid - 1
          end if
      end do

  end function intBinarySearch

!!
!!
!! Takes an integer array, finds all unique numbers, and returns a new
!! array with those unique numbers sorted in descending order of their
!! frequency in the original array.
function sortByHighestFrequency(input_array) result(output_array)
    integer(shortInt), dimension(:), intent(in)  :: input_array
    integer(shortInt), dimension(:), allocatable :: output_array
    integer(shortInt)                            :: array_size, num_unique, i, j
    integer(shortInt), dimension(:), allocatable :: sorted_temp_array
    type(num_freq), dimension(:), allocatable    :: freq_pairs

    array_size = size(input_array)

    ! Handle an empty input array
    if (array_size == 0) then
        if (allocated(output_array)) deallocate(output_array)
        allocate(output_array(0))
        return
    end if

    ! --- Step 1: Create a sorted copy of the input array to make counting easy ---
    sorted_temp_array = input_array
    call quicksort(sorted_temp_array, 1, array_size)

    ! --- Step 2: First pass to count the number of unique elements ---
    num_unique = 1
    do i = 2, array_size
        if (sorted_temp_array(i) /= sorted_temp_array(i-1)) then
            num_unique = num_unique + 1
        end if
    end do

    ! --- Step 3: Second pass to populate the frequency pairs array ---
    allocate(freq_pairs(num_unique))
    j = 1
    freq_pairs(j)%number = sorted_temp_array(1)
    freq_pairs(j)%frequency = 1
    do i = 2, array_size
        if (sorted_temp_array(i) == sorted_temp_array(i-1)) then
            freq_pairs(j)%frequency = freq_pairs(j)%frequency + 1
        else
            j = j + 1
            freq_pairs(j)%number = sorted_temp_array(i)
            freq_pairs(j)%frequency = 1
        end if
    end do

    deallocate(sorted_temp_array)

    ! --- Step 4: Sort the frequency pairs in descending order of frequency ---
    call quicksort_freq_pairs(freq_pairs, 1, num_unique)

    ! --- Step 5: Construct the final output array from the sorted pairs ---
    allocate(output_array(num_unique))
    do i = 1, num_unique
        output_array(i) = freq_pairs(i)%number
    end do

    deallocate(freq_pairs)

end function sortByHighestFrequency

!!
!! An efficient, recursive sorting algorithm for the 'num_freq' type. It
!! sorts the array in descending order based on the 'frequency' field.
recursive subroutine quicksort_freq_pairs(pairs, left, right)
    type(num_freq), dimension(:), intent(inout) :: pairs
    integer, intent(in)                         :: left, right
    integer                                     :: i, j, pivot_freq
    type(num_freq)                              :: temp

    if (left >= right) return

    ! Use the frequency of the middle element as the pivot
    pivot_freq = pairs((left + right) / 2)%frequency
    i = left
    j = right

    do
        ! Sort descending: find elements > pivot on the left
        do while (pairs(i)%frequency > pivot_freq)
            i = i + 1
        end do
        ! Sort descending: find elements < pivot on the right
        do while (pairs(j)%frequency < pivot_freq)
            j = j - 1
        end do

        if (i <= j) then
            temp = pairs(i)
            pairs(i) = pairs(j)
            pairs(j) = temp
            i = i + 1
            j = j - 1
        end if

        if (i > j) exit
    end do

    if (left < j) call quicksort_freq_pairs(pairs, left, j)
    if (i < right) call quicksort_freq_pairs(pairs, i, right)
end subroutine quicksort_freq_pairs

!!
!!
!!
subroutine remove_elements(a, idxs)
    implicit none
    integer(shortInt), allocatable, intent(inout) :: a(:)
    integer(shortInt),            intent(in)     :: idxs(:)

    integer(shortInt), allocatable :: iv(:)   ! valid, unique, sorted indices to drop
    integer(shortInt), allocatable :: tmp(:)
    integer :: n, m, r, i, j, k, out, t

    n = size(a)
    m = size(idxs)
    if (n == 0 .or. m == 0) return

    ! Collect in-range indices, remove duplicates
    allocate(iv(m))
    r = 0
    do i = 1, m
        k = idxs(i)
        if (k >= 1 .and. k <= n) then
            if (r == 0) then
                r = 1
                iv(1) = k
            else
                if (.not. any(iv(1:r) == k)) then
                    r = r + 1
                    iv(r) = k
                end if
            end if
        end if
    end do
    if (r == 0) then
        deallocate(iv)
        return
    end if

    ! Sort iv(1:r); r is tiny so a simple O(r^2) sort is optimal in practice
    if (r > 1) then
        do i = 1, r - 1
            do j = i + 1, r
                if (iv(j) < iv(i)) then
                    t = iv(i); iv(i) = iv(j); iv(j) = t
                end if
            end do
        end do
    end if

    ! Compact copy in one pass
    allocate(tmp(n - r))
    out = 0
    j = 1
    do k = 1, n
        if (j <= r .and. k == iv(j)) then
            j = j + 1
        else
            out = out + 1
            tmp(out) = a(k)
        end if
    end do

    call move_alloc(tmp, a)
    deallocate(iv)
end subroutine remove_elements


end module cartesianGenericProcedures