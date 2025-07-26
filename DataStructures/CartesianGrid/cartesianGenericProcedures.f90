module cartesianGenericProcedures

  use vertexShelf_class,            only : vertexShelf
  use edgeShelf_class,              only : edgeShelf
  use faceShelf_class,              only : faceShelf
  use elementShelf_class,           only : elementShelf
  use genericProcedures,            only : findCommon
  use numPrecision   

  implicit none

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
  pure function binarySearchAngle(array, value) result(idx)
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

end module cartesianGenericProcedures