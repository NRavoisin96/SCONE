module cartesianInitProcedures

  use universalVariables,           only : ZERO, INF
  use vertexShelf_class,            only : vertexShelf
  use edgeShelf_class,              only : edgeShelf
  use faceShelf_class,              only : faceShelf
  use elementShelf_class,           only : elementShelf
  use genericProcedures,            only : findCommon, append, eraseAt, fatalError
  use numPrecision   
  use cartesianGenericProcedures,   only : testIntervalIntersection

  implicit none

contains

  !!
  !!
  !!
  pure subroutine testEdgeIntersection(vertices, edges, edgeIdx, circumscribedBallRadius, targetDistance, centroid, &
                                       currEdgeVector, currVertexIdxs, a, cellToVertexIdx, cellToEdgeIdx)
    class(vertexShelf), intent(in)              :: vertices
    class(edgeShelf), intent(in)                :: edges
    integer(shortInt), intent(in)               :: edgeIdx
    integer(shortInt), intent(inout)            :: cellToVertexIdx, cellToEdgeIdx
    real(defReal), intent(in)                   :: circumscribedBallRadius, targetDistance, a
    real(defReal), dimension(3), intent(in)     :: centroid, currEdgeVector
    integer(shortInt), dimension(2), intent(in) :: currVertexIdxs
    real(defReal), dimension(3)                 :: dummyVector
    real(defReal)                               :: b, c, discriminant, edgeLength, inverseA, sqrtDiscriminant
    real(defReal), dimension(2)                 :: t, targetDistanceRatio
    integer(shortInt)                           :: i

    ! calculate constants
    edgeLength = edges % getEdgeLength(edgeIdx)

    ! calculate dummyVector = vertex1 coordinate - centroid of the cube
    dummyVector = vertices % getVertexCoordinates(currVertexIdxs(1)) - centroid

    ! calculate constants for quadratic formula
    b = dot_product(dummyVector, currEdgeVector)
    c = dot_product(dummyVector, dummyVector) - circumscribedBallRadius * circumscribedBallRadius
    discriminant = b * b - a * c
  
    ! calculate t and leave the subroutine if t is invalid
    if (discriminant < ZERO .or. (discriminant == ZERO .and. ZERO < b)) return
    sqrtDiscriminant = sqrt(discriminant)
    inverseA = ONE / a
    t = [-b - sqrtDiscriminant, -b + sqrtDiscriminant] * inverseA
    if (ONE < t(1) .or. t(2) < ZERO) return

    ! If passed to this point, t is valid. Hence, construct phi mapping
    cellToVertexIdx = merge(currVertexIdxs(1), currVertexIdxs(2), t(2) < HALF)

    ! test and constuct phiCapital mapping
    targetDistanceRatio(1) = targetDistance / edgeLength
    targetDistanceRatio(2) = ONE - targetDistanceRatio(1)
    do i = 1, 2
      if (targetDistanceRatio(1) < t(i) .and. t(i) < targetDistanceRatio(2)) then
        cellToEdgeIdx = edgeIdx
        return

      end if

    end do

  end subroutine testEdgeIntersection

  !!
  !!
  !!
  subroutine testPolyhedronInclusion(faces, currElementFaceIdxs, centroid, &
                                     faceNormalSigns, elementIdx, chi)
    class(faceShelf), intent(in)                        :: faces
    integer(shortInt), dimension(:), intent(in)         :: currElementFaceIdxs
    real(defReal), dimension(3), intent(in)             :: centroid
    real(defReal), dimension(:,:), intent(in)           :: faceNormalSigns
    integer(shortInt), intent(in)                       :: elementIdx
    integer(shortInt), intent(inout)                    :: chi
    integer(shortInt)                                   :: i
    real(defReal), dimension(3)                         :: furthestVertexCoord

    do i = 1, size(currElementFaceIdxs)
      furthestVertexCoord = centroid + faceNormalSigns(:, i) 
      if (dot_product(faces % getFaceNormal(currElementFaceIdxs(i)), furthestVertexCoord) + &
          faces % getFaceConst(currElementFaceIdxs(i)) > ZERO) return

    end do

    ! if survived to this point, then the cell is entired enclosed by the polyhedron. 
    ! Hence, set chi mapping to the current element index.
    chi = elementIdx

  end subroutine testPolyhedronInclusion

  !!
  !!
  !!
  subroutine testPolyhedronInclusion2Coarsest(faces, currElementFaceIdxs, centroid, &
                                              elementIdx, chi)
    class(faceShelf), intent(in)                        :: faces
    integer(shortInt), dimension(:), intent(in)         :: currElementFaceIdxs
    real(defReal), dimension(3), intent(in)             :: centroid
    integer(shortInt), intent(in)                       :: elementIdx
    integer(shortInt), intent(inout)                    :: chi
    integer(shortInt)                                   :: i

    do i = 1, size(currElementFaceIdxs)

      if (dot_product(faces % getFaceNormal(currElementFaceIdxs(i)), centroid) &
                + faces % getFaceExtraDistanceArr(currElementFaceIdxs(i), 1) &
                + faces % getFaceConst(currElementFaceIdxs(i)) > ZERO) then
          ! if TRUE, then at least a part of this cell lies outside of the polyhedron

          return
      end if 

    end do

    ! if survived to this point, then the cell is entired enclosed by the polyhedron. 
    ! Hence, set chi mapping to the current element index.
    chi = elementIdx

  end subroutine testPolyhedronInclusion2Coarsest

  !!
  !!
  !!
  subroutine testPolyhedronInclusion2(faces, currElementFaceIdxs, centroid, &
                                     faceNormalSigns, elementIdx, chi, &
                                     removedFaceIdxsInArr, isOut)
    class(faceShelf), intent(in)                              :: faces
    integer(shortInt), dimension(:), intent(in)               :: currElementFaceIdxs
    real(defReal), dimension(3), intent(in)                   :: centroid
    real(defReal), dimension(:,:), intent(in)                 :: faceNormalSigns
    integer(shortInt), intent(in)                             :: elementIdx
    integer(shortInt), intent(inout)                          :: chi
    integer(shortInt), dimension(:), allocatable, intent(out) :: removedFaceIdxsInArr
    logical(defBool), intent(out)                             :: isOut
    integer(shortInt)                                         :: i, j
    real(defReal), dimension(3)                               :: furthestVertexCoord
    logical(defBool)                                          :: isInside

    isOut = .FALSE.
    isInside = .TRUE.
    ! Add new test to set isInside False

    do i = 1, size(currElementFaceIdxs)
      
      do j = 1, 3
        furthestVertexCoord(j) = centroid(j) + faceNormalSigns(j,i) 
      end do

      if (dot_product(faces % getFaceNormal(currElementFaceIdxs(i)), furthestVertexCoord) &
          + faces % getFaceConst(currElementFaceIdxs(i)) > ZERO) then

        isInside = .FALSE.

      else

        call append(removedFaceIdxsInArr, i)

      end if 

    end do

    ! if survived to this point, then the cell is entired enclosed by the polyhedron. 
    ! Hence, set chi mapping to the current element index.
    if (isInside) chi = elementIdx



  end subroutine testPolyhedronInclusion2

  !!
  !!
  !!
  subroutine testPolyhedronInclusion2New(faces, currElementFaceIdxs, centroid, &
                                     elementIdx, chi, &
                                     removedFaceIdxsInArr, isOut, currLayer)
    class(faceShelf), intent(in)                              :: faces
    integer(shortInt), dimension(:), intent(in)               :: currElementFaceIdxs
    real(defReal), dimension(3), intent(in)                   :: centroid
    integer(shortInt), intent(in)                             :: elementIdx, currLayer
    integer(shortInt), intent(inout)                          :: chi
    integer(shortInt), dimension(:), allocatable, intent(out) :: removedFaceIdxsInArr
    logical(defBool), intent(out)                             :: isOut
    integer(shortInt)                                         :: i
    logical(defBool)                                          :: isInside
    ! real(defReal)                                             :: vectorDotCentroid, extraDistance, const
    !integer(shortInt), dimension(:), allocatable                  :: abc

    isOut = .FALSE.
    isInside = .TRUE.
    !print*, size(currElementFaceIdxs)
    ! Add new test to set isInside False

    do i = 1, size(currElementFaceIdxs)
      if (dot_product(faces % getFaceNormal(currElementFaceIdxs(i)), centroid) &
          + faces % getFaceExtraDistanceArr(currElementFaceIdxs(i), currLayer) &
          + faces % getFaceConst(currElementFaceIdxs(i)) > ZERO) then

        isInside = .FALSE.                                          

      else

        call append(removedFaceIdxsInArr, i)

      end if 

    end do

    ! if (size(currelementFaceIdxs) == size(removedFaceIdxsInArr)) then
    !   if (.NOT. isInside) then
    !     print*, "---------------------------------------"
    !     print*, currElementFaceIdxs
    !     print*, abc
    !   end if
    ! end if

    ! if survived to this point, then the cell is entired enclosed by the polyhedron. 
    ! Hence, set chi mapping to the current element index.
    if (isInside) chi = elementIdx



  end subroutine testPolyhedronInclusion2New

  !!
  !!
  !!
  function testPolyhedronInclusionNew(faces, currElementFaceIdxs, centroid, &
                                     faceNormalSigns) result(isIncluded)
    class(faceShelf), intent(in)                        :: faces
    integer(shortInt), dimension(:), intent(in)         :: currElementFaceIdxs
    real(defReal), dimension(3), intent(in)             :: centroid
    real(defReal), dimension(:,:), intent(in)           :: faceNormalSigns
    integer(shortInt)                                   :: i, j
    real(defReal), dimension(3)                         :: furthestVertexCoord
    logical(defBool)                                    :: isIncluded

    isIncluded = .FALSE.

    do i = 1, size(currElementFaceIdxs)
      
      do j = 1, 3
        furthestVertexCoord(j) = centroid(j) + faceNormalSigns(j,i) 
      end do

      if (dot_product(faces % getFaceNormal(currElementFaceIdxs(i)), furthestVertexCoord) &
          + faces % getFaceConst(currElementFaceIdxs(i)) > 0) then

          ! if TRUE, then at least a part of this cell lies outside of the polyhedron
          return

      end if 

    end do

    ! if survived to this point, then the cell is entired enclosed by the polyhedron. 
    isIncluded = .TRUE.

  end function testPolyhedronInclusionNew

  !!
  !!
  !!
  subroutine testFaceIntersection(vertices, edges, faces, currVertexIdxs, extraDistance, currFaceNormal, centroid, &
                                  cellSpacing, faceIdx, currFaceEdgeIdxs, targetDistance, chi, phi, phiCapital, faceIdxs)
    class(vertexShelf), intent(in)                 :: vertices
    class(edgeShelf), intent(in)                   :: edges
    class(faceShelf), intent(in)                   :: faces
    integer(shortInt), dimension(:), intent(in)    :: currVertexIdxs, currFaceEdgeIdxs
    real(defReal), dimension(3), intent(in)        :: currFaceNormal, centroid
    real(defReal), intent(in)                      :: extraDistance, cellSpacing, targetDistance
    integer(shortInt), intent(in)                  :: faceIdx, chi
    integer(shortInt), intent(inout)               :: phi, phiCapital
    integer(shortInt), dimension(2), intent(inout) :: faceIdxs
    integer(shortInt)                              :: i, j
    real(defReal), dimension(3)                    :: currVertexCoords, min1, max1, currEdgeUnitVector
    real(defReal)                                  :: faceConst, currValue, vectorDotCentroid, extraDistance2 

    !------------------------------------------------------------------------------------------------
    !if cell is contained within a polyhedron, the cell cannot intersect with the polyhedron's faces
    !Or, if phi and phiCapital mapping informations are assigned already from edge interesetion, no need to find new one.
    !------------------------------------------------------------------------------------------------
    if (chi /= 0) return
    if (phi /= 0 .and. phiCapital /= 0) return

    !------------------------------------------------------------------------------------------------
    ! Testing along face normal
    !------------------------------------------------------------------------------------------------
    vectorDotCentroid = dot_product(currFaceNormal, centroid)
    faceConst = faces % getFaceConst(faceIdx)
    if (.not. testIntervalIntersection(vectorDotCentroid - extraDistance, vectorDotCentroid + extraDistance, &
        -faceConst, -faceConst)) return

    !------------------------------------------------------------------------------------------------
    ! Testing along crossProduct(each of edgeUnitVector, three coordinate basis)
    !------------------------------------------------------------------------------------------------
    ! loop over all edgesUnitVectors of the current face
    do i = 1, size(currFaceEdgeIdxs)
      currEdgeUnitVector = edges % getEdgeUnitvector(currFaceEdgeIdxs(i))
      
      ! loop over all vertices of the current face to find the min and max of polygon interval
      do j = 1, size(currVertexIdxs)
        currVertexCoords = vertices % getVertexCoordinates(currVertexIdxs(j))

        ! for coordinate basis = (1,0,0)
        currValue = currEdgeUnitVector(3) * currVertexCoords(2) - currEdgeUnitVector(2) * currVertexCoords(3)
        if (j == 1) then
          min1(1) = currValue
          max1(1) = currValue

        else
          min1(1) = min(min1(1), currValue)
          max1(1) = max(max1(1), currValue)

        end if

        ! for coordinate basis = (0,1,0)
        currValue = -currEdgeUnitVector(3) * currVertexCoords(1) + currEdgeUnitVector(1) * currVertexCoords(3)
        if (j == 1) then
          min1(2) = currValue
          max1(2) = currValue

        else
          min1(2) = min(min1(2), currValue)
          max1(2) = max(max1(2), currValue)

        end if

        ! for coordinate basis = (0,0,1)
        currValue = currEdgeUnitVector(2) * currVertexCoords(1) - currEdgeUnitVector(1) * currVertexCoords(2)
        if (j == 1) then
          min1(3) = currValue
          max1(3) = currValue

        else
          min1(3) = min(min1(3), currValue)
          max1(3) = max(max1(3), currValue)

        end if

      end do

      ! test intersections of interval
      vectorDotCentroid = currEdgeUnitVector(3) * centroid(2) - currEdgeUnitVector(2) * centroid(3)
      extraDistance2 = HALF * sum(abs(currEdgeUnitVector([2, 3]))) * cellSpacing
      if (.not. testIntervalIntersection(vectorDotCentroid - extraDistance2, vectorDotCentroid + extraDistance2, &
                                         min1(1), max1(1))) return

      vectorDotCentroid = currEdgeUnitVector(1) * centroid(3) - currEdgeUnitVector(3) * centroid(1)
      extraDistance2 = HALF * sum(abs(currEdgeUnitVector([1, 3]))) * cellSpacing
      if (.not. testIntervalIntersection(vectorDotCentroid - extraDistance2, vectorDotCentroid + extraDistance2, &
                                         min1(2), max1(2))) return

      vectorDotCentroid = currEdgeUnitVector(2) * centroid(1) - currEdgeUnitVector(1) * centroid(2)
      extraDistance2 = HALF * sum(abs(currEdgeUnitVector([1, 2]))) * cellSpacing
      if (.not. testIntervalIntersection(vectorDotCentroid - extraDistance2, vectorDotCentroid + extraDistance2, &
                                         min1(3), max1(3))) return

    end do

    !------------------------------------------------------------------------------------------------
    ! if survived to this point, then there is no separating axis. Hence, construct mapping accordingly
    !------------------------------------------------------------------------------------------------

    ! (needs to be changed) (due to memory, taking shortCut)
    ! (originally, we just have to add faceIdx to an array of intersected faces)
    if (faceIdxs(1) == 0) then
      faceIdxs(1) = faceIdx

    else
      faceIdxs(2) = faceIdx
      ! tests if the two faces have a common edge. If yes, assign phi and phiCapital mappings.
      ! if there are more than two faces intersected with the cell, and phi and phicaptial have been assigned already,
      ! then, this cell is complete in terms of mapping construction, and this subroutine is terminated at the beginning
      ! of this subroutine.
      call testTwoIntersectedFaces(vertices, edges, faces, centroid, targetDistance, phi, phiCapital, faceIdxs)

    end if

  end subroutine testFaceIntersection

  !!
  !!
  !!
  subroutine testFaceIntersection2(vertices, edges, faces, currVertexIdxs, currFaceNormal, centroid, &
                                  cellSpacing, faceIdx, currFaceEdgeIdxs, targetDistance, chi, phi, phiCapital, faceIdxs, &
                                  currLayer)
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(in)                        :: edges
    class(faceShelf), intent(in)                        :: faces
    integer(shortInt), dimension(:), intent(in)         :: currVertexIdxs, currFaceEdgeIdxs
    real(defReal), dimension(3), intent(in)             :: currFaceNormal, centroid
    real(defReal), intent(in)                           :: cellSpacing, targetDistance
    integer(shortInt), intent(in)                       :: faceIdx, chi, currLayer
    integer(shortInt), intent(inout)                    :: phi, phiCapital
    integer(shortInt), dimension(2), intent(inout)      :: faceIdxs
    integer(shortInt)                                   :: i, j
    real(defReal), dimension(3)                         :: currVertexCoords, min1, max1, currEdgeUnitVector
    real(defReal)                                       :: faceConst, currValue, vectorDotCentroid, extraDistance2, &
                                                          extraDistance 
    extraDistance = faces % getFaceExtraDistanceArr(faceIdx, currLayer)  
    !------------------------------------------------------------------------------------------------
    !if cell is contained within a polyhedron, the cell cannot intersect with the polyhedron's faces
    !Or, if phi and phiCapital mapping informations are assigned already from edge interesetion, no need to find new one.
    !------------------------------------------------------------------------------------------------
    if (chi /= 0) then
      return
    elseif (phi /= 0 .AND. phiCapital /= 0) then
      return
    end if

    !------------------------------------------------------------------------------------------------
    ! Testing along face normal
    !------------------------------------------------------------------------------------------------
    vectorDotCentroid = dot_product(currFaceNormal, centroid)
    faceConst = faces % getFaceConst(faceIdx)
    ! extradistance2 = faces % getFaceExtraDistanceArr(faceIdx,4)
    ! if (extraDistance /= extraDistance2) then
    !   print*, "hahaha"
    ! end if

    if (.NOT. testIntervalIntersection(vectorDotCentroid - extraDistance, vectorDotCentroid + extraDistance, &
        -faceConst, -faceConst)) then
          !!!!!
          ! print*, "faceNormal"
          ! print*, vectorDotCentroid - extraDistance
          ! print*, vectorDotCentroid + extraDistance
          ! print*, -faceConst
          !!!!!
      return
    end if

    !------------------------------------------------------------------------------------------------
    ! Testing along crossProduct(each of edgeUnitVector, three coordinate basis)
    !------------------------------------------------------------------------------------------------
    ! loop over all edgesUnitVectors of the current face
    do i = 1, size(currFaceEdgeIdxs)
      currEdgeUnitVector = edges % getEdgeUnitvector(currFaceEdgeIdxs(i))
      
      ! loop over all vertices of the current face to find the min and max of polygon interval
      do j = 1, size(currVertexIdxs)
        currVertexCoords = vertices % getVertexCoordinates(currVertexIdxs(j))

        ! for coordinate basis = (1,0,0)
        currValue = currEdgeUnitVector(3)*currVertexCoords(2) - currEdgeUnitVector(2)*currVertexCoords(3)
        if (j == 1) then
          min1(1) = currValue
          max1(1) = currValue
        else
          if (currValue < min1(1)) then
            min1(1) = currValue
          elseif (currValue > max1(1)) then
            max1(1) = currValue
          end if 
        end if

        ! for coordinate basis = (0,1,0)
        currValue = - currEdgeUnitVector(3)*currVertexCoords(1) + currEdgeUnitVector(1)*currVertexCoords(3)
        if (j == 1) then
          min1(2) = currValue
          max1(2) = currValue
        else
          if (currValue < min1(2)) then
            min1(2) = currValue
          elseif (currValue > max1(2)) then
            max1(2) = currValue
          end if 
        end if

        ! for coordinate basis = (0,0,1)
        currValue =  currEdgeUnitVector(2)*currVertexCoords(1) - currEdgeUnitVector(1)*currVertexCoords(2)
        if (j == 1) then
          min1(3) = currValue
          max1(3) = currValue
        else
          if (currValue < min1(3)) then
            min1(3) = currValue
          elseif (currValue > max1(3)) then
            max1(3) = currValue
          end if 
        end if

      end do

      ! test intersections of interval
      vectorDotCentroid = currEdgeUnitVector(3)*centroid(2)-currEdgeUnitVector(2)*centroid(3)
      extraDistance2 = (abs(currEdgeUnitVector(3)) + abs(currEdgeUnitVector(2)))*cellSpacing*0.5
      if (.NOT. testIntervalIntersection(vectorDotCentroid - extraDistance2, vectorDotCentroid + extraDistance2, &
                                                min1(1), max1(1))) then
          !!!!!
          ! print*, "Crossx"
          !!!!!
          return
      end if

      vectorDotCentroid = currEdgeUnitVector(1)*centroid(3)-currEdgeUnitVector(3)*centroid(1)
      extraDistance2 = (abs(currEdgeUnitVector(3)) + abs(currEdgeUnitVector(1)))*cellSpacing*0.5
      if (.NOT. testIntervalIntersection(vectorDotCentroid - extraDistance2, vectorDotCentroid + extraDistance2, &
                                                min1(2), max1(2))) then
          !!!!!
          ! print*, "Crossy"
          !!!!!
          return
      end if

      vectorDotCentroid = currEdgeUnitVector(2)*centroid(1)-currEdgeUnitVector(1)*centroid(2)
      extraDistance2 = (abs(currEdgeUnitVector(2)) + abs(currEdgeUnitVector(1)))*cellSpacing*0.5
      if (.NOT. testIntervalIntersection(vectorDotCentroid - extraDistance2, vectorDotCentroid + extraDistance2, &
                                                min1(3), max1(3))) then
          !!!!!
          ! print*, "Crossz"
          !!!!!
          return
      end if

    end do

    !------------------------------------------------------------------------------------------------
    ! if survived to this point, then there is no separating axis. Hence, construct mapping accordingly
    !------------------------------------------------------------------------------------------------

    ! (needs to be changed) (due to memory, taking shortCut)
    ! (originally, we just have to add faceIdx to an array of intersected faces)
    if (faceIdxs(1) == 0) then
      faceIdxs(1) = faceIdx
    else
      faceIdxs(2) = faceIdx
      ! tests if the two faces have a common edge. If yes, assign phi and phiCapital mappings.
      ! if there are more than two faces intersected with the cell, and phi and phicaptial have been assigned already,
      ! then, this cell is complete in terms of mapping construction, and this subroutine is terminated at the beginning
      ! of this subroutine.
      call testTwoIntersectedFaces(vertices, edges, faces, centroid, targetDistance, phi, phiCapital, faceIdxs)
    end if
    
    !!!!!
    ! print*, "PASSESD"
    !!!!!

  end subroutine testFaceIntersection2

  !!
  !!
  !!
  subroutine testFaceIntersectionCoarsest(vertices, edges, faces, currVertexIdxs, extraDistance, currFaceNormal, & 
                                           centroid, cellSpacing, faceIdx, currFaceEdgeIdxs, intersectedFaceIdxs)
    class(vertexShelf), intent(in)                              :: vertices
    class(edgeShelf), intent(in)                                :: edges
    class(faceShelf), intent(in)                                :: faces
    integer(shortInt), dimension(:), intent(in)                 :: currVertexIdxs, currFaceEdgeIdxs
    real(defReal), dimension(3), intent(in)                     :: currFaceNormal, centroid
    real(defReal), intent(in)                                   :: extraDistance, cellSpacing
    integer(shortInt), intent(in)                               :: faceIdx
    integer(shortInt), dimension(:), allocatable, intent(inout) :: intersectedFaceIdxs
    integer(shortInt)                                           :: i, j
    real(defReal), dimension(3)                                 :: currVertexCoords, min1, max1, currEdgeUnitVector
    real(defReal)                                               :: faceConst, currValue, vectorDotCentroid, extraDistance2 

    !------------------------------------------------------------------------------------------------
    ! Testing along face normal
    !------------------------------------------------------------------------------------------------
    vectorDotCentroid = dot_product(currFaceNormal, centroid)
    faceConst = faces % getFaceConst(faceIdx)

    if (.NOT. testIntervalIntersection(vectorDotCentroid - extraDistance, vectorDotCentroid + extraDistance, &
        -faceConst, -faceConst)) then
          !!!!!
          ! print*, "faceNormal"
          ! print*, vectorDotCentroid - extraDistance
          ! print*, vectorDotCentroid + extraDistance
          ! print*, -faceConst
          !!!!!
      return
    end if

    !------------------------------------------------------------------------------------------------
    ! Testing along crossProduct(each of edgeUnitVector, three coordinate basis)
    !------------------------------------------------------------------------------------------------
    ! loop over all edgesUnitVectors of the current face
    do i = 1, size(currFaceEdgeIdxs)
      currEdgeUnitVector = edges % getEdgeUnitvector(currFaceEdgeIdxs(i))
      
      ! loop over all vertices of the current face to find the min and max of polygon interval
      do j = 1, size(currVertexIdxs)
        currVertexCoords = vertices % getVertexCoordinates(currVertexIdxs(j))

        ! for coordinate basis = (1,0,0)
        currValue = currEdgeUnitVector(3)*currVertexCoords(2) - currEdgeUnitVector(2)*currVertexCoords(3)
        if (j == 1) then
          min1(1) = currValue
          max1(1) = currValue
        else
          if (currValue < min1(1)) then
            min1(1) = currValue
          elseif (currValue > max1(1)) then
            max1(1) = currValue
          end if 
        end if

        ! for coordinate basis = (0,1,0)
        currValue = - currEdgeUnitVector(3)*currVertexCoords(1) + currEdgeUnitVector(1)*currVertexCoords(3)
        if (j == 1) then
          min1(2) = currValue
          max1(2) = currValue
        else
          if (currValue < min1(2)) then
            min1(2) = currValue
          elseif (currValue > max1(2)) then
            max1(2) = currValue
          end if 
        end if

        ! for coordinate basis = (0,0,1)
        currValue =  currEdgeUnitVector(2)*currVertexCoords(1) - currEdgeUnitVector(1)*currVertexCoords(2)
        if (j == 1) then
          min1(3) = currValue
          max1(3) = currValue
        else
          if (currValue < min1(3)) then
            min1(3) = currValue
          elseif (currValue > max1(3)) then
            max1(3) = currValue
          end if 
        end if

      end do

      ! test intersections of interval
      vectorDotCentroid = currEdgeUnitVector(3)*centroid(2)-currEdgeUnitVector(2)*centroid(3)
      extraDistance2 = (abs(currEdgeUnitVector(3)) + abs(currEdgeUnitVector(2)))*cellSpacing*0.5
      if (.NOT. testIntervalIntersection(vectorDotCentroid - extraDistance2, vectorDotCentroid + extraDistance2, &
                                                min1(1), max1(1))) then
          !!!!!
          ! print*, "Crossx"
          !!!!!
          return
      end if

      vectorDotCentroid = currEdgeUnitVector(1)*centroid(3)-currEdgeUnitVector(3)*centroid(1)
      extraDistance2 = (abs(currEdgeUnitVector(3)) + abs(currEdgeUnitVector(1)))*cellSpacing*0.5
      if (.NOT. testIntervalIntersection(vectorDotCentroid - extraDistance2, vectorDotCentroid + extraDistance2, &
                                                min1(2), max1(2))) then
          !!!!!
          ! print*, "Crossy"
          !!!!!
          return
      end if

      vectorDotCentroid = currEdgeUnitVector(2)*centroid(1)-currEdgeUnitVector(1)*centroid(2)
      extraDistance2 = (abs(currEdgeUnitVector(2)) + abs(currEdgeUnitVector(1)))*cellSpacing*0.5
      if (.NOT. testIntervalIntersection(vectorDotCentroid - extraDistance2, vectorDotCentroid + extraDistance2, &
                                                min1(3), max1(3))) then
          !!!!!
          ! print*, "Crossz"
          !!!!!
          return
      end if

    end do

    !------------------------------------------------------------------------------------------------
    ! if survived to this point, then there is no separating axis. Hence, add the face index to the list
    !------------------------------------------------------------------------------------------------
    call append(intersectedFaceIdxs, faceIdx)
    
  end subroutine testFaceIntersectionCoarsest

  !!
  !!
  !!
  subroutine testFaceIntersectionNonCoarsest(vertices, edges, faces, currVertexIdxs, extraDistance, currFaceNormal, & 
                                             centroid, cellSpacing, faceIdx, currFaceEdgeIdxs, intersectedFaceIdxs, &
                                             arrIdx4Face, extraDistanceArr, removedFaceIdxs)
    class(vertexShelf), intent(in)                              :: vertices
    class(edgeShelf), intent(in)                                :: edges
    class(faceShelf), intent(in)                                :: faces
    integer(shortInt), dimension(:), intent(in)                 :: currVertexIdxs, currFaceEdgeIdxs
    real(defReal), dimension(3), intent(in)                     :: currFaceNormal, centroid
    real(defReal), intent(in)                                   :: extraDistance, cellSpacing
    integer(shortInt), intent(in)                               :: faceIdx, arrIdx4Face
    integer(shortInt), dimension(:), allocatable, intent(inout) :: intersectedFaceIdxs, removedFaceIdxs
    real(defReal), dimension(:), allocatable, intent(inout)     :: extraDistanceArr
    real(defReal)                                               :: faceConst, vectorDotCentroid 

    !------------------------------------------------------------------------------------------------
    ! Testing along face normal
    !------------------------------------------------------------------------------------------------
    vectorDotCentroid = dot_product(currFaceNormal, centroid)
    faceConst = faces % getFaceConst(faceIdx)

    if (.NOT. testIntervalIntersection(vectorDotCentroid - extraDistance, vectorDotCentroid + extraDistance, &
        -faceConst, -faceConst)) then

      ! Update relerant indices and return
      call eraseAt(intersectedFaceIdxs,arrIdx4Face)
      call eraseAt(extraDistanceArr,arrIdx4Face)
      call append(removedFaceIdxs, faceIdx)
      return
    end if

    !------------------------------------------------------------------------------------------------
    ! Testing along crossProduct(each of edgeUnitVector, three coordinate basis)
    !------------------------------------------------------------------------------------------------
    ! loop over all edgesUnitVectors of the current face
    ! do i = 1, size(currFaceEdgeIdxs)
    !   currEdgeUnitVector = edges % getEdgeUnitvector(currFaceEdgeIdxs(i))
      
    !   ! loop over all vertices of the current face to find the min and max of polygon interval
    !   do j = 1, size(currVertexIdxs)
    !     currVertexCoords = vertices % getVertexCoordinates(currVertexIdxs(j))

    !     ! for coordinate basis = (1,0,0)
    !     currValue = currEdgeUnitVector(3)*currVertexCoords(2) - currEdgeUnitVector(2)*currVertexCoords(3)
    !     if (j == 1) then
    !       min1(1) = currValue
    !       max1(1) = currValue
    !     else
    !       if (currValue < min1(1)) then
    !         min1(1) = currValue
    !       elseif (currValue > max1(1)) then
    !         max1(1) = currValue
    !       end if 
    !     end if

    !     ! for coordinate basis = (0,1,0)
    !     currValue = - currEdgeUnitVector(3)*currVertexCoords(1) + currEdgeUnitVector(1)*currVertexCoords(3)
    !     if (j == 1) then
    !       min1(2) = currValue
    !       max1(2) = currValue
    !     else
    !       if (currValue < min1(2)) then
    !         min1(2) = currValue
    !       elseif (currValue > max1(2)) then
    !         max1(2) = currValue
    !       end if 
    !     end if

    !     ! for coordinate basis = (0,0,1)
    !     currValue =  currEdgeUnitVector(2)*currVertexCoords(1) - currEdgeUnitVector(1)*currVertexCoords(2)
    !     if (j == 1) then
    !       min1(3) = currValue
    !       max1(3) = currValue
    !     else
    !       if (currValue < min1(3)) then
    !         min1(3) = currValue
    !       elseif (currValue > max1(3)) then
    !         max1(3) = currValue
    !       end if 
    !     end if

    !   end do

    !   ! test intersections of interval
    !   vectorDotCentroid = currEdgeUnitVector(3)*centroid(2)-currEdgeUnitVector(2)*centroid(3)
    !   extraDistance2 = (abs(currEdgeUnitVector(3)) + abs(currEdgeUnitVector(2)))*cellSpacing*0.5
    !   if (.NOT. testIntervalIntersection(vectorDotCentroid - extraDistance2, vectorDotCentroid + extraDistance2, &
    !                                             min1(1), max1(1))) then

    !       ! Update relerant indices and return
    !       call eraseAt(intersectedFaceIdxs,arrIdx4Face)
    !       call eraseAt(extraDistanceArr,arrIdx4Face)
    !       call append(removedFaceIdxs, faceIdx)
    !       return
    !   end if

    !   vectorDotCentroid = currEdgeUnitVector(1)*centroid(3)-currEdgeUnitVector(3)*centroid(1)
    !   extraDistance2 = (abs(currEdgeUnitVector(3)) + abs(currEdgeUnitVector(1)))*cellSpacing*0.5
    !   if (.NOT. testIntervalIntersection(vectorDotCentroid - extraDistance2, vectorDotCentroid + extraDistance2, &
    !                                             min1(2), max1(2))) then

    !       ! Update relerant indices and return
    !       call eraseAt(intersectedFaceIdxs,arrIdx4Face)
    !       call eraseAt(extraDistanceArr,arrIdx4Face)
    !       call append(removedFaceIdxs, faceIdx)
    !       return
    !   end if

    !   vectorDotCentroid = currEdgeUnitVector(2)*centroid(1)-currEdgeUnitVector(1)*centroid(2)
    !   extraDistance2 = (abs(currEdgeUnitVector(2)) + abs(currEdgeUnitVector(1)))*cellSpacing*0.5
    !   if (.NOT. testIntervalIntersection(vectorDotCentroid - extraDistance2, vectorDotCentroid + extraDistance2, &
    !                                             min1(3), max1(3))) then

    !       ! Update relerant indices and return
    !       call eraseAt(intersectedFaceIdxs,arrIdx4Face)
    !       call eraseAt(extraDistanceArr,arrIdx4Face)
    !       call append(removedFaceIdxs, faceIdx)
    !       return
    !   end if

    ! end do

    !------------------------------------------------------------------------------------------------
    ! if survived to this point, then there is no separating axis. Hence, do not remove the face index
    !------------------------------------------------------------------------------------------------
    
  end subroutine testFaceIntersectionNonCoarsest

  !!
  !!
  !! 
  subroutine testTwoIntersectedFaces(vertices, edges, faces, centroid, targetDistance, phi, phiCapital, faceIdxs)
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(in)                        :: edges
    class(faceShelf), intent(in)                        :: faces
    real(defReal), dimension(3), intent(in)             :: centroid
    real(defReal), intent(in)                           :: targetDistance
    integer(shortInt), intent(inout)                    :: phi, phiCapital
    integer(shortInt), dimension(2), intent(in)         :: faceIdxs
    integer(shortInt)                                   :: commonEdgeIdx
    integer(shortInt), dimension(:), allocatable        :: commonEdgeVertexIdxs
    real(defReal)                                       :: dummyConstant, dotProduct1, dotProduct2
    real(defReal), dimension(3)                         :: targetFaceNormal, edgeVertexCoord, c_1, &
                                                           edgeUnitVector, edgeVertex2Coord, &
                                                           tempVector1, tempVector2, chiCoord

    ! find the common edge index of the given two faces intersected by the cell 
    commonEdgeIdx = faces % findCommonEdgeIdx(faceIdxs(1), faceIdxs(2))

    ! if there is no common edge, exit the subroutine early (other combinations of faces to be tried later)
    if (commonEdgeIdx == 0) return

    ! calculate centre of intersection between the circumscribed ball and the plane parallel to the polygon
    commonEdgeVertexIdxs = edges % getEdgeVertexIdxs(commonEdgeIdx)
    targetFaceNormal = faces % getFaceNormal(faceIdxs(1))
    edgeVertexCoord = vertices % getVertexCoordinates(commonEdgeVertexIdxs(1))
    dummyConstant = dot_product(centroid - edgeVertexCoord, targetFaceNormal)
    c_1 = centroid - dummyConstant * targetFaceNormal

    ! calculate chi
    edgeUnitVector = edges % getEdgeUnitVector(commonEdgeIdx)
    dummyConstant = dot_product(c_1 - edgeVertexCoord, edgeUnitVector)
    chiCoord = edgeVertexCoord + dummyConstant * edgeUnitVector

    ! calculate distances to be compared
    edgeVertex2Coord = vertices % getVertexCoordinates(commonEdgeVertexIdxs(2))
    tempVector1 = chiCoord - edgeVertexCoord
    tempVector2 = chiCoord - edgeVertex2Coord
    dotProduct1 = dot_product(tempVector1, tempVector1)
    dotProduct2 = dot_product(tempVector2, tempVector2)

    ! make comparison between the two distances (dotProduct1/2)
    phi = merge(commonEdgeVertexIdxs(1), commonEdgeVertexIdxs(2), dotProduct1 < dotProduct2)
    if (dotProduct1 > targetDistance .and. dotProduct2 > targetDistance) phiCapital = commonEdgeIdx
    
  end subroutine testTwoIntersectedFaces

  !!
  !!
  !!
  pure subroutine constructMapSingleFace(faces, faceIdxs, phiCapital)
    class(faceShelf), intent(in)                        :: faces
    integer(shortInt), dimension(2), intent(in)         :: faceIdxs
    integer(shortInt), intent(inout)                    :: phiCapital
    integer(shortInt), dimension(:), allocatable        :: currFaceEdgeIdxs

    ! testing if the cell intersects with a single face
    ! (necessary to test this because this subroutine is called for all cells)
    ! no need to test if chi or both (phi and phiCapital) have been assigned from edge intersection or 
    ! polyhedron inclusion because these cells are not tested against face intersection (self % faceIdxs(1) = 0)
    if(faceIdxs(1) /= 0 .and. faceIdxs(2) == 0) then 
      ! (needs to be changed) (rather than saving candidate indices, write a subroutine that directly returns the first index of the array)
      ! (Does this subroutine has to save the array locally anyways?)
      currFaceEdgeIdxs = faces % getFaceEdgeIdxs(faceIdxs(1))

      ! (needs to be changed) (can find one in currFaceEdgeIdxs that are already been assigned to other cells)
      ! (Hence, we reduce the number of edges about which angles are calculated as well as memory if virtual grid is used)
      phiCapital = currFaceEdgeIdxs(1)

    end if

  end subroutine constructMapSingleFace

  !!
  !!
  !!
  subroutine coverFinitePrecision(faces, elements, centroid, elementIdx, chi)
    class(faceShelf), intent(in)                        :: faces
    class(elementShelf), intent(in)                     :: elements
    real(defReal), dimension(3), intent(in)             :: centroid
    integer(shortInt), intent(in)                       :: elementIdx
    integer(shortInt), intent(inout)                    :: chi
    integer(shortInt), dimension(:), allocatable        :: currElementFaceIdxs
    integer(shortInt)                                   :: i

    currElementFaceIdxs = elements % getElementFaceIdxs(elementIdx)

    do i = 1, size(currElementFaceIdxs)

      if (dot_product(faces % getFaceNormal(currElementFaceIdxs(i)), centroid) &
          + faces % getFaceConst(currElementFaceIdxs(i)) > 0) then
          ! if TRUE, then centroid of this cell lies outside of the polyhedron
          return
      end if 

    end do

    !print*, "TESTING", centroid
    ! If survived to this point, this cell is contained within a single mesh element, but
    ! chi mapping info was not assigned due to floating-point-finite-precision.
    chi = elementIdx

    !print*, "hahahoho"!!!!!

  end subroutine coverFinitePrecision

  !!
  !!
  !!
  subroutine setFaceParameters(faces)
    class(faceShelf), intent(inout)               :: faces
    integer(shortInt)                             :: i, j
    real(defReal)                                 :: extraDistance
    real(defReal), dimension(3)                   :: currFaceNormal
    integer(shortInt), dimension(3)               :: faceNormalSigns

    do i = 1, faces % getSize()

      ! Calculate and set extraDistance (to be multiplied by grid spacing)
      currFaceNormal = faces % getFaceNormal(i)
      extraDistance = (abs(currFaceNormal(1)) + abs(currFaceNormal(2)) + abs(currFaceNormal(3)))*0.5  
      call faces % setFaceExtraDistance(i, extraDistance)
      
      ! Calculate and set faceNormalSigns (to be multiplied by a half of grid spacing)
      do j = 1, 3
        if (currFaceNormal(j) > 0) then
          faceNormalSigns(j) = 1
        else
          faceNormalSigns(j) = -1
        end if
      end do
      call faces % setFaceNormalSigns(i, faceNormalSigns)

      ! Calculate and set const (constant = dot(any point on the plane ⊥ the face, face normal))
      call faces % setFaceConst(i, dot_product(faces % getFaceNormal(i), faces % getFaceCentroid(i))*(-1))

    end do

  end subroutine setFaceParameters

  !!
  !!
  !!
  subroutine setFaceParameters2(faces, spacing, n_layers)
    class(faceShelf), intent(inout)           :: faces
    real(defReal), dimension(:), intent(in)   :: spacing
    integer(shortInt)                         :: n_layers
    real(defReal), dimension(3)               :: currFaceNormal
    integer(shortInt)                         :: i, j
    real(defReal), dimension(:), allocatable  :: extraDistanceArr

    if (.NOT. allocated(extraDistanceArr)) allocate(extraDistanceArr(n_layers))

    do i = 1, faces % getSize()
      currFaceNormal = faces % getFaceNormal(i)
      extraDistanceArr(:) = abs(currFaceNormal(1)) + abs(currFaceNormal(2)) + abs(currFaceNormal(3))

      do j = 1, n_layers
        extraDistanceArr(j) = extraDistanceArr(j)*(0.5*spacing(j))
      end do

      call faces % setFaceExtraDistanceArr(i, extraDistanceArr, n_layers)
    end do

  end subroutine

  !!
  !! returns true if (dot_product(faceNormal, pointOfInterest-centroidOfFace) > 0)
  !! i.e. it returns true when the pointOfInterest lies in the non-owner element of the face
  function faceHalfSpaceTest(currFaceNormal, pointOfInterest, currFaceConst) result(outcome)
    real(defReal), dimension(3), intent(in)   :: currFaceNormal, pointOfInterest
    real(defReal), intent(in)                 :: currFaceConst
    logical(defBool)                          :: outcome 
  
    if (dot_product(currFaceNormal, pointOfInterest) + currFaceConst > 0) then
      outcome = .TRUE.
    else
      outcome = .FALSE.
    end if

  end function faceHalfSpaceTest

  !!
  !! Changes the order of element indices of a give face so that
  !! the order is [owner element, non-owner element]
  subroutine fixFaceElementIdxsOrder(faces, elements)
    class(faceShelf), intent(inout)               :: faces
    class(elementShelf), intent(in)               :: elements
    integer(shortInt)                             :: i
    integer(shortInt), dimension(:), allocatable  :: currFaceElementIdxs

    do i = 1, faces % getSize()

      if (allocated(currFaceElementIdxs)) deallocate(currFaceElementIdxs)
      currFaceElementIdxs = faces % getFaceElementIdxs(i)

      if (faceHalfSpaceTest(faces % getFaceNormal(i), &
          elements % getElementCentroid(currFaceElementIdxs(1)), faces % getFaceConst(i))) then

        ! Swap the order of indices if the first index is for the non-owner element of the given face.
        call faces % swapFaceElementIdxsOrder(i)

      end if

    end do

  end subroutine fixFaceElementIdxsOrder

  !!
  !!
  !!
  subroutine setGridFirstDimension(gridBounds_min, gridBounds_max, spacing, n_layers, n_xyz, &
                                   ratioFinest2Coarsest)
    real(defReal), dimension(3), intent(in)            :: gridBounds_min, gridBounds_max
    real(defReal), dimension(:), intent(inout)         :: spacing
    integer(shortInt), intent(in)                      :: n_layers
    integer(shortInt), dimension(:,:), intent(out)     :: n_xyz 
    integer(shortInt), intent(out)                     :: ratioFinest2Coarsest
    real(defReal)                                      :: gridSize, targetRatio
    integer(shortInt)                                  :: targetNFinest, targetNCoarsest, i, j  
    integer(shortInt), dimension(:), allocatable       :: candidateRatios, primeFactors
    
    ! Calculate basic parameters
    gridSize = gridBounds_max(1) - gridBounds_min(1)
    targetNFinest = ceiling(gridSize/spacing(n_layers))
    targetNCoarsest = ceiling(gridSize/spacing(1))
    targetRatio = targetNFinest/targetNCoarsest

    ! Depending on the number of layers used, retrieve candidate ratios between the number of cells in the finest and coarsest
    if (n_layers == 2) then
      ! Any integers
      candidateRatios = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
    elseif (n_layers == 3) then
      ! Integers that are produces of at least two prime numbers
      candidateRatios = [4,6,8,9,10,12,14,15,16,18,20,21,22,24,25]
    elseif (n_layers == 4) then
      ! Integers that are produces of at least three prime numbers
      candidateRatios = [8,12,16,18,20,24]
    end if

    ! Choose relevant ratio from candidateRatios
    if (targetratio <= 2) then
      ratioFinest2Coarsest = 2
    else

      loop: do i = size(candidateRatios), 1, -1
        if (targetRatio >= candidateRatios(i)) then
          exit loop
        end if
      end do loop
      ratioFinest2Coarsest = candidateRatios(i)

    end if

    ! Calculate the number of cells and cell spacing for the finest layer
    n_xyz(n_layers,1) = targetNFinest + &
                      (ratioFinest2Coarsest - mod(targetNFinest, ratioFinest2Coarsest))
    spacing(n_layers) = gridSize / n_xyz(n_layers,1)                                          

    ! Calculate the number of cells and cell spacings for the coarser layers
    if (n_layers > 2) then
      primeFactors = prime_factors(ratioFinest2Coarsest)

      j = 1
      do i = n_layers-1, 1, -1
        n_xyz(i,1) = n_xyz(i+1,1)/primeFactors(j)
        spacing(i) = gridSize / n_xyz(i,1)
        j = j + 1
      end do

    else ! i.e. if n_layers = 2

      n_xyz(1,1) = n_xyz(2,1)/ratioFinest2Coarsest
      spacing(1) = gridSize / n_xyz(1,1)

    end if

  end subroutine setGridFirstDimension

  !!
  !!
  !!
  subroutine setGridOtherDimensions(gridBounds_min, gridBounds_max, spacing, n_layers, n_xyz, &
                                   ratioFinest2Coarsest, dimension)
    real(defReal), dimension(3), intent(in)            :: gridBounds_min
    real(defReal), dimension(3), intent(inout)         :: gridBounds_max
    real(defReal), dimension(:), intent(in)            :: spacing
    integer(shortInt), intent(in)                      :: n_layers, ratioFinest2Coarsest, dimension
    integer(shortInt), dimension(:,:), intent(inout)   :: n_xyz 
    real(defReal)                                      :: gridSize, targetRatio, extraDistance
    integer(shortInt)                                  :: targetNFinest, targetNCoarsest, i, j  
    integer(shortInt), dimension(:), allocatable       :: primeFactors
    
    ! Calculate basic parameters
    gridSize = gridBounds_max(dimension) - gridBounds_min(dimension)
    targetNFinest = ceiling(gridSize/spacing(n_layers))
    targetNCoarsest = ceiling(gridSize/spacing(1))
    targetRatio = targetNFinest/targetNCoarsest

    ! Calculate the number of cells and cell spacing for the finest layer
    n_xyz(n_layers,dimension) = targetNFinest + &
                      (ratioFinest2Coarsest - mod(targetNFinest, ratioFinest2Coarsest))                                        

    ! Calculate the number of cells and cell spacings for the coarser layers
    if (n_layers > 2) then
      primeFactors = prime_factors(ratioFinest2Coarsest)

      j = 1
      do i = n_layers-1, 1, -1
        n_xyz(i,dimension) = n_xyz(i+1,dimension)/primeFactors(j)
        j = j + 1
      end do

    else ! i.e. if n_layers = 2

      n_xyz(1,dimension) = n_xyz(2,dimension)/ratioFinest2Coarsest

    end if

    ! Calculate the extra distance to be added to host the n_xyz(n_layers,dimension)
    extraDistance = n_xyz(n_layers,dimension)*spacing(n_layers) - gridSize
    ! if (extraDistance < 0) then
    !   call fatalError("SetGridOtherDimensions, cartesianInitProcedures.f90", &
    !   "the extra distance to be added to host the n_xyz(n_layers,dimension) is negative")
    ! elseif (extraDistance == 0) then
    !   return
    ! end if

    print*, extraDistance
    gridBounds_max(dimension) = gridBounds_max(dimension) + extraDistance


  end subroutine setGridOtherDimensions

  !!
  !!
  !!
  function prime_factors(n) result(factors)
      implicit none
      integer(shortInt), intent(in)  :: n
      integer(shortInt), allocatable :: factors(:)
      integer(shortInt)              :: num, i, count

      num = n

      do i = 2, num
          do while (mod(num, i) == 0)
              count = size(factors)
              call append(factors, i)
              num = num / i
          end do
      end do
  end function prime_factors

end module cartesianInitProcedures