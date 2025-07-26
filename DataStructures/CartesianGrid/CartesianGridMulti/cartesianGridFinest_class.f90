module cartesianGridFinest_class

  use vertexShelf_class,                       only : vertexShelf
  use edgeShelf_class,                         only : edgeShelf
  use faceShelf_class,                         only : faceShelf
  use elementShelf_class,                      only : elementShelf
  use cartesianGridSubLayer_inter,             only : cartesianGridSubLayer
  use cartesianCellFinest_class,               only : cartesianCellFinest
  use numPrecision   
  use genericProcedures,                       only : crossProduct
  use cartesianGenericProcedures

  implicit none
  private


  !!
  !!
  !!
  type, public, extends(cartesianGridSubLayer)                :: cartesianGridFinest
    private
    type(cartesianCellFinest), dimension(:,:,:), allocatable  :: grid

  contains

    ! Build procedures
    procedure                                    :: init
    procedure                                    :: constructMapping
    procedure                                    :: sortAngles
    procedure                                    :: setGridIsOutsideMesh
    ! Runtime procedures
    procedure                                    :: getGridChi
    procedure                                    :: getGridPhi
    procedure                                    :: getGridPhiCapital
  end type cartesianGridFinest

contains

  !!
  !!
  !!
  pure subroutine init(self, vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                       currLayer, candidateElementIdxs, gridBoundsMin, alpha, wStar)
    class(cartesianGridFinest), intent(inout)           :: self
    class(vertexShelf), intent(in)                      :: vertices
    class(edgeShelf), intent(inout)                     :: edges
    class(faceShelf), intent(inout)                     :: faces
    class(elementShelf), intent(in)                     :: elements
    real(defReal), dimension(:), intent(in)             :: spacing, spacingInv
    integer(shortInt), dimension(:,:), intent(in)       :: n_xyz
    integer(shortInt), intent(in)                       :: n_layers, currLayer
    integer(shortInt), dimension(:), intent(in)         :: candidateElementIdxs
    real(defReal), dimension(3), intent(in)             :: gridBoundsMin
    real(defReal), intent(in)                           :: alpha, wStar
    integer(shortInt), dimension(3)                     :: localNxyz
    integer(shortInt)                                   :: i

    ! calculate local number of cells 
    do i = 1, 3
      localNxyz(i) = n_xyz(currLayer,i)/n_xyz(currLayer-1,i)
    end do

    ! construct mapping and start constructing the full mapping
    allocate(self % grid(localNxyz(1), localNxyz(2), localNxyz(3)))
    call self % constructMapping(vertices, edges, faces, elements, spacing, spacingInv, n_xyz, n_layers, &
                                 currLayer, candidateElementIdxs, gridBoundsMin, localNxyz, alpha, wStar)
    call self % sortAngles(edges, faces, localNxyz)
    call self % setGridIsOutsideMesh(localNxyz)

  end subroutine

  !!
  !!
  !!
  pure subroutine constructMapping(self, vertices, edges, faces, elements, spacing, spacingInv, n_xyz, &
                      n_layers, currLayer, candidateElementIdxs, gridBoundsMin, localNxyz, alpha, wStar)
    class(cartesianGridFinest), intent(inout)             :: self
    class(vertexShelf), intent(in)                        :: vertices
    class(edgeShelf), intent(inout)                       :: edges
    class(faceShelf), intent(inout)                       :: faces
    class(elementShelf), intent(in)                       :: elements
    real(defReal), dimension(:), intent(in)               :: spacing, spacingInv
    integer(shortInt), dimension(:,:), intent(in)         :: n_xyz
    integer(shortInt), intent(in)                         :: n_layers, currLayer
    integer(shortInt), dimension(:), intent(in)           :: candidateElementIdxs
    real(defReal), dimension(3), intent(in)               :: gridBoundsMin
    integer(shortInt), dimension(3), intent(in)           :: localNxyz
    real(defReal), intent(in)                             :: alpha, wStar
    integer(shortInt)                                     :: h, i, j, k, l
    integer(shortInt), dimension(:), allocatable          :: currVertexIdxs, currElementFaceIdxs, currFaceEdgeIdxs, &
                                                             candElementEdgeIdxs, candElementFaceIdxs
    real(defReal)                                         :: circumscribedBallRadius, targetDistance, a, &
                                                             extraDistance
    real(defReal), dimension(3)                           :: centroid, currEdgeVector, currFaceNormal
    real(defReal), dimension(:,:), allocatable            :: faceNormalSigns

    !----------------------------------------------------------------------------------------------
    ! edge interesection tests
    !----------------------------------------------------------------------------------------------
    !calculate constants for edge interesection tests
    circumscribedBallRadius = sqrt(3.0d0)*(spacing(n_layers))/2
    targetDistance = (wStar) / (1 + SIN(alpha))

    do h = 1, size(candidateElementIdxs)
      candElementEdgeIdxs = elements % getElementEdgeIdxs(candidateElementIdxs(h))

      do i = 1, size(candElementEdgeIdxs)

          ! calculate edge-only-dependent properties
          currVertexIdxs = edges % getEdgeVertexIdxs(candElementEdgeIdxs(i))
          currEdgeVector = (edges % getEdgeUnitvector(candElementEdgeIdxs(i)))* &
                           (edges % getEdgeLength(candElementEdgeIdxs(i)))
          a = dot_product(currEdgeVector, currEdgeVector)

          !Loop over all cartesian cells in the box and test if each cell intersects with the edge
          !(needs to be changed) (k and l can be a function of j e.g. k = datum + slope*j so that box is narrowed down)
          do j = 1, localNxyz(1)
              do k = 1, localNxyz(2)
                  do l = 1, localNxyz(3)

                      ! (needs to be changed) (store centroid info)
                      centroid(1) = (gridBoundsMin(1)) + (spacing(n_layers)) * (j-0.5)
                      centroid(2) = (gridBoundsMin(2)) + (spacing(n_layers)) * (k-0.5)
                      centroid(3) = (gridBoundsMin(3)) + (spacing(n_layers)) * (l-0.5)

                      call self % grid(j,k,l) % cellTestEdgeIntersection(vertices, edges, candElementEdgeIdxs(i), &
                               circumscribedBallRadius, targetDistance, centroid, currEdgeVector, currVertexIdxs, a)
                      
                  end do 
              end do    
          end do

      end do

    end do

    !----------------------------------------------------------------------------------------------
    ! polyhedron inclusion tests
    !----------------------------------------------------------------------------------------------
    allocate(faceNormalSigns(3, 2))
    do i = 1, size(candidateElementIdxs)

        ! calculate element-only-dependent properties
        currElementFaceIdxs = elements % getElementFaceIdxs(candidateElementIdxs(i))

        deallocate(faceNormalSigns)
        allocate(faceNormalSigns(3, size(currElementFaceIdxs)))
        do j = 1, size(currElementFaceIdxs)
          currFaceNormal = faces % getFaceNormal(currElementFaceIdxs(j))
          do k = 1, 3
            if (currFaceNormal(k) > 0) then
              faceNormalSigns(k, j) = 1
            else
              faceNormalSigns(k, j) = -1
            end if
          end do
        end do
        faceNormalSigns = faceNormalSigns * (spacing(n_layers))/2


        !Loop over all cartesian cells in the box and test if each cell is entirely included in the polyhedron
        !(needs to be changed) (k and l can be a function of j e.g. k = datum + slope*j so that box is narrowed down)
        do j = 1, localNxyz(1)
            do k = 1, localNxyz(2)
                do l = 1, localNxyz(3)

                    ! (needs to be changed) (store centroid info)
                    centroid(1) = (gridBoundsMin(1)) + (spacing(n_layers)) * (j-0.5)
                    centroid(2) = (gridBoundsMin(2)) + (spacing(n_layers)) * (k-0.5)
                    centroid(3) = (gridBoundsMin(3)) + (spacing(n_layers)) * (l-0.5)

                    call self % grid(j,k,l) % cellTestPolyhedronInclusion(faces, currElementFaceIdxs, centroid, &
                                                                          faceNormalSigns, candidateElementIdxs(i))

                end do 
            end do    
        end do

    end do

    !----------------------------------------------------------------------------------------------
    ! face intersection tests
    !----------------------------------------------------------------------------------------------
    targetDistance = targetDistance**2

    do h = 1, size(candidateElementIdxs)
      candElementFaceIdxs = abs(elements % getElementFaceIdxs(candidateElementIdxs(h)))

      do i = 1, size(candElementFaceIdxs)

        ! calculate face-only-dependent properties
        currVertexIdxs = faces % getFaceVertexIdxs(candElementFaceIdxs(i))
        currFaceEdgeIdxs = faces % getFaceEdgeIdxs(candElementFaceIdxs(i))
        currFaceNormal = faces % getFaceNormal(candElementFaceIdxs(i))
        extraDistance = (abs(currFaceNormal(1)) + abs(currFaceNormal(2)) + abs(currFaceNormal(3))) &
                        * (spacing(n_layers)) * 0.5               

        !Loop over all cartesian cells in the box and test if each cell intersect with the current face
        !(needs to be changed) (k and l can be a function of j e.g. k = datum + slope*j so that box is narrowed down)
        do j = 1, localNxyz(1)
            do k = 1, localNxyz(2)
                do l = 1, localNxyz(3)

                    ! (needs to be changed) (store centroid info)
                    centroid(1) = (gridBoundsMin(1)) + (spacing(n_layers)) * (j-0.5)
                    centroid(2) = (gridBoundsMin(2)) + (spacing(n_layers)) * (k-0.5)
                    centroid(3) = (gridBoundsMin(3)) + (spacing(n_layers)) * (l-0.5)

                    call self % grid(j,k,l) % cellTestFaceIntersection(vertices, edges, faces, &
                                              currVertexIdxs, extraDistance, currFaceNormal, &
                                              centroid, spacing(n_layers), candElementFaceIdxs(i), &
                                              currFaceEdgeIdxs, targetDistance)

                end do 
            end do    
        end do

      end do

    end do

    ! for cells that intersec more than one face, mappings constructions are performed during face intersection test
    ! for those that intersect exactly one face, mappings constructions are performed here.

    ! loop over all cartesian cells and call relevant subroutine
    do i = 1, localNxyz(1)
      do j = 1, localNxyz(2)
        do k = 1, localNxyz(3)
          call self % grid(i,j,k) % cellConstructMapSingleFace(faces)
        end do
      end do 
    end do

  end subroutine constructMapping

  !!
  !!
  !!
  pure subroutine sortAngles(self, edges, faces, localNxyz)
    class(cartesianGridFinest), intent(inout)           :: self
    class(edgeShelf), intent(inout)                     :: edges
    class(faceShelf), intent(in)                        :: faces
    integer(shortInt), dimension(3), intent(in)         :: localNxyz
    integer(shortInt)                                   :: i, j, k, l, m, n, currPhiCapital, v_e, pointerIdx, &
                                                           outer2LoopSize, currEdgeVertex1Idx
    real(defReal), dimension(3)                         :: currEdgeUnitVector, localBasis1, localBasis2, currUnitVector
    integer(shortInt), dimension(:), allocatable        :: currEdgeFaceIdxs, currFaceEdgeIdxs, faceIdxsArray, &
                                                           elementIdxsArray
    integer(shortInt), dimension(2)                     :: currEdgeVertexIdxs, currVertexIdxs, face1ElementIdxs, &
                                                           face2ElementIdxs
    real(defReal)                                       :: x, y, thetaHat
    real(defReal), dimension(:), allocatable            :: anglesArray

    ! loop through all cartesian cells so that only [edge index]s, where there exists at least one [cell index] s.t. 
    ! phiCapital([cell index]) = [edge index], are used.
    ! (needs to be changed) (since all edges (most likely) are assigned for phiCapital mapping anyways, just loop through all edges?)
    do i = 1, localNxyz(1)
      do j = 1, localNxyz(2)
        do k = 1, localNxyz(3)

          currPhiCapital = self % grid(i,j,k) % getPhiCapital()

          ! continue only if the current cell contains valid edge index mapping of phiCapital
          if (currPhiCapital /= 0) then
            ! continue only if the current edge index (= phiCapital) has not been used for sorting angles yet
            if (.NOT. edges % isAllocatedEdgeAnglesArray(currPhiCapital)) then

              !---------------------------------------------------------------------------------------------------------------
              ! construct 2D local cooridnate system (localBasis1,localBasis2) on the plane whose normal is given as the current edge's unit vector
              ! and contains the second vertex of the edge.
              !---------------------------------------------------------------------------------------------------------------
              currEdgeUnitVector = edges % getEdgeUnitVector(currPhiCapital)
              
              ! construct localBasis1
              if (abs(currEdgeUnitVector(1)) <= abs(currEdgeUnitVector(2)) .AND. &
                  abs(currEdgeUnitVector(1)) <= abs(currEdgeUnitVector(3))) then
                    localBasis1 = [0.0d0, currEdgeUnitVector(3), -currEdgeUnitVector(2)]
              elseif (abs(currEdgeUnitVector(2)) <= abs(currEdgeUnitVector(3))) then 
                    localBasis1 = [-currEdgeUnitVector(3), 0.0d0, currEdgeUnitVector(1)]
              else
                    localBasis1 = [currEdgeUnitVector(2), -currEdgeUnitVector(1), 0.0d0]
              end if
              localBasis1 = localBasis1 / norm2(localBasis1)

              ! construct localBasis2 (cross product gives the normalised vector)
              localBasis2 = crossProduct(currEdgeUnitVector, localBasis1)

              ! store localBasis1 and localBasis2 to each associated edge
              call edges % setEdgeLocalBasis1(currPhiCapital, localBasis1)
              call edges % setEdgeLocalBasis2(currPhiCapital, localBasis2)

              !---------------------------------------------------------------------------------------------------------------
              ! construct (unsorted) arrays for angles and associated elementIdxs
              !---------------------------------------------------------------------------------------------------------------
              ! retrieve relevant information
              currEdgeFaceIdxs = edges % getEdgeFaceIdxs(currPhiCapital)
              currEdgeVertexIdxs = edges % getEdgeVertexIdxs(currPhiCapital)
              v_e = currEdgeVertexIdxs(2)
              currEdgeVertex1Idx = currEdgeVertexIdxs(1)

              ! initialise arrays for angle and face index
              if (allocated(anglesArray)) deallocate(anglesArray)
              if (allocated(faceIdxsArray)) deallocate(faceIdxsArray)
              if (allocated(elementIdxsArray)) deallocate(elementIdxsArray)
              allocate(anglesArray(size(currEdgeFaceIdxs)))
              allocate(faceIdxsArray(size(currEdgeFaceIdxs)))
              allocate(elementIdxsArray(size(currEdgeFaceIdxs)))

              ! loop through all faces attached to the current edge (= currPhiCapital)
              outer1: do l = 1, size(currEdgeFaceIdxs)
                currFaceEdgeIdxs = faces % getFaceEdgeIdxs(currEdgeFaceIdxs(l))
                
                ! loop through all edges of the current face (of the current edge = currPhiCapital)
                inner1: do m = 1, size(currFaceEdgeIdxs)
                  currVertexIdxs = edges % getEdgeVertexIdxs(currFaceEdgeIdxs(m))

                  ! test if this current edge of the current face is the one that is connected to edge = currPhiCapital
                  ! and that this current edge is not currPhiCapital itself
                  if (ANY(currVertexIdxs == v_e)) then
                    if (.NOT. ANY(currVertexIdxs == currEdgeVertex1Idx)) then

                      ! retrieve and correct the orientation of the unit vector of the edge
                      if (currVertexIdxs(1) == v_e) then
                        currUnitVector = edges % getEdgeUnitVector(currFaceEdgeIdxs(m))
                      else
                        currUnitVector = edges % getEdgeUnitVector(currFaceEdgeIdxs(m))*(-1)
                      end if

                      ! calculate the 2D local coordinates
                      x = dot_product(currUnitVector, localBasis1)
                      y = dot_product(currUnitVector, localBasis2)

                      ! calculate pseudo angle
                      thetaHat = SIGN(1 - (x / (abs(x) + abs(y))), y)

                      ! add the calculated angle and face index to each corresponding arrays
                      anglesArray(l) = thetaHat
                      faceIdxsArray(l) = currEdgeFaceIdxs(l)

                      ! once the calculation is performed for current face, move on to the next face of the current edge (= currPhiCapital)
                      ! (l = l + 1)
                      exit inner1

                    end if
                  end if

                end do inner1

              end do outer1

              ! perform index sorting on anglesArray and faceIdxsArray
              call sortPairs(anglesArray, faceIdxsArray)

              ! construct a sorted array for element index
              ! pointer of index for sorted array assignment
              pointerIdx = 1

              ! loop through all neighbouring faces
              outer2LoopSize = size(faceIdxsArray)
              outer2: do l = 1, outer2LoopSize
                face1ElementIdxs = faces % getFaceElementIdxs(faceIdxsArray(l))
                face2ElementIdxs = faces % getFaceElementIdxs(faceIdxsArray(mod(l, outer2LoopSize) + 1))

                ! if boundary face, set the second element index = 0, which can be used when constructing elementIdxsArray. 
                if (size(faces % getFaceElementIdxs(faceIdxsArray(l))) == 1) face1ElementIdxs(2) = 0
                if (size(faces % getFaceElementIdxs(faceIdxsArray(mod(l, outer2LoopSize) + 1))) == 1) face2ElementIdxs(2) = 0

                ! find the common element indicies from 2 X arrays of size 2
                ! (needs checking) (check and throw fatal error if there are two common element idxs-for boundary edge?)
                ! (There can be only one element attached to the face for boundary faces?)
                ! (If it is a boundary edge with only two faces attached, no need to calculate this)
                ! (During face intersection, if a cell intersects one or two boundary faces (of the same element))
                ! (then we can directly assign "chi" rather than phi and phiCapital)
                ! (and this removes the issue disscussed above because only [edge idx] s.t. there exists [cell idx] s.t.)
                ! (phiCapital([cell idx]) = [edgeIdx] are used for this subroutine)
                ! (this applies for boundary edges with any number of faces attached to it)
                ! (needs to be changed) (due to the possible acceleration above)
                middle2: do m = 1, 2
                  inner2: do n = 1, 2
                    if (face1ElementIdxs(m) == face2ElementIdxs(n)) then
                      elementIdxsArray(pointerIdx) = face1ElementIdxs(m)
                      pointerIdx = pointerIdx + 1
                      exit middle2
                    end if
                  end do inner2
                end do middle2

              end do outer2

              !!!
              ! if any of the faces is a boundary face, set the last element index = 0
              ! (needs to be changed) (it might not be the last element index in some cases?)
              ! (initially set all indices = 0 so that if face1ElementIdxs(m) != face2ElementIdxs(n))
              ! (for all m and n, space outside the unstructured mesh has element index = 0)
              ! (for this move "pointerIdx = pointerIdx + 1" out of loops) (have tried but not work?)
              ! if (flagBoundaryFace) elementIdxsArray(size(elementIdxsArray)) = 0
              ! solved problem by: if (size(faces % getFaceElementIdxs(faceIdxsArray(l))) == 1) face1ElementIdxs(2) = 0
              !!!

              ! pass and set elementIdxsArray and anglesArray to each corresponding edge
              call edges % setEdgeAnglesArray(currPhiCapital, anglesArray)
              call edges % setEdgeElementIdxsArray(currPhiCapital, elementIdxsArray)

            end if
          end if

        end do
      end do 
    end do

  end subroutine sortAngles

  !!
  !!
  !! 
  pure subroutine setGridIsOutsideMesh(self, localNxyz)
    class(cartesianGridFinest), intent(inout)              :: self
    integer(shortInt), dimension(3), intent(in)            :: localNxyz
    integer(shortInt)                                      :: i, j, k

    ! loop over all cartesian cells and call relevant subroutine
    do i = 1, localNxyz(1)
      do j = 1, localNxyz(2)
        do k = 1, localNxyz(3)
          call self % grid(i,j,k) % setIsOutsideMesh()
        end do
      end do 
    end do

  end subroutine setGridIsOutsideMesh

  !!
  !!
  !! 
  pure function getGridChi(self, baseIntegerCoord, shift, mask, currLayer) result(chi)
    class(cartesianGridFinest), intent(in)              :: self
    integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
    integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
    integer(shortInt), intent(in)                       :: currLayer
    integer(shortInt)                                   :: chi
    integer(shortInt), dimension(3)                     :: cellIdxs

    cellIdxs = getLocalIdx(baseIntegerCoord, shift(currLayer,:), mask(currLayer,:))

    chi = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getChi()

  end function getGridChi

  !!
  !!
  !! 
  pure function getGridPhi(self, baseIntegerCoord, shift, mask, currLayer) result(Phi)
    class(cartesianGridFinest), intent(in)              :: self
    integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
    integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
    integer(shortInt), intent(in)                       :: currLayer
    integer(shortInt)                                   :: phi
    integer(shortInt), dimension(3)                     :: cellIdxs

    cellIdxs = getLocalIdx(baseIntegerCoord, shift(currLayer,:), mask(currLayer,:))

    phi = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getphi()

  end function getGridPhi

  !!
  !!
  !! 
  pure function getGridPhiCapital(self, baseIntegerCoord, shift, mask, currLayer) result(phiCapital)
    class(cartesianGridFinest), intent(in)              :: self
    integer(shortInt), dimension(3), intent(in)         :: baseIntegerCoord
    integer(shortInt), dimension(:,:), intent(in)       :: shift, mask
    integer(shortInt), intent(in)                       :: currLayer
    integer(shortInt)                                   :: phiCapital
    integer(shortInt), dimension(3)                     :: cellIdxs

    cellIdxs = getLocalIdx(baseIntegerCoord, shift(currLayer,:), mask(currLayer,:))

    phiCapital = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getPhiCapital()

  end function getGridPhiCapital

end module cartesianGridFinest_class
