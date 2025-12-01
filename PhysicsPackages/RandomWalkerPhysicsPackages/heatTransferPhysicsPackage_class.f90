module heatTransferPhysicsPackage_class

  use coordList_class,            only : coordList
  use dictionary_class,           only : dictionary
  use element_class,              only : element, elementBox, elementIntersectionTestPayload, elementIntersectionTestResult, &
                                         inclusionTestResult, newElementIntersectionTestPayload
  use errors_mod,                 only : fatalError
  use face_class,                 only : face, orientatedFaceBox
  use genericProcedures,          only : append, areEqual, numToChar, rotateVector
  use geometry_inter,             only : geometry
  use geometryMesh_class,         only : getGeometryMeshPtr, geometryMesh
  use mesh_inter,                 only : mesh
  use numPrecision
  use outputFile_class,           only : outputFile
  use physicsPackage_inter,       only : init_super => init, initPhysicsPackagePayload, kill_super => kill, physicsPackage
  use randomWalker_class,         only : newRandomWalker, randomWalker
  use RNG_class,                  only : RNG
  use scalarField_inter,          only : getScalarFieldValue
  use tallyAdmin_class,           only : tallyAdmin
  use topologicalObject_inter,    only : topologicalObjectBox
  use transportOperatorWoS_class, only : transportOperatorWoS
  use universalVariables
  use unstructuredMesh_inter,     only : getCastUnstructuredMeshPtr, unstructuredMesh
  use vertex_class,               only : vertex, vertexBox

  implicit none
  private

  ! Parameters.
  integer(shortInt), parameter :: DEFAULT_N_WALKS_PER_BATCH = 1000
  real(defReal), parameter     :: DEFAULT_CONVERGENCE_CRITERION = 1.0e-3_defReal, DEFAULT_SURFACE_TOLERANCE = 1.0e-3_defReal

  ! Parameters (for now).
  real(defReal), parameter :: conductivity = 27.0e-2_defReal ! W cm⁻¹ K⁻¹

  !!
  !!
  !!
  type, public, extends(physicsPackage)      :: heatTransferPhysicsPackage
    private
    class(unstructuredMesh), pointer         :: unstructuredMeshPtr => null()
    integer(shortInt)                        :: nRuns = 0, nWalksPerBatch = 0
    real(defReal)                            :: convergenceCriterion = ZERO, surfaceTolerance = ZERO
    real(defReal), dimension(:), allocatable :: means, M2, parentErrors, parentSumOfScores, parentSumOfScoresSquared, variances
    type(RNG), pointer                       :: RNGPtr => null()
    type(tallyAdmin), pointer                :: tallyPtr => null()
    type(transportOperatorWoS)               :: transportOperator
  contains
    procedure :: collectSpecificResults
    procedure :: flushResults
    procedure :: getMeans
    procedure :: init
    procedure :: kill
    procedure :: run
    procedure :: walk
  end type heatTransferPhysicsPackage

contains
  !!
  !!
  !!
  subroutine collectSpecificResults(self, out)
    class(heatTransferPhysicsPackage), intent(in) :: self
    type(outputFile), intent(inout)               :: out

  end subroutine collectSpecificResults

  !!
  !!
  !!
  subroutine flushResults(self)
    class(heatTransferPhysicsPackage), intent(inout) :: self

    self % means = ZERO
    self % variances = ZERO
    self % nRuns = 0

  end subroutine flushResults

  !!
  !!
  !!
  pure function getMeans(self) result(means)
    class(heatTransferPhysicsPackage), intent(in) :: self
    real(defReal), dimension(:), allocatable      :: means

    if (allocated(self % means)) then
      means = self % means

    else
      allocate(means(0))

    end if

  end function getMeans

  !!
  !!
  !!
  subroutine init(self, payload)
    class(heatTransferPhysicsPackage), intent(inout) :: self
    class(initPhysicsPackagePayload), intent(in)     :: payload
    class(geometry), pointer                         :: geometryPtr
    class(mesh), pointer                             :: meshPtr
    class(unstructuredMesh), pointer                 :: unstructuredMeshPtr
    integer(shortInt)                                :: nParentElements
    character(*), parameter                          :: here = 'init (heatTransferPhysicsPackage_class.f90)'

    ! Initialise superclass.
    call init_super(self, payload)

    ! Load parameters.
    call payload % dict % getOrDefault(self % nWalksPerBatch, 'walksPerBatch', DEFAULT_N_WALKS_PER_BATCH)
    call payload % dict % getOrDefault(self % convergenceCriterion, 'convergenceCriterion', DEFAULT_CONVERGENCE_CRITERION)
    call payload % dict % getOrDefault(self % surfaceTolerance, 'surfaceTolerance', DEFAULT_SURFACE_TOLERANCE)

    ! Retrieve pointer to mesh geometry (hardcoded for now).
    geometryPtr => self % getGeometryPtr()

    ! Retrieve number of parent elements in the mesh geometry and allocate memory.
    meshPtr => geometryPtr % getMeshPtr(1)
    unstructuredMeshPtr => getCastUnstructuredMeshPtr(meshPtr)
    if (.not. associated(unstructuredMeshPtr)) call fatalError(here, 'Unable to retrieve unstructured mesh pointer.')
    self % unstructuredMeshPtr => unstructuredMeshPtr

    nParentElements = self % unstructuredMeshPtr % getParentElementsNumber()
    allocate(self % means(nParentElements), self % M2(nParentElements), self % parentErrors(nParentElements), &
             self % parentSumOfScores(nParentElements), self % parentSumOfScoresSquared(nParentElements), &
             self % variances(nParentElements))

    ! Initialise variables.
    self % means = ZERO
    self % M2 = ZERO
    self % parentErrors = INF
    self % parentSumOfScores = ZERO
    self % parentSumOfScoresSquared = ZERO
    self % variances = ZERO

    ! Initialise RNG.
    allocate(self % RNGPtr)
    call self % RNGPtr % init(self % getInitialSeed())

  end subroutine init

  !!
  !!
  !!
  subroutine kill(self)
    class(heatTransferPhysicsPackage), intent(inout) :: self

    ! Superclass.
    call kill_super(self)

    ! Local.
    self % unstructuredMeshPtr => null()
    self % nWalksPerBatch = 0
    self % convergenceCriterion = ZERO
    self % surfaceTolerance = ZERO
    if (allocated(self % means)) deallocate(self % means)
    if (allocated(self % M2)) deallocate(self % M2)
    if (allocated(self % parentErrors)) deallocate(self % parentErrors)
    if (allocated(self % parentSumOfScores)) deallocate(self % parentSumOfScores)
    if (allocated(self % parentSumOfScoresSquared)) deallocate(self % parentSumOfScoresSquared)
    if (associated(self % RNGPtr)) deallocate(self % RNGPtr)
    self % tallyPtr => null()
    call self % transportOperator % kill()

  end subroutine kill

  !!
  !!
  !!
  subroutine run(self)
    class(heatTransferPhysicsPackage), intent(inout) :: self
    integer(shortInt)                                :: elementIdx, i, j, nWalks, nWalksBatchStart
    integer(shortInt), dimension(:), allocatable     :: childrenIdxs
    real(defReal)                                    :: accumulatedValue, batchMean, batchVariance, previousMean, randomNumber
    type(coordList), pointer                         :: coordsPtr
    type(elementBox)                                 :: box, childBox
    type(randomWalker)                               :: walker
    type(RNG)                                        :: wRNG
    character(*), parameter                          :: here = 'run (heatTransferPhysicsPackage_class.f90)'

    print *, repeat("<>", 50)
    print *, "/\/\ HEAT TRANSFER CALCULATION /\/\"

    self % parentErrors = INF
    self % parentSumOfScores = ZERO
    self % parentSumOfScoresSquared = ZERO

    ! Loop over all regions.
    self % nRuns = self % nRuns + 1
    ! Initialise walks at centroid of each mesh elements (this should be handled by tally map).
    do i = 1, self % unstructuredMeshPtr % getElementsNumber()
      ! Retrieve current element. Check if it is a parent element.
      box = self % unstructuredMeshPtr % getElementBox(i)
      if (box % ptr % getParentIdx() == 0) then
        ! Loop until the error for this parent is below the convergence criterion.
        elementIdx = box % ptr % getIdx()
        nWalks = 0
        previousMean = ZERO
        do while(self % convergenceCriterion < self % parentErrors(i))
          
          nWalksBatchStart = nWalks
          ! Initialise parallel region here.
          !$omp parallel do &
          !$omp private(accumulatedValue, childBox, childrenIdxs, coordsPtr, j, randomNumber, wRNG, walker) &
          !$omp shared(box, elementIdx, nWalks, nWalksBatchStart, self) schedule(dynamic)
          
          ! Loop for a fixed number of walks.
          do j = 1, self % nWalksPerBatch
            ! Generate a new random walker and check if the current parent element is active.
            walker = newRandomWalker()
            wRNG = self % RNGPtr
            call walker % setRNGPtr(wRNG)
            call walker % strideRNG(nWalksBatchStart + j)
            coordsPtr => walker % getCoordsPtr()
            
            call coordsPtr % setNesting(1)
            call coordsPtr % setMeshIdx(1, 1)
            if (box % ptr % getIsActive()) then
              ! Current parent is active. Initialise the coordinates to the centroid of this element.
              call coordsPtr % setPosition(box % ptr % getCentroid(), 1)
              call coordsPtr % setElementIdx(elementIdx, 1)

            else
              ! Pick a child element at random and set the position to its centroid.
              childrenIdxs = box % ptr % getChildrenIdxs()

              call walker % generateRandomNumber(randomNumber)

              childBox = self % unstructuredMeshPtr % getElementBox(childrenIdxs(int(randomNumber * size(childrenIdxs)) + 1))
              call coordsPtr % setPosition(childBox % ptr % getCentroid(), 1)
              call coordsPtr % setElementIdx(childBox % ptr % getIdx(), 1)

            end if
            ! Walk now. Move this to transport operator.
            ! call self % transportOperator % transport(walker, self % tallyPtr)
            call self % walk(walker)

            ! Update statistics here. Retrieve accumulated value and update scores.
            accumulatedValue = walker % getAccumulatedValue()

            ! Increment number of walks.
            !$omp atomic update
            nWalks = nWalks + 1
            !$omp end atomic

            if (accumulatedValue == ZERO) cycle
            associate(sum => self % parentSumOfScores(elementIdx), sumOfSquares => self % parentSumOfScoresSquared(elementIdx))
              !$omp atomic update
              sum = sum + accumulatedValue

              !$omp atomic update
              sumOfSquares = sumOfSquares + accumulatedValue * accumulatedValue

            end associate

          end do
          !$omp end parallel do

          ! Update standard error for current parent.
          if (2000 < nWalks) then
            associate(sum => self % parentSumOfScores(elementIdx), sumOfSquares => self % parentSumOfScoresSquared(elementIdx))
              ! Compute mean and variance for the current batch.
              batchMean = sum / nWalks
              batchVariance = (sumOfSquares - sum * batchMean) / (nWalks - 1)

              ! Update global statistics.
              previousMean = self % means(i)
              self % means(i) = self % means(i) + (batchMean - self % means(i)) / self % nRuns

              if (self % nRuns == 1) then
                self % variances(i) = batchVariance / nWalks

              else
                self % M2(i) = self % M2(i) + (batchMean - previousMean) * (batchMean - self % means(i))
                self % variances(i) = self % M2(i) / (self % nRuns - 1)

              end if

              ! Convergence check.
              if (self % means(i) /= ZERO) then
                if (self % nRuns == 1) then
                  self % parentErrors(i) = sqrt(batchVariance / nWalks) / abs(self % means(i))

                else
                  self % parentErrors(i) = sqrt(self % variances(i) / self % nRuns) / abs(self % means(i))

                end if

              end if

            end associate

          end if

        end do

        print *, 'Temperature:', self % means(i), '+/-', sqrt(self % variances(i))
        print *, 'Number of walks:', nWalks

      end if

    end do

    print *
    print *, "\/\/ END OF HEAT TRANSFER CALCULATION \/\/"
    print *

  end subroutine run

  !!
  !!
  !!
  subroutine walk(self, walker)
    class(heatTransferPhysicsPackage), intent(inout)      :: self
    type(randomWalker), intent(inout)                     :: walker
    integer(shortInt)                                     :: boundaryCondition, i, j, nSharingElements
    real(defReal)                                         :: minDistance, mu, phi, remainingDistance, transmissionProbability, &
                                                             valueToAccumulate, dist
    real(defReal), dimension(2)                           :: coefficients
    real(defReal), dimension(3)                           :: outwardNormal, r, u
    type(coordList), pointer                              :: coordsPtr
    type(face), pointer                                   :: internalFacePtr
    type(element), pointer                                :: chosenElementPtr, elementPtr, neighbourElementPtr
    type(elementBox)                                      :: box
    type(elementIntersectionTestResult)                   :: faceIntersectionResults
    type(orientatedFaceBox), dimension(:), allocatable    :: faceBoxes
    type(RNG), pointer                                    :: RNGPtr
    type(topologicalObjectBox), dimension(:), allocatable :: sharingElements
    character(*), parameter                               :: here = 'walk (heatTransferPhysicsPackage_class.f90)'

    ! Initialise isDead = .false.
    coordsPtr => walker % getCoordsPtr()

    walkLoop: do
      ! Sample initial direction on the unit sphere and set it.
      RNGPtr => walker % getRNGPtr()
      call RNGPtr % generateMu(mu)
      call RNGPtr % generatePhi(phi)
      u = rotateVector([ONE, ZERO, ZERO], mu, phi)
      call coordsPtr % setDirection(u, 1)
      
      ! Compute distances to all faces of the current element.
      box = self % unstructuredMeshPtr % getElementBox(coordsPtr % getLowestElementIdx())
      faceBoxes = box % ptr % getOrientatedFaces()
      
      r = coordsPtr % getPosition(1)
      minDistance = INF
      elementPtr => box % ptr
      internalFacePtr => null()
      chosenElementPtr => elementPtr
      do i = 1, size(faceBoxes)
        ! Get the boundary condition for the current face. Do not include distance if it is a Neumann boundary condition.
        boundaryCondition = faceBoxes(i) % face % ptr % getBoundaryCondition(TEMPERATURE_BCs)
        if (boundaryCondition == ZERO_TEMPERATURE_GRADIENT_BC) cycle
        dist = dot_product(faceBoxes(i) % face % ptr % getCentroid() - r, faceBoxes(i) % outwardNormal)

        ! If walker is below surface tolerance for this face, snap it onto the face and apply boundary condition.
        if (dist < self % surfaceTolerance) then
          select case(boundaryCondition)
            case(FIXED_TEMPERATURE_BC)
              ! Retrieve the boundary value of the face and kill walker.
              call walker % accumulateValue(faceBoxes(i) % face % ptr % getBoundaryValue(TEMPERATURE_BCs))
              return

            case(INTERNAL_TEMPERATURE_BC)
              ! Associate internal face pointer if distance is smaller than minimum found.
              if (dist < minDistance) then
                internalFacePtr => faceBoxes(i) % face % ptr
                outwardNormal = faceBoxes(i) % outwardNormal

              end if

          end select

        end if
        minDistance = min(minDistance, dist)

      end do

      ! If the internal face pointer is associated, we need to recompute the minimum distance.
      if (associated(internalFacePtr)) then
        ! Snap walker on the closest internal face.
        r = r + minDistance * outwardNormal
        call coordsPtr % setPosition(r, 1)

        ! Get the elements sharing the face.
        sharingElements = internalFacePtr % getSharingElements()
        nSharingElements = size(sharingElements)
        if (nSharingElements /= 2) call fatalError(here, 'Internal face is not associated with two elements.')
        do j = 1, 2
          select type(ptr => sharingElements(j) % ptr)
            type is(element)
              if (.not. associated(elementPtr, ptr)) neighbourElementPtr => ptr
            
            class default
              call fatalError(here, 'Element with index: '//numToChar(ptr % getIdx())//' is not an element.')

          end select

        end do

        ! Compute transmission probability.
        call RNGPtr % generate(transmissionProbability)

        if (transmissionProbability < HALF) then
          ! Walker is reflected back into original element. Check that it points in the correct direction.
          chosenElementPtr => elementPtr
          do while(ZERO <= dot_product(u, outwardNormal))
            call RNGPtr % generateMu(mu)
            call RNGPtr % generatePhi(phi)
            u = rotateVector([ONE, ZERO, ZERO], mu, phi)

          end do

        else
          ! Walker transmits into the neighbouring element. Check that it points in the correct direction.
          chosenElementPtr => neighbourElementPtr
          do while(dot_product(u, outwardNormal) <= ZERO)
            call RNGPtr % generateMu(mu)
            call RNGPtr % generatePhi(phi)
            u = rotateVector([ONE, ZERO, ZERO], mu, phi)

          end do

        end if
        ! Exit loop.
        call coordsPtr % setDirection(u, 1)
        call coordsPtr % setElementIdx(chosenElementPtr % getIdx(), 1)

        minDistance = INF
        faceBoxes = elementPtr % getOrientatedFaces()
        do i = 1, size(faceBoxes)
          boundaryCondition = faceBoxes(i) % face % ptr % getBoundaryCondition(TEMPERATURE_BCs)
          if (boundaryCondition == ZERO_TEMPERATURE_GRADIENT_BC .or. associated(faceBoxes(i) % face % ptr, internalFacePtr)) cycle
          dist = dot_product(faceBoxes(i) % face % ptr % getCentroid() - r, faceBoxes(i) % outwardNormal)
          minDistance = min(minDistance, dist)

        end do

        faceBoxes = neighbourElementPtr % getOrientatedFaces()
        do i = 1, size(faceBoxes)
          boundaryCondition = faceBoxes(i) % face % ptr % getBoundaryCondition(TEMPERATURE_BCs)
          if (boundaryCondition == ZERO_TEMPERATURE_GRADIENT_BC .or. associated(faceBoxes(i) % face % ptr, internalFacePtr)) cycle
          dist = dot_product(faceBoxes(i) % face % ptr % getCentroid() - r, faceBoxes(i) % outwardNormal)
          minDistance = min(minDistance, dist)

        end do

      end if

      if (minDistance < ZERO) call fatalError(here, 'Minimum distance is negative.')

      ! If reached here, simply transport the walker until it travels minDistance.
      remainingDistance = minDistance
      transportLoop: do
        ! Compute the distance to the next intersection.
        call chosenElementPtr % intersects_Ray(newElementIntersectionTestPayload(coordsPtr % getPosition(1), &
                                                                                 coordsPtr % getDirection(1), &
                                                                                 remainingDistance, .true., .true.), &
                                               faceIntersectionResults)

        ! If no intersection is detected, simply transport the walker to its end destination and exit.
        if (.not. faceIntersectionResults % intersects) then
          call coordsPtr % moveLocal(remainingDistance, coordsPtr % getNesting())
          exit transportLoop

        else
          ! Get the boundary condition associated with the face that was hit.
          call coordsPtr % moveLocal(faceIntersectionResults % d, coordsPtr % getNesting())
          boundaryCondition = faceIntersectionResults % intersectedFace % ptr % getBoundaryCondition(TEMPERATURE_BCs)
          if (boundaryCondition == ZERO_TEMPERATURE_GRADIENT_BC) then
            u = coordsPtr % getDirection(1)
            call faceIntersectionResults % intersectedFace % ptr % flipDirection(u)
            call coordsPtr % setDirection(u, 1)

          end if
          remainingDistance = remainingDistance - faceIntersectionResults % d
          if (areEqual(remainingDistance, ZERO)) exit transportLoop

        end if

      end do transportLoop

      ! Accumulate heat source. Set valueToAccumulate = ZERO in case the heat source field does not exist.
      valueToAccumulate = getScalarFieldValue(nameHeatSource, ZERO, coordsPtr, SIXTH * minDistance * minDistance / conductivity)
      call walker % accumulateValue(valueToAccumulate)

    end do walkLoop

  end subroutine walk

end module heatTransferPhysicsPackage_class