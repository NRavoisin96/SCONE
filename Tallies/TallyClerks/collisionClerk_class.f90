module collisionClerk_class

  use dictionary_class,            only : dictionary
  use errors_mod,                  only : fatalError
  use nuclearDatabase_inter,       only : nuclearDatabase
  use numPrecision
  use outputFile_class,            only : outputFile
  use physicalParticle_inter,      only : physicalParticle
  use scoreMemory_class,           only : scoreMemory
  use tallyCodes
  use tallyClerk_inter,            only : tallyClerk, kill_super => kill
  use tallyFilter_inter,           only : tallyFilter
  use tallyFilterFactory_func,     only : new_tallyFilter
  use tallyMap_inter,              only : tallyMap
  use tallyMapFactory_func,        only : new_tallyMap
  use tallyResponseSlot_class,     only : tallyResponseSlot
  use tallyResult_class,           only : tallyResult, tallyResultArrays
  use transportObjectState_class,  only : transportObjectState
  use universalVariables

  implicit none
  private

  !!
  !! Collision estimator of reaction rates
  !! Calculates flux weighted integral from collisions
  !!
  !! Private Members:
  !!   filter   -> Space to store tally Filter
  !!   map      -> Space to store tally Map
  !!   response -> Array of responses
  !!   width    -> Number of responses (# of result bins for each map position)
  !!
  !! Interface
  !!   tallyClerk Interface
  !!
  !! SAMPLE DICTIONARY INPUT:
  !!
  !! myCollisionClerk {
  !!   type collisionClerk;
  !!   # filter { <tallyFilter definition> } #
  !!   # map    { <tallyMap definition>    } #
  !!   response (resName1 #resName2 ... #)
  !!   resName1 { <tallyResponse definition> }
  !!   #resNamew { <tallyResponse definition #
  !! }
  !!
  type, public, extends(tallyClerk) :: collisionClerk
    private
    ! Filter, Map & Vector of Responses
    character(nameLen), dimension(:), allocatable      :: responseNames
    class(tallyFilter), allocatable                    :: filter
    class(tallyMap), allocatable                       :: map
    type(tallyResponseSlot), dimension(:), allocatable :: responses

    ! Useful data
    integer(shortInt)  :: width = 0

    ! Settings
    logical(defBool)   :: handleVirtual = .true.

  contains
    ! Procedures used during build
    procedure :: init
    procedure :: kill
    procedure :: validReports
    procedure :: getResult
    procedure :: getSize

    ! File reports and check status -> run-time procedures
    procedure :: computeVolumeWeightedSum
    procedure :: flush
    procedure :: reportInColl

    ! Output procedures
    procedure :: display
    procedure :: print

  end type collisionClerk

contains
  !!
  !!
  !!
  function computeVolumeWeightedSum(self, memory) result(volumeWeightedSum)
    class(collisionClerk), intent(in) :: self
    type(scoreMemory), intent(in)     :: memory
    integer(longInt)                  :: address, baseAddress
    integer(shortInt)                 :: i
    real(defReal)                     :: volumeWeightedSum
    character(*), parameter           :: here = 'computeVolumeWeightedSum (collisionClerk_class.f90)'

    ! Initialise volumeWeightedSum = ZERO.
    volumeWeightedSum = ZERO

    ! Call fatalError if map is not allocated.
    if (.not. allocated(self % map)) call fatalError(here, 'Tally map is not allocated.')

    ! Accumulate sum.
    baseAddress = self % getMemAddress()
    do i = 1, self % map % bins(0)
      ! Calculate bin address
      address = baseAddress + self % width * (i - 1)
      volumeWeightedSum = volumeWeightedSum + memory % getScore(address) * self % map % getBinVolume(i)

    end do

  end function computeVolumeWeightedSum

  !!
  !! Initialise clerk from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(collisionClerk), intent(inout)          :: self
    class(dictionary), intent(in)                 :: dict
    character(nameLen), intent(in)                :: name
    character(nameLen), dimension(:), allocatable :: responseNames
    integer(shortInt)                             :: i

    ! Assign name
    call self % setName(name)

    ! Load filter
    if (dict % isPresent('filter')) call new_tallyFilter(self % filter, dict % getDictPtr('filter'))

    ! Load map
    if (dict % isPresent('map')) call new_tallyMap(self % map, dict % getDictPtr('map'))

    ! Get names of response dictionaries
    call dict % get(responseNames, 'response')

    ! Set width
    self % width = size(responseNames)

    ! Load responses
    allocate(self % responses(self % width))
    do i = 1, self % width
      call self % responses(i) % init(dict % getDictPtr(responseNames(i)))

    end do

    ! Handle virtual collisions
    call dict % getOrDefault(self % handleVirtual, 'handleVirtual', .true.)

    ! Load responseNames.
    self % responseNames = responseNames

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(collisioNClerk), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    if (allocated(self % responseNames)) deallocate(self % responseNames)

    ! Kill and deallocate filter
    if (allocated(self % filter)) then
      deallocate(self % filter)
    end if

    ! Kill and deallocate map
    if (allocated(self % map)) then
      call self % map % kill()
      deallocate(self % map)
    end if

    ! Kill and deallocate responses
    if (allocated(self % responses)) deallocate(self % responses)

    self % width   = 0
    self % handleVirtual = .true.

  end subroutine kill

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(collisionClerk), intent(in)           :: self
    integer(shortInt), dimension(:), allocatable :: validCodes

    validCodes = [inColl_CODE]

  end function validReports

  !!
  !!
  !!
  pure subroutine getResult(self, res, mem)
    class(collisionClerk), intent(in)              :: self
    class(tallyResult), allocatable, intent(inout) :: res
    type(scoreMemory), intent(in)                  :: mem
    character(nameLen)                             :: name
    integer(longInt)                               :: addr
    integer(shortInt)                              :: i, j, nBins
    real(defReal)                                  :: val, STD
    type(tallyResultArrays), pointer               :: resultsPtr

    ! Allocate result to tallyResultArrays
    ! Do not deallocate if already allocated to FMresult
    if (allocated(res)) then
      select type(res)
        class is (tallyResultArrays)
          ! Do nothing.

        class default
          ! Deallocate.
          deallocate(res)

      end select

    end if
    if (.not. allocated(res)) allocate(tallyResultArrays :: res)

    select type(ptr => res)
      type is (tallyResultArrays)
        resultsPtr => ptr

    end select

    ! Get name of clerk.
    name = self % getName()

    ! Get number of bins.
    nBins = 1
    if (allocated(self % map)) nBins = self % map % bins(0)

    ! Enforce shape of results array.
    if (allocated(resultsPtr % results)) then
      if (size(resultsPtr % results) /= self % width) then
        deallocate(resultsPtr % results)
        allocate(resultsPtr % results(self % width))

      end if

    else
      allocate(resultsPtr % results(self % width))

    end if

    ! Enforce shape of inner arrays.
    do i = 1, self % width
      associate(currentResults => resultsPtr % results(i))
        if (allocated(currentResults % values)) then
          if (size(currentResults % values) /= nBins) then
            deallocate(currentResults % values)
            allocate(currentResults % values(nBins))

          end if

        else
          allocate(currentResults % values(nBins))

        end if

        if (allocated(currentResults % standardDeviations)) then
          if (size(currentResults % standardDeviations) /= nBins) then
            deallocate(currentResults % standardDeviations)
            allocate(currentResults % standardDeviations(nBins))

          end if

        else
          allocate(currentResults % standardDeviations(nBins))

        end if

      end associate

    end do

    ! Load entries.
    addr = self % getMemAddress() - 1
    do i = 1, self % width
      associate(currentResults => resultsPtr % results(i))
        currentResults % clerkName = name
        currentResults % responseName = self % responseNames(i)
        do j = 1, nBins
          addr = addr + 1
          call mem % getResult(val, STD, addr)
          currentResults % values(j) = val
          currentResults % standardDeviations(j) = STD

        end do

      end associate

    end do

  end subroutine getResult

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(collisionClerk), intent(in) :: self
    integer(shortInt)                 :: S

    S = size(self % responses)
    if (allocated(self % map)) S = S * self % map % bins(0)

  end function getSize

  !!
  !!
  !!
  subroutine flush(self, memory)
    class(collisionClerk), intent(in) :: self
    type(scoreMemory), intent(inout)  :: memory
    integer(longInt)                  :: addr
    integer(shortInt)                 :: i, j, nBins

    nBins = 1
    if (allocated(self % map)) nBins = self % map % bins(0)

    ! Flush entries.
    addr = self % getMemAddress() - 1
    do i = 1, self % width
      do j = 1, nBins
        addr = addr + 1
        call memory % flush(addr)

      end do

    end do

  end subroutine flush

  !!
  !! Process incoming collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportInColl(self, p, virtual, xsData, mem)
    class(collisionClerk), intent(inout)     :: self
    class(physicalParticle), intent(in)      :: p
    logical(defBool), intent(in)             :: virtual
    class(nuclearDatabase), intent(inout)    :: xsData
    type(scoreMemory), intent(inout)         :: mem
    class(transportObjectState), pointer     :: currentStatePtr
    integer(shortInt)                        :: binIdx, i
    integer(longInt)                         :: addr
    real(defReal)                            :: scoreVal, flux
    character(*), parameter                  :: Here = 'reportInColl (collisionClerk_class.f90)'

    ! Return if collision is virtual but virtual collision handling is off
    if (virtual .and. .not. self % handleVirtual) return

    ! Get current particle state
    currentStatePtr => p % updateAndGetCurrentStatePtr()

    ! Check if within filter
    if (allocated(self % filter)) then
      if (self % filter % isFail(currentStatePtr)) return

    end if

    ! Find bin index
    binIdx = 1
    if (allocated(self % map)) binIdx = self % map % map(currentStatePtr)

    ! Return if invalid bin index
    if (binIdx == 0) return

    ! Calculate flux with the right cross section according to virtual collision handling
    if (self % handleVirtual) then
      flux = p % getWeight() / xsData % getTrackingXS(p, p % getMaterialIdx(), TRACKING_XS)

    else
      flux = p % getWeight() / xsData % getTotalMatXS(p, p % getMaterialIdx())

    end if

    ! Calculate bin address
    addr = self % getMemAddress() + self % width * (binIdx - 1) - 1

    ! Append all bins
    do i = 1, self % width
      call self % responses(i) % get(p, scoreVal, xsData)
      call mem % score(scoreVal * flux, addr + i)

    end do

  end subroutine reportInColl

  !!
  !! Display convergence progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(collisionClerk), intent(in)  :: self
    type(scoreMemory), intent(in)      :: mem

    print *, 'collisionClerk does not support display yet'

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(collisionClerk), intent(in)          :: self
    class(outputFile), intent(inout)           :: outFile
    type(scoreMemory), intent(in)              :: mem
    real(defReal)                              :: val, std
    integer(shortInt)                          :: i
    integer(shortInt), dimension(:), allocatable :: resArrayShape
    character(nameLen)                         :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    ! If collision clerk has map print map information. Then, write results.
    if (allocated(self % map)) then
      call self % map % print(outFile)
      resArrayShape = [size(self % responses), self % map % binArrayShape()]

    else
      resArrayShape = [size(self % responses)]

    end if

    ! Start array
    name ='Res'
    call outFile % startArray(name, resArrayShape)

    ! Print results to the file
    do i = 1, product(resArrayShape)
      call mem % getResult(val, std, self % getMemAddress() - 1 + i)
      call outFile % addResult(val, std)

    end do

    call outFile % endArray()
    call outFile % endBlock()

  end subroutine print

end module collisionClerk_class
