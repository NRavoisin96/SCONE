module temperatureClerk_class

  use dictionary_class,     only : dictionary
  use numPrecision
  use outputFile_class,     only : outputFile
  use scoreMemory_class,    only : scoreMemory
  use tallyClerk_inter,     only : tallyClerk
  use tallyCodes
  use tallyMap_inter,       only : tallyMap
  use tallyMapFactory_func, only : new_tallyMap

  implicit none
  private

  !!
  !!
  !!
  type, public, extends(tallyClerk) :: temperatureClerk
    private
    class(tallyMap), allocatable    :: map
  contains
    procedure :: display
    procedure :: getSize
    procedure :: init
    procedure :: print
    procedure :: validReports
  end type temperatureClerk

contains
  !!
  !!
  !!
  subroutine display(self, mem)
    class(temperatureClerk), intent(in) :: self
    type(scoreMemory), intent(in)       :: mem

    print *, 'temperatureClerk does not support display yet.'

  end subroutine display

  !!
  !!
  !!
  elemental function getSize(self) result(S)
    class(temperatureClerk), intent(in) :: self
    integer(shortInt)                   :: S

    S = 1
    if (allocated(self % map)) S = self % map % bins(0)

  end function getSize

  !!
  !!
  !!
  subroutine init(self, dict, name)
    class(temperatureClerk), intent(inout) :: self
    class(dictionary), intent(in)          :: dict
    character(nameLen), intent(in)         :: name

    ! Assign name.
    call self % setName(name)

    ! Load map.
    if (dict % isPresent('map')) call new_tallyMap(self % map, dict % getDictPtr('map'))

  end subroutine init

  !!
  !!
  !!
  subroutine print(self, outFile, mem)
    class(temperatureClerk), intent(in)          :: self
    class(outputFile), intent(inout)             :: outFile
    type(scoreMemory), intent(in)                :: mem
    real(defReal)                                :: val, std
    integer(shortInt)                            :: i
    integer(shortInt), dimension(:), allocatable :: resArrayShape
    character(nameLen)                           :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    ! If collision clerk has map print map information. Then, write results.
    if (allocated(self % map)) then
      call self % map % print(outFile)
      resArrayShape = [1, self % map % binArrayShape()]

    else
      resArrayShape = [1]

    end if

    ! Start array
    name = 'Res'
    call outFile % startArray(name, resArrayShape)

    ! Print results to the file
    do i = 1, product(resArrayShape)
      call mem % getResult(val, std, self % getMemAddress() - 1 + i)
      call outFile % addResult(val,std)

    end do

    call outFile % endArray()
    call outFile % endBlock()

  end subroutine print

  !!
  !!
  !!
  function validReports(self) result(validCodes)
    class(temperatureClerk), intent(in)          :: self
    integer(shortInt), dimension(:), allocatable :: validCodes

    validCodes = [temperature_CODE]

  end function validReports

end module temperatureClerk_class