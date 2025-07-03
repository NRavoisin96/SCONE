module longIntMap_test
  
  use funit
  use longIntMap_class, only : longIntMap
  use numPrecision

  implicit none


@testCase
  type, extends(TestCase) :: test_longIntMap
    private
    type(longIntMap)      :: map
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_longIntMap

  ! Parameters
  integer(shortInt), parameter :: VAL1 = 9
  integer(shortInt), parameter :: VAL2 = 1
  integer(shortInt), parameter :: VAL3 = -1
  integer(shortInt), parameter :: VAL4 =  huge(VAL1)
  integer(shortInt), parameter :: VAL5 = -huge(VAL1)
  integer(shortInt), parameter :: VAL6 =  133316666 ! Number of demons of hell
  integer(shortInt), parameter :: VAL7 = -133316666

  integer(longInt), parameter :: KEY1 = 0
  integer(longInt), parameter :: KEY2 = 1
  integer(longInt), parameter :: KEY3 = -1
  integer(longInt), parameter :: KEY4 =  huge(VAL1)
  integer(longInt), parameter :: KEY5 = -huge(VAL1)
  integer(longInt), parameter :: KEY6 =  133316666 ! Number of demons of hell
  integer(longInt), parameter :: KEY7 = -133316666

contains

  !!
  !! Sets up test_intMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_longIntMap), intent(inout) :: this

    ! Load entries
     call this % map % add(KEY1, VAL1)
     call this % map % add(KEY2, VAL2)
     call this % map % add(KEY3, VAL3)
     call this % map % add(KEY4, VAL4)
     call this % map % add(KEY5, VAL5)
     call this % map % add(KEY6, VAL6)
     call this % map % add(KEY7, VAL7)

  end subroutine setUp

  !!
  !! Kills test_intMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_longIntMap), intent(inout) :: this

    call this % map % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test retrieving values from the map.
  !! All values are present
  !!
@Test
  subroutine testRetrieval(this)
    class(test_longIntMap), intent(inout) :: this
    integer(shortInt)                     :: temp

    temp = this % map % get(KEY1)
    @assertEqual(VAL1, temp)

    temp = this % map % get(KEY2)
    @assertEqual(VAL2, temp)

    temp = this % map % get(KEY3)
    @assertEqual(VAL3, temp)

    temp = this % map % get(KEY4)
    @assertEqual(VAL4, temp)

    temp = this % map % get(KEY5)
    @assertEqual(VAL5, temp)

    temp = this % map % get(KEY6)
    @assertEqual(VAL6, temp)

    temp = this % map % get(KEY7)
    @assertEqual(VAL7, temp)

  end subroutine testRetrieval

  !!
  !! Test Retrivial with deletions
  !!
@Test
  subroutine testRetrievalWithDel(this)
    class(test_longIntMap), intent(inout) :: this
    integer(shortInt)                     :: temp

    ! Delate Elements
    call this % map % del(KEY1)
    call this % map % del(KEY2)
    call this % map % del(KEY3)
    call this % map % del(KEY4)
    call this % map % del(KEY5)

    ! Obtain some elements
    temp = this % map % get(KEY6)
    @assertEqual(VAL6, temp)

    temp = this % map % get(KEY7)
    @assertEqual(VAL7, temp)

  end subroutine testRetrievalWithDel

  !!
  !! Test getting entery with default for absent keys
  !!
@Test
  subroutine testGetOrDefault(this)
    class(test_longIntMap), intent(inout) :: this
    integer(shortInt)                     :: temp
    integer(shortInt)                     :: default = 8

    ! Google says its Chinese Lucky number. We do need luck
    default = 8

    ! Retrieve present entries. Do all cases to spot possible
    ! small differences between get and getOrDefault implementations
    temp = this % map % getOrDefault(KEY1, default)
    @assertEqual(VAL1, temp)

    temp = this % map % getOrDefault(KEY2, default)
    @assertEqual(VAL2, temp)

    temp = this % map % getOrDefault(KEY3, default)
    @assertEqual(VAL3, temp)

    temp = this % map % getOrDefault(KEY4, default)
    @assertEqual(VAL4, temp)

    temp = this % map % getOrDefault(KEY5, default)
    @assertEqual(VAL5, temp)

    temp = this % map % getOrDefault(KEY6, default)
    @assertEqual(VAL6, temp)

    temp = this % map % getOrDefault(KEY7, default)
    @assertEqual(VAL7, temp)

    ! Get default. (Yes I love Tolkien - MAK)
    @assertEqual(default, this % map % getOrDefault(1973_longInt, default))

  end subroutine testGetOrDefault

  !!
  !! Test deletions
  !!
@Test
  subroutine testDel(this)
    class(test_longIntMap), intent(inout) :: this
    integer(shortInt)                     :: temp

    ! Before deletion
    temp = this % map % get(KEY6)
    @assertEqual(VAL6, temp)

    call this % map % del(KEY6)

    ! After deletion
    temp = this % map % getOrDefault(KEY6, VAL1)
    @assertEqual(VAL1, temp)

  end subroutine testDel

  !!
  !! Test Looping over map
  !!
@Test
  subroutine testLooping(this)
    class(test_longIntMap), intent(inout)      :: this
    integer(shortInt), dimension(6), parameter :: VALS = [VAL1, VAL2, VAL3, VAL4, VAL5, VAL7]
    integer(longInt), dimension(6), parameter  :: KEYS = [KEY1, KEY2, KEY3, KEY4, KEY5, KEY7]
    integer(longInt), dimension(6)             :: KEYS_PAST
    integer(shortInt)                          :: counter, it, tVal
    integer(longInt)                           :: tKey

    ! Initialise parameters
    KEYS_PAST = 7
    counter = 0

    ! Delete Entry 6 from map
    call this % map % del(KEY6)

    ! Loop over remaining elements
    it = this % map % begin()
    do while(it /= this % map % end())
      tVal = this % map % atVal(it)
      tKey = this % map % atKey(it)

      @assertTrue(any(VALS == tVal),"Wrong Value")
      @assertTrue(any(KEYS == tKey),"Wrong Key")

      ! Change value and read it
      call this % map % atSet(7, it)
      @assertEqual(7, this % map % atVal(it), "Ups. Wrong Value after change")

      ! Make shure that KEY is not getting repeated
      @assertFalse(any(KEYS_PAST(1:counter) == tKey),"Repeated KEY")
      counter = counter + 1
      KEYS_PAST(counter) = tKey

      it = this % map % next(it)
    end do

    @assertEqual(6, counter)

  end subroutine testLooping

  !!
  !! Test Looping Edge Cases
  !!
@Test
  subroutine testLoopingEdgeCases(this)
    class(test_longIntMap), intent(inout) :: this
    type(longIntMap)                      :: locMap
    integer(shortInt)                     :: it

    ! Loop over uninitialised map
    it = locMap % begin()
    do while(it /= locMap % end())
      @assertTrue(.false., "Should not enter the loop")
      it = locMap % next(it)
    end do

    ! Loop over empty map
    call locMap % init(8)
    it = locMap % begin()
    do while(it /= locMap % end())
      @assertTrue(.false., "Should not enter the loop")
      it = locMap % next(it)
    end do

    ! Loop over map with deleted elements
    call locMap % add(1_longInt, 3)
    call locMap % add(2_longInt, 3)
    call locMap % add(7_longInt, 3)

    call locMap % del(1_longInt)
    call locMap % del(2_longInt)
    call locMap % del(7_longInt)

    it = locMap % begin()
    do while(it /= locMap % end())
      @assertTrue(.false., "Should not enter the loop")
      it = locMap % next(it)
    end do

  end subroutine testLoopingEdgeCases

  !!
  !! Test putting new entry under exoisting key
  !!
@Test
  subroutine testOverwriting(this)
    class(test_longIntMap), intent(inout) :: this
    integer(shortInt)                     :: temp

    temp = 7

    ! Store over existing keyword
    call this % map % add(KEY4, temp)

    ! Verify correctness
    @assertEqual(temp, this % map % get(KEY4))

  end subroutine testOverwriting

  !!
  !! Test getting length
  !!
@Test
  subroutine testGetLength(this)
    class(test_longIntMap), intent(inout) :: this

    @assertEqual(7, this % map % length())

    ! Deleate some elements
    call this % map % del(KEY1)
    call this % map % del(KEY6)

    @assertEqual(5, this % map % length())

  end subroutine testGetLength

  !!
  !! Test getOrDefaoult from uninitialised map
  !!
@Test
  subroutine testGetOrDefaultUninitialised(this)
    class(test_longIntMap), intent(inout) :: this
    type(longIntMap)                      :: locMap

    @assertEqual(7, locMap % getOrDefault(3_longInt, 7))

  end subroutine testGetOrDefaultUninitialised

  !!
  !! Test getOrDefault from initialised but empty map
  !!
@Test
  subroutine testGetOrDefaultEmpty(this)
    class(test_longIntMap), intent(inout) :: this
    type(longIntMap)                      :: locMap

    ! Pure empty map
    call locMap % init(2)
    @assertEqual(7, locMap % getOrDefault(3_longInt, 7))

    ! Map with deleted element
    call locMap % add(3_longInt, 7)
    call locMap % del(3_longInt)
    call locMap % init(2)
    @assertEqual(2, locMap % getOrDefault(3_longInt, 2))

  end subroutine testGetOrDefaultEmpty

end module longIntMap_test