module RNG_test

  use numPrecision
  use RNG_class, only : RNG
  use funit

  implicit none

contains

  !!
  !! Test random numbers
  !!
@Test
  subroutine testRN()
    integer(shortInt), parameter :: N = 1000
    type(RNG)                    :: pRNG
    real(defReal), dimension(N)  :: randomNumbers
    integer(shortInt)            :: i

    call pRNG % init(int(z'5c3a84c9', longInt))
    call pRNG % generate(randomNumbers)

    ! Check correcness
    @assertGreaterThanOrEqual(ONE, randomNumbers)
    @assertLessThanOrEqual(ZERO, randomNumbers)

  end subroutine testRN

  !!
  !! Test skip forward and backwards
  !!
@Test
  subroutine testSkip()
    type(RNG)         :: rand1
    type(RNG)         :: rand2
    integer(longInt)  :: seed
    real(defReal)     :: r_start, r2_start, r_end, r2_end
    integer(shortInt) :: i, N

    !! Initialise both RNGs to a nice number
    seed = int(z'5c3a84c9', longInt)
    call rand1 % init(seed)
    call rand2 % init(seed)

    !! Get initial random number
    call rand1 % generate(r_start)

    !! Move forward by 13456757 steps
    N = 13456757
    do i= 1, N
      call rand1 % generate(r_end)
    end do

    ! Skip 2nd generator forward
    call rand2 % skip(int(N, longInt))
    call rand2 % generate(r2_end)

    ! Skip 2nd generator backwards. Must be 1 more becouse we drew a RN from generator
    call rand2 % skip(-int(N + 1, longInt))
    call rand2 % generate(r2_start)

    ! Verify values
    @assertEqual(r_end, r2_end)
    @assertEqual(r_start, r2_start)

  end subroutine testSkip

end module RNG_test