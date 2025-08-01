module hashFunctions_func

  use numPrecision
  use genericProcedures, only : fatalError
  use iso_fortran_env,   only : int8,  int32, int64

  implicit none
  private


  !!
  !! Public interface for FNV_1 hash function
  !!
  public :: FNV_1
  interface FNV_1
    module procedure FNV_1_int32
    !module procedure FNV_1_int64
  end interface

  !!
  !! Public interface for Knuth Hash
  !!
  public :: knuthHash
  interface knuthHash
    module procedure knuthHash_int32
    module procedure knuthHash_int64
  end interface

contains

  !!
  !! Multiplicative Hash outlined in D. Knuth "The Art of Computer Programming"
  !!
  !! Implementation uses a 32-bit hash. To avoid complier warning associated with integer
  !! overflow on multiplication a temporary integer of higher precision is used and
  !! then converted to shortInt by ignoring higher bits.
  !!
  !! Ratio was changed from the popular "golden ratio" to a value based on sphere
  !! maximum packing fraction in Euclidian geometry = pi/3/sqrt(2) [3180339487 for 32 bits]
  !! NOTE: Ratio was changed for fun -> no particular reason
  !!
  !!
  elemental function knuthHash_int32(key, m) result(hash)
    integer(int32), intent(in) :: key
    integer(int32), intent(in) :: m
    integer(int32)             :: hash, m_loc
    integer(int32), parameter  :: constant = int(Z'9E3779B1', int32)

    ! Constrain m to range 1-31
    m_loc = max(1, m)
    m_loc = min(31, m_loc)

    ! Calculate constant * key.
    hash = constant * key

    ! Keep m uppermost bits.
    hash = shiftr(hash, 32 - m_loc)

  end function knuthHash_int32

  !!
  !! Multiplicative Hash outlined in D. Knuth "The Art of Computer Programming"
  !!
  !! Implementation uses a 64-bit hash. To avoid complier warning associated with integer
  !! overflow on multiplication a temporary integer of higher precision is used and
  !! then converted to shortInt by ignoring higher bits.
  !!
  !! Ratio was changed from the popular "golden ratio" to a value based on sphere
  !! maximum packing fraction in Euclidian geometry = pi/3/sqrt(2) [3180339487 for 32 bits]
  !! NOTE: Ratio was changed for fun -> no particular reason
  !!
  !!
  elemental function knuthHash_int64(key, m) result(hash)
    integer(int64), intent(in) :: key
    integer(int32), intent(in) :: m
    integer(int32)             :: hash, m_loc
    integer(int64)             :: hash_loc
    integer(int64), parameter  :: constant = int(Z'9E3779B97F4A7C15', int64), mask = int(Z'7FFFFFFFFFFFFFFF', int64)
    ! Constrain m to range 1-31
    m_loc = max(1, m)
    m_loc = min(31, m_loc)

    ! Calculate constant * key.
    hash_loc = iand(constant * key, mask)

    ! Extract uppermost bits.
    hash_loc = shiftr(hash_loc, 64 - m_loc)

    ! Convert back to a 32-bit integer.
    hash = int(hash_loc, int32)

  end function knuthHash_int64

  !!
  !! Implementation of Fowler-Noll-Vo 1 hash function for 32 bit integer
  !! Prime and offset taken from:
  !!   http://www.isthe.com/chongo/tech/comp/fnv/index.html
  !!
  !! Note:
  !!   Since we cannot have unsigned integer we need to set the bits to signed integer.
  !!   To avoid compiler warnings the conversion was done externally.
  !!
  !! TODO: Hash functions should be offloaded to C
  !!
  pure subroutine FNV_1_int32(key, hash)
    character(*), intent(in)    :: key
    integer(int32), intent(out) :: hash
    integer(int32), parameter   :: FNV_prime = 16777619_int32
    integer(int32), parameter   :: FNV_offset = -2128831035_int32 !int(z'811c9dc5',int32)
    integer(int32)              :: bajt
    integer(shortInt)           :: i

    ! There is a problem in gfortran 8.3, where if the loop
    ! starts from 2 (1st iteration is unrolled) incorrect code
    ! is generated
    hash = FNV_offset
    do i= 1, len(key)
      bajt = iachar(key(i:i), int32)
      hash = hash * FNV_prime
      hash = ieor(hash,bajt)
    end do

  end subroutine FNV_1_int32

!  !!
!  !! Implementation of Fowler-Noll-Vo 1 hash function for 64 bit integer
!  !! Prime and offset taken from: 14695981039346656037
!  !!   http://www.isthe.com/chongo/tech/comp/fnv/index.html
!  !!
!  subroutine FNV_1_int64(key, hash)
!    character(*), intent(in)    :: key
!    integer(int64), intent(out) :: hash
!    integer(int64), parameter   :: FNV_prime  = 1099511628211_int64
!    integer(int64), parameter   :: FNV_offset = transfer(z'cbf29ce484222325',int64)
!    integer(int64)             :: bajt
!    integer(shortInt)          :: i
!
!    bajt = ichar(key(1:1),int64)
!    hash = FNV_offset
!    hash = hash * FNV_prime
!    hash = ieor(hash,bajt)
!
!    do i=2,len(key)
!      bajt = ichar(key(i:i),int64)
!      hash = hash * FNV_prime
!      hash = ieor(hash, bajt)
!
!    end do
!
!  end subroutine FNV_1_int64


end module hashFunctions_func
