module ratint
    
  use limb_class
  use numPrecision
  use, intrinsic :: iso_fortran_env
  use, intrinsic :: ieee_arithmetic

  implicit none 

  type ratint_t
    type(limb_t) :: p, q
  end type ratint_t

  interface convert_ieee 
    module procedure convert_ieee64
    module procedure convert_ieee64Vector
  end interface convert_ieee

  interface isZero 
    module procedure isZero_flat 
    module procedure isZero_Vector
  end interface isZero

  interface operator (+)
    module procedure addpure 
    module procedure addmixedL
    module procedure addmixedR
  end interface operator (+)

  interface operator (-)
    module procedure subtractpure
    module procedure subtractmixedL
    module procedure subtractmixedR
  end interface operator (-)

  interface operator (*)
    module procedure multiplypure
    module procedure multiplymixedL
    module procedure multiplymixedR
    module procedure multiplyVectorL
    module procedure multiplyVectorR
  end interface operator (*)

  interface operator (/)
    module procedure dividepure
    module procedure dividemixedL
    module procedure dividemixedR
  end interface operator (/)

  interface assignment (=)
    module procedure assignpure
  end interface assignment (=)

  interface operator (==)
    module procedure equality
  end interface operator (==)

  interface operator (>)
    module procedure gneq
  end interface operator (>)

  interface operator (>=)
    module procedure geq
  end interface operator (>=)

  interface dot_product 
    module procedure dot_product_ratint
  end interface dot_product

  interface swapSign 
    module procedure swapSign_flat 
    module procedure swapSign_vector 
  end interface swapSign
  
  contains

  !! Enforces the denominator to always be positive by moving any negative
  !! sign into the numerator. All comparison operators cross-multiply by q 
  !! and are only valid while this holds.
  pure subroutine normaliseSign(r)
    type(ratint_t), intent(inout) :: r

    if(checkInvalidRatint(r)) return
    if(r % q % sign == -1) then
      r % q % sign = 1
      if(.not. limbiszero(r % p)) r % p % sign = -r % p % sign

    end if

  end subroutine normaliseSign

  pure function def_ratint(n, d, s) result(r)
      integer(longInt), intent(in) :: n, d, s
      type(ratint_t) :: r 

      r%p = initlimb(n) * s
      r%q = initlimb(d)

      !r = simplify(r)

  end function def_ratint

  pure function def_ratint_large() result (r)
      type(ratint_t) :: r 

      r%p = initlimb() 
      r%p%limbs(1) = 0 
      r%p%limbs(2) = 0
      r%p%limbs(3) = 0
      r%p%limbs(4) = 1 
      r%p%front = 4 
      r%q = initlimb(1_8)

  end function def_ratint_large

  pure function initratint_vector() result(r)
      type(ratint_t), dimension(3) :: r 
      integer(shortInt) :: i 

      do i =1, 3 
          r(i) = convert_int(0_8)
      end do


  end function initratint_vector


  pure function convert_int(n) result(r)
      integer(longInt), intent(in) :: n 
      type(ratint_t) :: r 

      r%q = initlimb(1_8) 
      r%p = initlimb(n)

  end function convert_int

  pure subroutine swapSign_flat(a)
      type(ratint_t), intent(inout) :: a 

      a%p%sign = a%p%sign * (-1)
      
  end subroutine swapSign_flat


  pure subroutine swapSign_vector(a)
      type(ratint_t), dimension(:), intent(inout) :: a 
      integer(shortInt) :: i

      do i =1, size(a)
          call swapSign_flat(a(i))
      end do
      
  end subroutine swapSign_vector

  
  pure function convert_ieee64(n) result(r)
      real(defReal), intent(in) :: n 
      type(ratint_t) :: r 
      real(defReal) :: n1
      integer(longInt) :: i, shift

      real(defReal) :: frac, exp
      type(ratint_t) :: expratint, fracratint

      !print *, n
      frac = fraction(n)
      exp = exponent(n)


      shift = 0


      
      do i=0, 52

          if (int(frac, 8)*1_defReal == frac) then 
              shift = i 
              exit 
          end if 
          shift = i+1

          frac = frac * 2
      end do
    
      
      ! Simplification based on powers of 2
      if (exp > 0) then 
          if (exp >= shift) then 
              exp = exp - shift 
              shift = 0 
          else 
              shift = shift - exp 
              exp = 0 
          end if 
      end if


      fracratint%p = initlimb(int(frac, 8)*1_8) 
      fracratint%q = initlimb(2**shift)

      

      ! print * , sign(1, floor(exp))
      if (exp < 0) then 
          expratint%p = initlimb(1_8)
          expratint%q = initlimb(2**int(abs(floor(exp)), 8))
      else 

          expratint%p = initlimb(2**int(exp, 8))
          expratint%q = initlimb(1_8)
      end if 



      r = expratint * fracratint

      if (sign(1.0_defReal,n) == -1.0_defReal) then 
          r%p%sign = -1
      end if




  end function convert_ieee64

  pure function convert_ieee64Vector(n) result(r)
      real(defReal), dimension(:), intent(in) :: n 
      type(ratint_t), allocatable:: r(:)
      integer(shortInt) :: i, height 

      height = size(n)
      allocate(r(height))

      do i=1, height 
          r(i) = convert_ieee64(n(i))
      end do


  end function convert_ieee64Vector

  pure function checkInvalidRatint(r) result(x)
      type(ratint_t), intent(in) :: r 
      logical :: x 

      x = .false.

      if (checkInvalid(r%p) .or. checkInvalid(r%q)) then 
          x = .true. 
          return
      end if 

  end function checkInvalidRatint


  pure subroutine setInvalidRatint(r)
      type(ratint_t), intent(inout) :: r 

      call setInvalid(r%p)

      call setINvalid(r%q)

  end subroutine setInvalidRatint

  pure function isZero_flat(r) result(x)
      type(ratint_t), intent(in) :: r 
      logical :: x

      x = .false.

      if (limbiszero(r%p)) then 
          x = .true.
          return 
      end if 

  end function isZero_flat


  pure function isZero_Vector(r) result(x)
      type(ratint_t), dimension(:), intent(in) :: r 
      logical :: x
      integer(shortInt) :: i, height

      x = .true.

      height = size(r)

      do i=1, height
          if (.not. isZero_flat(r(i))) then 
              x = .false. 
              return 
          end if 
      end do 

  end function isZero_Vector



  pure function get_numerator(r) result(n)
      type(ratint_t), intent(in) :: r 
      type(limb_t) :: n 

      n = r%p 
  end function get_numerator

  pure function get_denominator(r) result(d)
      type(ratint_t), intent(in) :: r 
      type(limb_t) :: d

      d = r%q
  end function get_denominator

  pure function evaluate(r) result(v)
    type(ratint_t), intent(in) :: r 
    real(defReal) :: v 

    v = r%p / r%q

  end function evaluate


  elemental type(ratint_t) function addpure(r1, r2)
    type(ratint_t), intent(in) :: r1, r2
    type(ratint_t)             :: r1t, r2t
    type(ratint_t)             :: r3

    if(checkInvalidRatint(r1) .or. checkInvalidRatint(r2)) then 
      call setInvalidRatint(r3)
      addpure = r3
      return

    end if

    ! Modify numerators so that denominators are the same
    ! NOTE: Because of lcm calculation this is guaranteed to be a whole number
    r1t%p = r1%p * r2%q
    r2t%p = r2%p * r1%q

    r3%p = r1t%p + r2t%p

    ! Sets the denominator to be the lowest common multiple
    r3%q = r1%q * r2%q 
    addpure = r3

    ! Normalise signs.
    call normaliseSign(addpure)

  end function addpure

  ! Allows addition between : int + ratint
  elemental type(ratint_t) function addmixedL(n, r1)
    integer(longInt), intent(in) :: n
    type(ratint_t), intent(in) :: r1
    type(ratint_t) :: rn 

    rn = convert_int(n)
    addmixedL = addpure(rn, r1)

  end function addmixedL


  ! Allows addition between : ratint + int
  elemental type(ratint_t) function addmixedR(r1, n)
    integer(longInt), intent(in) :: n
    type(ratint_t), intent(in) :: r1
    type(ratint_t) :: rn 

    rn = convert_int(n)
    addmixedR = addpure(r1, rn)

  end function



  ! negates the second value and adds the results
  elemental type(ratint_t) function subtractpure(r1, r2)
    type(ratint_t), intent(in) :: r1, r2 
    type(ratint_t) :: r2t

    r2t%p = r2%p
    r2t%p%sign = r2%p%sign * (-1)
    r2t%q = r2%q

    subtractpure = addpure(r1, r2t)
  
  end function subtractpure


  ! Allows subtraction between : int - ratint
  elemental type(ratint_t) function subtractmixedL(n, r1)
    integer(longInt), intent(in) :: n
    type(ratint_t), intent(in) :: r1
    type(ratint_t) :: rn 

    rn = convert_int(n)

    subtractmixedL = subtractpure(rn, r1)

  end function subtractmixedL

  
  ! Allows subtraction between : ratint - int
  elemental type(ratint_t) function subtractmixedR(r1, n)
      integer(longInt), intent(in) :: n
      type(ratint_t), intent(in) :: r1
      type(ratint_t) :: rn 

      rn =convert_int(n)

      subtractmixedR = subtractpure(r1, rn)
  end function subtractmixedR

  
  !! multiplies numerator and denominator then simplifies the fraction
  elemental type(ratint_t) function multiplypure(r1, r2)
    type(ratint_t), intent(in) :: r1, r2 
    type(ratint_t) :: r3 

    if (checkInvalidRatint(r1) .or. checkInvalidRatint(r2)) then 
      call setInvalidRatint(r3)
      multiplypure = r3
      return

    end if 

    r3%p = r1%p * r2%p 
    r3%q = r1%q * r2%q

    multiplypure = r3

    ! Normalise signs.
    call normaliseSign(multiplypure)

  end function multiplypure



  ! Allows multiplication between : int * ratint
  elemental type(ratint_t) function multiplymixedL(n, r1)
      integer(longInt), intent(in) :: n
      type(ratint_t), intent(in) :: r1
      type(ratint_t) :: rn 

      rn = convert_int(n)

      multiplymixedL = multiplypure(rn, r1)

  end function multiplymixedL

  

  ! Allows multiplication between : ratint * int
  elemental type(ratint_t) function multiplymixedR(r1, n)
      integer(longInt), intent(in) :: n
      type(ratint_t), intent(in) :: r1
      type(ratint_t) :: rn 

      rn =convert_int(n)

      multiplymixedR = multiplypure(r1, rn)
  end function multiplymixedR


  pure function multiplyVectorL(n, r1) result(rn)
      integer(longInt), intent(in) :: n
      type(ratint_t), dimension(:), intent(in) :: r1
      type(ratint_t), dimension(size(r1)) :: rn 
      integer :: i 

      do i=1, size(r1)
          rn(i) = r1(i) * n 

      end do


  end function multiplyVectorL

  pure function multiplyVectorR(r1, n) result(rn)
      integer(longInt), intent(in) :: n
      type(ratint_t), dimension(:), intent(in) :: r1
      type(ratint_t), dimension(size(r1)) :: rn 
      integer :: i 

      do i=1, size(r1)
          rn(i) = r1(i) * n 

      end do


  end function multiplyVectorR


  

  ! Follows keep, change, flip rule, then applies multiplication
  ! NOTE: division by 0 causes NaN via modulo() call in gcd

  elemental type(ratint_t) function dividepure(r1,r2)
    type(ratint_t), intent(in) :: r1,r2 
    type(ratint_t) :: r3
    type(limb_t) :: temp 

    temp = r2%p 
    r3%p = r2%q 
    r3%q = temp 

    dividepure = multiplypure(r1, r3)

    ! Normalise signs.
    call normaliseSign(dividepure)

  end function dividepure


  ! Allows division between : int / ratint
  elemental type(ratint_t) function dividemixedL(n, r1)
      integer(longInt), intent(in) :: n
      type(ratint_t), intent(in) :: r1
      type(ratint_t) :: rn 

      rn = convert_int(n)

      dividemixedL = dividepure(rn, r1)

  end function dividemixedL

  
  ! Allows division between : ratint / int
  elemental type(ratint_t) function dividemixedR(r1, n)
      integer(longInt), intent(in) :: n
      type(ratint_t), intent(in) :: r1
      type(ratint_t) :: rn 

      rn = convert_int(n)

      dividemixedR = dividepure(r1, rn)
  end function dividemixedR


  ! Copies over the values from rin (r input) into rout (r output)
  pure subroutine assignpure(rout, rin)
      type(ratint_t), intent(out) :: rout 
      type(ratint_t), intent(in) :: rin

      if (checkInvalidRatint(rin)) then 
          call setInvalidRatint(rout)
          return 
      end if

      rout%p = rin%p 
      rout%q = rin%q

  end subroutine assignpure

  
  ! Simplifies the input via the gcd method
  pure function simplify(r) result(rs)
      type(ratint_t), intent(in) ::  r 
      type(ratint_t) :: rs 
      type(limb_t) :: rp 
      type(limb_t)  :: rq 
      type(limb_t) :: gcdval


      rp = r%p 
      rq = r%q

      ! Gets gcd between numerator and denominator (p,q)
      gcdval = gcd(rp, rq)
      ! division by zero check
      if (limbiszero(gcdVal)) then 
          rs%p = initlimb4(0)
          rs%q = initlimb4(0)
      else
          ! Sets the new numerator and denominator 
          rp = initlimb4(floor(rp / gcdVal))
          rq = initlimb4(floor(rq / gcdVal))


          rs%p = rp 
          rs%q = rq
      end if 


  end function simplify





  ! gcd via modulus version as all numerator/denominator are positive
  pure function gcd (a,b) result(v)
      type(limb_t), intent(in) :: a,b
      type(limb_t) :: at, bt
      integer(shortInt) :: temp
      type(limb_t) :: v
      integer(shortInt) :: asign, bsign 

      at = a 
      bt = b

      at%sign = 1
      bt%sign = 1


      do while (.not. (at == bt))

          if (at > bt) then 
              at = at - bt 
          else 
              bt = bt - at
          end if 
      end do 

      v = at 
      
  end function gcd


  pure function equality(a, b) result(r)
      type(ratint_t), intent(in) :: a,b 
      type(limb_t) :: prod1, prod2
      logical :: r 

      r = .false.

      if (checkInvalidRatint(a) .or. checkInvalidRatint(b)) then 
          return 
      end if

      if (a%q == b%q) then 
          if (a%p == b%p) then 
              r = .true. 
              return 
          end if 
      else 
          prod1 = a%p * b%q 
          prod2 = b%p * a%q 
          if (prod1 == prod2) then 
              r = .true.
              return 
          end if 
      end if

  
  end function equality


  pure function gneq(a, b) result(r)
      type(ratint_t), intent(in) :: a,b 
      type(limb_t) :: prod1, prod2
      logical :: r 

      r = .false.

      if (checkInvalidRatint(a) .or. checkInvalidRatint(b)) then 
          return 
      end if


      if (a%q == b%q) then 
          if (a%p > b%p) then 
              r = .true. 
              return 
          end if 
      else 
          prod1 = a%p * b%q 
          prod2 = b%p * a%q 
          if (prod1 > prod2) then 
              r = .true.
              return 
          end if 
      end if
  end function gneq 




  pure function geq(a,b) result(r)
      type(ratint_t), intent(in) :: a,b 
      type(limb_t) :: prod1, prod2
      logical :: r 

      r = .false.

      if (checkInvalidRatint(a) .or. checkInvalidRatint(b)) then 
          return 
      end if


      if (a%q == b%q) then 
          if (a%p >= b%p) then 
              r = .true. 
              return 
          end if 
      else 
          prod1 = a%p * b%q 
          prod2 = b%p * a%q 
          if (prod1 >= prod2) then 
              r = .true.
              return 
          end if 
      end if

  end function geq

  pure subroutine setZero(r)
      type(ratint_t), intent(inout) :: r 
      r%p = initlimb(0_8)
      r%q = initlimb(1_8)
  end subroutine setZero



  subroutine printRatInt(a)
      type(ratint_t), intent(in) :: a 
      print *, 'Numerator'
      call printlimb(a%p)
      print *, '/////'
      print *, 'Denominator'
      call printlimb(a%q)
  end subroutine printRatInt



  pure function dot_product_ratint(a,b) result(r)
      type(ratint_t), intent(in) :: a(:), b(:) 
      type(ratint_t) :: r
      integer(shortInt) :: height, i

      call setZero(r)

      if (size(a) /= size(b)) then 
          call setInvalidRatint(r)
          return 
      end if 

      height = size(a)

      do i = 1, height
          r = r + (a(i) * b(i))
      end do

  end function dot_product_ratint

end module



! program test 
!     use limb_class
!     use ratint
!     use, intrinsic :: iso_fortran_env
!     use, intrinsic :: ieee_arithmetic

!     implicit none 

!     real(defReal) :: v1, v2, v3, eval
!     type(ratint_t) :: ratint1, ratint2, vres
!     logical :: result


!     v1 = 1.0_real64 / 75.0_real64
!     print *, '----'
!     print *, v1

!     ratint1 = convert_ieee64(v1)
!     call printRatInt(ratint1)
!     eval = evaluate(ratint1)
!     print *, 'evaluated:'
!     print *, eval 

    
    
!     result = v1 == eval



! end program test