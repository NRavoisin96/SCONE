module ratint
    
    use limb_class
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    implicit none 

    type ratint_t
        type(limb_t) :: p
        type(limb_t) :: q
    end type ratint_t

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
    end interface operator (*)

    interface operator (/)
        module procedure dividepure
        module procedure dividemixedL
        module procedure dividemixedR
    end interface operator (/)

    interface assignment (=)
        module procedure assignpure
    end interface assignment (=)


    contains 

        function def_ratint(n, d, s) result(r)
            integer(8), intent(in) :: n, d, s
            type(ratint_t) :: r 

            r%p = initlimb(n)
            r%q = initlimb(d)

            !r = simplify(r)

        end function def_ratint


        pure function convert_int(n) result(r)
            integer(8), intent(in) :: n 
            type(ratint_t) :: r 

            r%q = initlimb(1_8) 
            r%p = initlimb(n)

        end function convert_int

        function convert_ieee64(n) result(r)
            real(real64), intent(in) :: n 
            type(ratint_t) :: r 
            real(real64) :: n1
            integer(8) :: i, shift

            real(real64) :: frac, exp
            type(ratint_t) :: expratint, fracratint

            !print *, n
            frac = fraction(n)
            exp = exponent(n)

            ! !print *, frac * (2**exp)
            ! print *, 'in func'
            ! print *, frac 
            ! print *, exp
            shift = 0


            
            do i=0, 52

                if (int(frac, 8)*1_real64 == frac) then 
                    shift = i 
                    exit 
                end if 
                shift = i+1

                frac = frac * 2
            end do
            ! print *, shift
            ! print *, frac
            
            
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

            ! print *, '////'
            ! call printlimb(fracratint%p)
            ! print *, '--'
            ! call printlimb(fracratint%q)
            ! print *, '==='
            ! call printlimb(expratint%p)
            ! print *, '--'
            ! call printlimb(expratint%q)

            r = expratint * fracratint

            if (sign(1.0_real64,n) == -1.0_real64) then 
                r%p%sign = -1
            end if

            ! print *, '-!!!---'
            ! call printRatInt(r)

            ! print *, evaluate(r)


        end function convert_ieee64


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


   
        
        function evaluate(r) result(v)
            type(ratint_t), intent(in) :: r 
            real(8) :: v 

            v = r%p / r%q

        end function evaluate


        type(ratint_t) function addpure(r1, r2)
            type(ratint_t), intent(in) :: r1, r2
            type(ratint_t) :: r1t, r2t
            type(ratint_t) :: r3
            type(limb_t) :: gcdVal
            type(limb_t) :: lcm
  
            ! Get the greatest common divisor and least common multiple
            !gcdVal = gcd(r1%q, r2%q)
            
            !lcm = initlimb(1_8)
            ! NOTE: the result of this division should be an exact value, so its safe to round down
            !lcm = floor(r1%q/gcdVal)*1_8 * r2%q
            ! print *, 'lcm'
            ! call printlimb(lcm)


            ! Modify numerators so that denominators are the same
            ! NOTE: Because of lcm calculation this is guaranteed to be a whole number
            !r1t%p = r1%p * initlimb(int(lcm/r1%q , 8))
            r1t%p = r1%p * r2%q
            ! print *, '----'
            ! call printRatInt(r1)
            ! call printRatInt(r2)
            ! call printlimb(r1t%p)
            ! call printlimb(r2t%p)
            !print *, 'lkfdkf'
            !call printlimb(r1t%p)
            !r2t%p = r2%p * initlimb(int(lcm/r2%q , 8))
            r2t%p = r2%p * r1%q
            ! call printlimb(r2t%p)


            r3%p = r1t%p + r2t%p

            ! Sets the denominator to be the lowest common multiple
            r3%q = r1%q * r2%q 

            !addpure = simplify(r3)


            addpure = r3


        end function addpure

        ! Allows addition between : int + ratint
        type(ratint_t) function addmixedL(n, r1)
            integer(8), intent(in) :: n
            type(ratint_t), intent(in) :: r1
            type(ratint_t) :: rn 

            rn = convert_int(n)

            addmixedL = addpure(rn, r1)

        end function addmixedL


        ! Allows addition between : ratint + int
        type(ratint_t) function addmixedR(r1, n)
            integer(8), intent(in) :: n
            type(ratint_t), intent(in) :: r1
            type(ratint_t) :: rn 

            rn = convert_int(n)

            addmixedR = addpure(r1, rn)
        end function



        ! negates the second value and adds the results
        type(ratint_t) function subtractpure(r1, r2)
            type(ratint_t), intent(in) :: r1, r2 
            type(ratint_t) :: r2t

            r2t%p = r2%p
            r2t%p%sign = r2%p%sign * (-1)
            r2t%q = r2%q

            subtractpure = addpure(r1, r2t)
        
        end function subtractpure


        ! Allows subtraction between : int - ratint
        type(ratint_t) function subtractmixedL(n, r1)
            integer(8), intent(in) :: n
            type(ratint_t), intent(in) :: r1
            type(ratint_t) :: rn 

            rn = convert_int(n)

            subtractmixedL = subtractpure(rn, r1)

        end function subtractmixedL

        
        ! Allows subtraction between : ratint - int
        type(ratint_t) function subtractmixedR(r1, n)
            integer(8), intent(in) :: n
            type(ratint_t), intent(in) :: r1
            type(ratint_t) :: rn 

            rn =convert_int(n)

            subtractmixedR = subtractpure(r1, rn)
        end function subtractmixedR

        
        !! multiplies numerator and denominator then simplifies the fraction
        type(ratint_t) function multiplypure(r1, r2)
            type(ratint_t), intent(in) :: r1, r2 
            type(ratint_t) :: r3 

            r3%p = r1%p * r2%p 
            r3%q = r1%q * r2%q

            multiplypure = r3
            !multiplypure = simplify(r3)

        end function multiplypure



        ! Allows multiplication between : int * ratint
        type(ratint_t) function multiplymixedL(n, r1)
            integer(8), intent(in) :: n
            type(ratint_t), intent(in) :: r1
            type(ratint_t) :: rn 

            rn = convert_int(n)

            multiplymixedL = multiplypure(rn, r1)

        end function multiplymixedL

        

        ! Allows multiplication between : ratint * int
        type(ratint_t) function multiplymixedR(r1, n)
            integer(8), intent(in) :: n
            type(ratint_t), intent(in) :: r1
            type(ratint_t) :: rn 

            rn =convert_int(n)

            multiplymixedR = multiplypure(r1, rn)
        end function multiplymixedR


        

        ! Follows keep, change, flip rule, then applies multiplication
        ! NOTE: division by 0 causes NaN via modulo() call in gcd

        type(ratint_t) function dividepure(r1,r2)
            type(ratint_t), intent(in) :: r1,r2 
            type(ratint_t) :: r3
            type(limb_t) :: temp 

            temp = r2%p 
            r3%p = r2%q 
            r3%q = temp 

            dividepure = multiplypure(r1, r3)

        end function dividepure


        ! Allows division between : int / ratint
        type(ratint_t) function dividemixedL(n, r1)
            integer(8), intent(in) :: n
            type(ratint_t), intent(in) :: r1
            type(ratint_t) :: rn 

            rn = convert_int(n)

            dividemixedL = dividepure(rn, r1)

        end function dividemixedL

        
        ! Allows division between : ratint / int
        type(ratint_t) function dividemixedR(r1, n)
            integer(8), intent(in) :: n
            type(ratint_t), intent(in) :: r1
            type(ratint_t) :: rn 

            rn = convert_int(n)

            dividemixedR = dividepure(r1, rn)
        end function dividemixedR


        ! Copies over the values from rin (r input) into rout (r output)
        subroutine assignpure(rout, rin)
            type(ratint_t), intent(out) :: rout 
            type(ratint_t), intent(in) :: rin

            rout%p = rin%p 
            rout%q = rin%q

        end subroutine assignpure

        
        ! Simplifies the input via the gcd method
        function simplify(r) result(rs)
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
        function gcd (a,b) result(v)
            type(limb_t), intent(in) :: a,b
            type(limb_t) :: at, bt
            integer :: temp
            type(limb_t) :: v
            integer :: asign, bsign 

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


        subroutine printRatInt(a)
            type(ratint_t), intent(in) :: a 
            print *, 'Numerator'
            call printlimb(a%p)
            print *, '/////'
            print *, 'Denominator'
            call printlimb(a%q)
        end subroutine printRatInt


        ! Potentially implement, but apparently speed difference isn't too big
        !pure function binarygcd(a,b) result(v)
        !    integer, intent(in) :: a,b 
        !    integer :: at, bt
        !    integer :: v

        
        !end function binarygcd


end module



! program test 
!     use limb_class
!     use ratint
!     use, intrinsic :: iso_fortran_env
!     use, intrinsic :: ieee_arithmetic

!     implicit none 

!     real(real64) :: v1, v2, v3, eval
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