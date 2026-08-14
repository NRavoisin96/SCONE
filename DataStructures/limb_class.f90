module limb_class
    
    use numPrecision

    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic

    implicit none
    
    integer, parameter :: maxsz = 1977

    type limb_t
    ! Indexed from 0 to allow the lowest index to act as a 'null' index for managing unfilled values
    ! Only use 1-1977 in calculations
      integer(shortInt)                     :: front = -1, sign = 1
      integer(shortInt), dimension(0:maxsz) :: limbs
    end type limb_t

    interface initlimb 
        module procedure initlimbempty, initlimbnumlong!, initlimbnumshort
    end interface initlimb

    interface operator (+)
        module procedure addlimbs
        module procedure addlimbmixedL
        module procedure addlimbmixedR
    end interface operator (+)

    interface operator (-)
        module procedure subtractlimbs
        module procedure subtractlimbmixedL
        module procedure subtractlimbmixedR
    end interface operator (-)

    interface operator (*)
        module procedure multiplylimbs
        module procedure multiplylimbmixedL
        module procedure multiplylimbmixedR
    end interface operator (*)

    interface operator (/)
        module procedure dividelimbs
    !    module procedure dividemixedL
        module procedure dividelimbsmixedR
    end interface operator (/)

    interface operator (==)
        module procedure limbequality
    end interface operator (==)

    interface operator (>)
        module procedure limbgneq
    end interface operator (>)

    interface operator (>=)
        module procedure limbgeq
    end interface operator (>=)

    interface assignment (=)
        module procedure assignlimbs
    end interface assignment (=)

    contains

        pure function initlimbempty() result(l)
            type(limb_t) :: l 
            integer(4) :: i

            l%front = 1
            l%sign = 1
            l%limbs(1) = 0
            l%limbs(0) = 0

            !do i=1, maxsz 
            !    l%limbs(i) = 0
            !end do

        end function initlimbempty


        pure function initlimbFull() result(l)
            type(limb_t) :: l 
            integer(4) :: i

            l%front = maxsz
            l%sign = 1
            l%limbs(1) = 0
            l%limbs(0) = 0

            do i=1, maxsz 
               l%limbs(i) = 0
            end do

        end function initlimbFull


        pure function initlimbSize(n) result(l)
            integer, intent(in) :: n 
            type(limb_t) :: l 
            integer(4) :: i

            if (n > 1977) then 
                call setInvalid(l)
                return
            end if

            l%front = n
            l%sign = 1
            l%limbs(0) = 0

            do i=1, n 
               l%limbs(i) = 0
            end do

        end function initlimbSize


!! CREATE A GENERAL FUNCTION TO ALLOW DISTINGUISHING BETWEEN SHORT/LONG, CURRENT JUST DEFAULTS TO LONG
        pure function initlimbnumlong(n) result(l)
            integer(longInt), intent(in) :: n
            type(limb_t) :: l 
            integer(longInt) :: st
            integer(longInt) :: overflow
            integer(4) :: i


            l = initlimbempty()
            !l%sign = 1
            l%sign = sign(int(1, 8), n)
            !do i=1, maxsz 
            !    l%limbs(i) = 0
            !end do

            overflow = abs(n)
            st = 0

            ! Do loop to deal with rolling carries, exits if none
            do i=1, maxsz
                st = overflow
                overflow = shiftl(shiftr(abs(st), 31), 31)
                l%limbs(i) = st - overflow

                overflow = shiftr(overflow, 31)
                
                if (overflow == 0) then 
                    l%front = i
                    exit
                end if 

            end do

    
        end function initlimbnumlong


        pure function initlimb4(n) result(l)
            integer(4), intent(in) :: n
            type(limb_t) :: l 
            integer(4) :: st
            integer(4) :: overflow
            integer(4) :: i


            l = initlimbempty()
            !l%sign = 1
            l%sign = sign(1_4, n)
            overflow = 0
            !do i=1, maxsz 
            !    l%limbs(i) = 0
            !end do

            overflow = abs(n)
            st = 0

            ! Do loop to deal with rolling carries, exits if none
            do i=1, maxsz
                st = overflow
                overflow = shiftl(shiftr(abs(st), 31), 31)
                l%limbs(i) = st - overflow

                overflow = shiftr(overflow, 31)
                
                if (overflow == 0) then 
                    l%front = i
                    exit
                end if 

            end do


        end function initlimb4


        pure subroutine setInvalid(a)
            type(limb_t), intent(inout) :: a 

            a%front = -1 
            a%sign = 0 

        end subroutine setInvalid


        pure function checkInvalid(a) result (r)
            type(limb_t), intent(in) :: a 
            logical :: r 

            r = a%front == -1 

        end function checkInvalid
        


        pure subroutine addOnSign(a,b, asig, bsig, s)
            type(limb_t), intent(in)      :: a,b 
            integer(shortInt), intent(in) :: asig, bsig
            type(limb_t), intent(inout)   :: s
            integer(longInt)              :: st, zeroIndA, zeroIndB, overflow, borrow, va, vb
            integer(shortInt)             :: i, maxfront

            if (a%front == 1977 .and. b%front == 1977 .and. (a%limbs(a%front)*1_8+b%limbs(b%front)*1_8) > 2_8**31) then 
                call setInvalid(s)
                return 
            end if

            st = 0

            maxfront = max(a%front, b%front)
            overflow = 0
            borrow = 0


            zeroIndA = 0
            zeroIndB = 0

            do i=1, maxsz
                if (i > maxfront .and. overflow == 0) then 
                    s%front = s%front - (1 - min(1_8, st))
                    exit
                end if



                ! ! Indicator for preventing the use of an index outside of limb size
                ! zeroIndA = max(min(a%front+1-i, 1), 0)
                ! zeroIndB = max(min(b%front+1-i, 1), 0)

                ! ! Add components based on sign
                ! st = (a%limbs(i*zeroInda))*asig + (b%limbs(i*zeroIndb))*bsig + overflow + borrow
                ! borrow = 0

                if (i > a%front) then 
                    va = 0
                else 
                    va = a%limbs(i) * asig 
                end if 

                if (i > b%front) then 
                    vb = 0
                else 
                    vb = b%limbs(i)*bsig 
                end if 
                
                st = va + vb + overflow + borrow
                borrow = 0
                

                ! Borrows from above to make the value positive, stores in borrow
                do while (st < 0) 
                    borrow = borrow - 1
                    st = (2_8)**31 + st 
                end do 

                ! Compute and remove overflow via 31-bit shifting
                overflow = shiftl(shiftr(st, 31), 31) 

                s%limbs(i) = st - overflow 
                ! Any overflow is by at most one bit
                overflow = shiftr(overflow, 31)
                
                s%front = i
                
                !if (i>=maxfront .and. overflow == 0) then 
                    ! Covers the case of if the final sum is 0 (prevents leading 0)
                    ! NOTE: If this doesn't have a max(..., 1) it WILL cause a zero index issue if the sum is 0
                !    s%front = max(int(i - (1 - min(1_8, st)), 4), 1)
                 !   exit
                !end if
                
            end do 


            call adjustFront(s)

        end subroutine addOnSign


        pure function addLimbs(a,b) result (s)
            type(limb_t), intent(in) :: a, b 
            type(limb_t) :: s 
            logical :: greaterThan, tr1, tr2, tr3, tr4, tr5, tr6
            integer(shortInt) cd1, cd2, cd3, cd4, cd5, cd6, greaterThanInt

            if (checkInvalid(a) .or. checkInvalid(b)) then 
                call setInvalid(s)
                return 
            end if

            ! Initialising s 
            s = initlimb()

            

            ! Alternative representation using max, min and integers to replace booleans 
            ! Less readable but im not sure if merge triggers branch prediction
            ! cds are equivalent to trs above
            ! greaterThanInt = absgneqInt(a,b)
            ! cd1 = 1 - max(min(a%sign*b%sign, 1_8), 0_8)
            ! cd2 = max(min(a%sign*(-1_8), 1_8), 0_8) * greaterThanInt
            ! cd3 = (1 - greaterThanInt) * max(min(b%sign*(-1_8), 1_8), 0_8)
            ! cd4 = max(min(a%sign*b%sign, 1_8), 0_8) * max(min(a%sign*(-1_8), 1_8), 0_8)
            ! cd5 = max(min(a%sign*(-1_8), 1_8), 0_8) * (1 - greaterThanInt) 
            ! cd6 = max(min(b%sign*(-1_8), 1_8), 0_8) * greaterThanInt

            ! s%sign = 1_8 + (-2_8) * (min(cd4 + (cd1 * (min(cd2 + cd3, 1))), 1))

            ! call addOnSign(a, &
            !                 b, &
            !                 1_8 + (-2_8)*(min((cd1 * cd3) + (cd1 * cd5), 1)), &
            !                 1_8 + (-2_8)*(min((cd1 * cd2) + (cd1 * cd6), 1)), s)

            ! Cases (for setting sign of result):
            ! (1) If signs are different and magnitudes mean their sum is negative, sets sign = -1
            ! (2) If both signs are the same and negative, sets sign = -1
            ! (3) Any other case only has positive signs, sets sign = 1

            ! Cases (for sign inputs):
            ! (1) If a,b are the same sign, both sides of 'or' fail, sets input sign to 1 (  sum)
            ! (2) If a,b signs different:
                ! (i)  a = negative with larger magnitude, tr3 and tr5 fail so sign for a = 1
                    ! because of above, tr6 fails but tr2 succeeds, so sign for b = -1
                    ! (subtracting smaller magnitude from bigger one)
                ! (ii) b = negative with larger magnitude, tr3 succeeds, so sign for a = -1
                    ! because of above, tr2 and tr6 fail, so sign for b = 1
                ! (iii) greater magitude is positive, smaller is negative, so signs are kept the same
                    ! if a is negative, tr5 succeeds and tr6 fails (signs remain the same)
                    ! is b is negative, tr6 succeeds and tr5 fails

            ! greaterthan = absgneq(a, b)

            ! tr1 = (a%sign /= b%sign)  
            ! tr2 = (greaterthan .and. a%sign == -1)
            ! tr3 = ((.not. greaterthan) .and. b%sign == -1)
            ! tr4 = (a%sign == b%sign) .and. a%sign == -1
            ! tr5 = (a%sign == -1) .and. (.not. greaterThan)
            ! tr6 = (b%sign == -1) .and. (greaterThan)

            ! s%sign = 1_8 + (-2_8) * merge(1, 0, (tr1 .and. (tr2 .or. tr3)) .or. tr4)

            ! call addOnSign(a, &
            !                b, &
            !                1_8 + (-2_8)*merge(1, 0, (tr1 .and. tr3) .or. (tr1 .and. tr5)), &
            !                1_8 + (-2_8)*merge(1, 0, (tr1 .and. tr2) .or. (tr1 .and. tr6)), s)

            

            if (a%sign /= b%sign) then 
                greaterthan = absgneq(a, b)

                if (greaterthan .and. a%sign == -1) then 
                    s%sign = -1
                    call addOnSign(a, b, 1, -1, s)
                    return

                else if ((.not. greaterthan) .and. b%sign == -1) then 
                    s%sign = -1
                    call addOnSign(a, b, -1, 1, s)
                    return

                end if
            else 
                s%sign = a%sign
                call addOnSign(a, b, 1, 1, s)
                return
                
            end if
            
            s%sign = 1
            call addOnSign(a, b, a%sign, b%sign, s)

        end function addLimbs


        pure function absgneq(a,b) result (r)
            type(limb_t), intent(in) :: a,b
            logical :: r 
            integer(shortInt) i

        
            r = .false.

            if (checkInvalid(a) .or. checkInvalid(b)) then 
                return 
            end if
   
            ! Catches any potential zero index front error (remove?? may mask issues)
            if (a%front /= b%front) then
                r = (a%front*min(1, a%limbs(a%front)) > b%front*min(1, b%limbs(b%front)))
                return 
            end if

            do i=a%front, 1, -1 
                if (a%limbs(i) > b%limbs(i)) then 
                    r = .true. 
                    exit
                else if (a%limbs(i) < b%limbs(i)) then
                    exit 
                end if 
            end do 
        end function absgneq



        pure function absgneqInt(a,b) result (r)
            type(limb_t), intent(in) :: a,b
            integer(shortInt) r 
            integer(shortInt) i

            
            r = 0

            if (checkInvalid(a) .or. checkInvalid(b)) then 
                return 
            end if
   
            ! min(1, ...) catches any potential zero index front error (remove?? may mask issues)
            if (a%front /= b%front) then
                r = max(min(a%front*min(1, a%limbs(a%front))-b%front*min(1, b%limbs(b%front)),1) , 0)
                !r = (a%front*min(1, a%limbs(a%front)) > b%front*min(1, b%limbs(b%front)))
                return 
            end if

            do i=a%front, 1, -1 
                if (a%limbs(i) > b%limbs(i)) then 
                    r = 1 
                    exit
                else if (a%limbs(i) < b%limbs(i)) then
                    exit 
                end if 
            end do 
        end function absgneqInt

        

        pure function addlimbmixedL(n,a) result(s)
            integer(longInt), intent(in) :: n
            type(limb_t), intent(in) :: a 
            type(limb_t) :: s, nlimb

            nlimb = initlimb(n)

            s = a + nlimb

        end function addlimbmixedL


        pure function addlimbmixedR(a,n) result(s)
            integer(longInt), intent(in) :: n
            type(limb_t), intent(in) :: a 
            type(limb_t) :: s, nlimb


            nlimb = initlimb(n)

            s = a + nlimb

        end function addlimbmixedR





        pure function subtractlimbs(a, b) result(s)
            type(limb_t), intent(in) :: a,b 
            type(limb_t) :: s, bt

            bt = b 
            bt%sign = (-1)*bt%sign 

            s = a + bt

        end function subtractlimbs



        pure function subtractlimbmixedL(n,a) result(s)
            integer(longInt), intent(in) :: n
            type(limb_t), intent(in) :: a 
            type(limb_t) :: s, nlimb 

            nlimb = initlimb(n)

            s = nlimb - a 

        end function subtractlimbmixedL 



        pure function subtractlimbmixedR(a,n) result(s)
            integer(longInt), intent(in) :: n
            type(limb_t), intent(in) :: a 
            type(limb_t) :: s, nlimb 

            nlimb = initlimb(n)

            s = a - nlimb

        end function subtractlimbmixedR



        pure function multiplylimbs(a,b) result(p)
            type(limb_t), intent(in) :: a,b 
            type(limb_t) :: p

            if (checkInvalid(a) .or. checkInvalid(b)) then 
                call setInvalid(p)
                return 
            end if

            p = initlimb()

            if (limbiszero(a) .or. limbiszero(b)) then 
                return 
            end if

            call recursiveMultLimbs(a, b, 1, a%front, 1, b%front, p)

            call adjustFront(p)


            p%sign = a%sign * b%sign
        end function multiplylimbs


        pure subroutine recursiveMultLimbs(a, b, asplitstart, asplitend, bsplitstart, bsplitend, p)
            type(limb_t), intent(in) :: a,b 
            integer, intent(in) :: asplitstart, bsplitstart, asplitend, bsplitend
            type(limb_t) , intent(inout) :: p

            if (checkInvalid(a) .or. checkInvalid(b)) then 
                call setInvalid(p)
                return 
            end if

            if ((asplitend-asplitstart+1) < 32 .or. (bsplitend-bsplitstart+1) < 32) then 
                call longMultiplication(a, b, asplitstart, asplitend, bsplitstart, bsplitend, p)

            else 
                call karatsubamult(a, b, asplitstart, asplitend, bsplitstart, bsplitend, p)
            end if

        end subroutine recursiveMultLimbs

        ! a * b, i.e.
        !   (a_m a_m-1 .... a_2 a_1)
        ! * (b_n b_n-1 .... b_2 b_1)
        pure subroutine longMultiplication(a, b, aStart, aEnd, bStart, bEnd, p)
            type(limb_t), intent(in) :: a, b 
            integer, intent(in) :: aStart, aEnd, bStart, bEnd 
            type(limb_t), intent(inout) :: p 
            integer(shortInt) i, j, pos, zeroIndA 
            integer(longInt) :: pt, overflow, va

            ! Placing the check here, however, n+m is only an estimation, ideally it should be in the loop,
            ! but this will affect performance
            ! It seems like past numbers beginning with 316... n+m applies (verify later)

            if (a%front + b%front > 1977 .or. (checkInvalid(a) .or. checkInvalid(b))) then 
                call setInvalid(p)
                return
            end if


            p = initLimbSize((aEnd-aStart+1)+(bEnd-bStart+1)+1)
            pos = 0 


            do i=bStart, bEnd
                pos = i - bStart +1
                if (b%limbs(i) == 0) then 
                    cycle 
                end if

                overflow = 0
                pt = 0


                do j=aStart, maxsz
                    if (j > aEnd .and. overflow == 0) then 
                        exit 
                    end if
                    ! zeroIndA = max(min(aEnd - j+1, 1), 0)

                    ! pt = p%limbs(pos) + (b%limbs(i)*1_8 * a%limbs(j*zeroIndA)*1_8) + overflow

                    if (j > aEnd) then 
                        va = 0 
                    else 
                        va = a%limbs(j) 
                    end if 
                    pt = p%limbs(pos) + (b%limbs(i)*1_8 * va* 1_8) + overflow

                    overflow = shiftl(shiftr(pt, 31), 31) 
                    p%limbs(pos) = pt - overflow 

                    ! Any overflow is by at most one bit
                    overflow = shiftr(overflow, 31)

                    pos = pos + 1

                    
                end do 
            end do


            call adjustFront(p)

        end subroutine longMultiplication



        ! NOTE: removed the p limb initialisation at the start of the subroutine, return if issues with accesss
        recursive pure subroutine karatsubamult(a,b, asplitstart, asplitend, bsplitstart, bsplitend, p)
            type(limb_t), intent(in) :: a,b 
            integer, intent(in) :: asplitstart, bsplitstart, asplitend, bsplitend
            type(limb_t) , intent(inout) :: p
            type(limb_t) :: p0, p1, p2, p3, a0a1, b0b1
            integer(longInt) :: z0, z1, z2, z3, a0, a1, b0, b1
            integer(shortInt) shiftmult, splita, splitb, split

            p0 = initlimb()
            p1 = initlimb()
            p2 = initlimb()
            p3 = initlimb()
            a0a1 = initlimb()
            b0b1 = initlimb()



            ! ! Base case
            ! ! Assumes each are at only one split (so can easily be split and multiplied)
            ! if (abs(asplitstart - asplitend) <= 1 .and. abs(bsplitstart - bsplitend) <= 1) then 

            !     ! If either one is at only size 1, then the LSB becomes start and MSB is 0
            !     ! NOTE: relies on end>start, this should hold (if not, signals deeper bug)
            !     a0 = a%limbs(asplitstart)
            !     a1 = a%limbs(asplitend*max(min(asplitend-asplitstart, 1), 0))
            !     b0 = b%limbs(bsplitstart)
            !     b1 = b%limbs(bsplitend*max(min(bsplitend-bsplitstart, 1), 0))


            !     ! Polynomial expansion steps
            !     z0 = a0 * b0
            !     z2 = a1 * b1

            !     ! Karatsuba simplification
            !     z3 = (a0 + a1) * (b0 + b1)
            !     z1 = (z3 - z2) - z0


            !     ! Adding these values to the limbs so that they can be shifted and combined into result p
            !     ! equivalent to: p = z2*2^(31*2) + z1*2^(31) + z0
            !     p0 = initlimb(z0)
            !     p1 = initlimb(z1)
            !     p2 = initlimb(z2)

            !     call tripleShiftSumBy(p0, p1, p2, 1, p)


            !     ! Adjusts the front of the result
            !     call adjustFront(p)


            ! !Recursive step
            ! else
                if ((asplitend-asplitstart+1)+(bsplitend-bsplitstart+1) > 1977 .or. (checkInvalid(a) .or. checkInvalid(b))) then 
                    call setInvalid(p)
                    return
                end if
                ! Calculate the splits to take
                ! NOTE : this is split LENGTH, not relative to the starting point

                splita = (asplitend-asplitstart+1)/2
                splitb = (bsplitend-bsplitstart+1)/2
                split = min(splita, splitb)



                ! Multiplies fronts and backs
                call recursiveMultLimbs(a, b, asplitstart, min(asplitend, (asplitstart + split)-1), & 
                                         bsplitstart, min((bsplitstart + split)-1, bsplitend), p0)
                

               

                call recursiveMultLimbs(a, b, min(asplitend, (asplitstart + split)), asplitend, & 
                                         min(bsplitend, (bsplitstart + split)), bsplitend, p2)


        

                call combineIntervals(a, asplitstart, split, asplitstart+split, (asplitend-asplitstart-split+1), a0a1)


                call combineIntervals(b, bsplitstart, split, bsplitstart+split, (bsplitend-bsplitstart-split+1), b0b1)


                call recursiveMultLimbs(a0a1, b0b1, 1, a0a1%front, 1, b0b1%front, p3)

                p1 = (p3 - p2) - p0


                ! Shift all the values according to the split taken in this recursive step
                ! equivalent to: p = p2*B^2 + p1*B + p0 where B depends on split
                shiftmult = split


                call tripleShiftSumBy(p0, p1, p2, shiftmult, p)
  
                call adjustFront(p)



            ! end if 

        end subroutine karatsubamult

        pure subroutine adjustFront(a)
            type(limb_t), intent(inout) :: a 
            integer(shortInt) i

            do i=a%front,1,-1
                if (a%limbs(i) /= 0) then 
                    a%front = i 
                    return 
                end if 
            end do 

            a%front = 1

        end subroutine

        ! NOTE: Assumes that at,bt are already initialised
        pure subroutine padFronts(a, b, at, bt)
            type(limb_t), intent(in) :: a,b 
            type(limb_t), intent(inout) :: at, bt
            integer(shortInt) i, groupfront, zeroIndA, zeroIndB, va, vb


            groupfront = max(a%front, b%front)

            do i=1, groupfront 
                ! zeroIndA = max(min(a%front+1-i, 1), 0)
                ! zeroIndB = max(min(b%front+1-i, 1), 0)

                ! at%limbs(i) = a%limbs(i*zeroIndA) 
                ! bt%limbs(i) = b%limbs(i*zeroIndB) 
                if (i > a%front) then 
                    va = 0 
                else 
                    va = a%limbs(i)
                end if 

                if (i > b%front) then 
                    vb = 0 
                else 
                    vb = a%limbs(i)
                end if

                at%limbs(i) = va 
                at%limbs(i) = vb
            end do

            at%front = groupfront
            bt%front = groupfront


        end subroutine padFronts



        pure subroutine combineIntervals(a, start1, len1, start2, len2, s)
            type(limb_t), intent(in) :: a 
            integer, intent(in) :: start1, len1, start2, len2
            type(limb_t), intent(inout) :: s 
            integer(shortInt) zeroIndS, zeroIndE, i, zeroIndAS, zeroIndAE
            integer(longInt) :: overflow, st, maxLen, v1, v2

            if (checkInvalid(a)) then 
                call setInvalid(s)
                return 
            end if

            s = initlimb()
            s%sign = a%sign


            overflow = 0
            st = 0
            maxLen = max(len1, len2)


            do i=1, maxsz
                if (i > maxLen) then 
                    s%limbs(i) = overflow
                    s%front = i - (1 - min(overflow, 1_8))
                    exit
                end if
                ! ! Indicator for preventing the use of an index outside of the length to add in the same limb
                ! zeroIndS = max(min(len1-i+1, 1), 0)     
                ! zeroIndE = max(min(len2-i+1, 1), 0)
                ! ! Indicator for preventing access outside of the limb
                ! zeroIndAS = max(min((a%front-start1+1)-i+1, 1), 0)     
                ! zeroIndAE = max(min((a%front-start2+1)-i+1, 1), 0)     
                
                ! ! Add components based on sign
                ! st = (a%limbs((start1+i-1)*zeroIndS*zeroIndAS))*1_8 + (a%limbs((start2 + i-1)*zeroIndE*zeroIndAE))*1_8 + overflow 

                if (i > len1) then 
                    v1 = 0
                else if (i > a%front-start1+1) then 
                    v1 = 0
                else 
                    v1 = a%limbs(start1+i-1)
                end if 

                if (i > len2) then 
                    v2 = 0
                else if (i > a%front-start2+1) then 
                    v2 = 0
                else 
                    v2 = a%limbs(start2+i-1)
                end if 

                st = v1*1_8 + v2*1_8 + overflow 




                ! Compute and remove overflow via 31-bit shifting
                overflow = shiftl(shiftr(st, 31), 31) 
                s%limbs(i) = st - overflow

                ! Any overflow is by at most one bit
                overflow = shiftr(overflow, 31)
                s%front = i
                
            end do

            call adjustFront(s)


        end subroutine combineIntervals



        pure subroutine tripleShiftSumBy(v1, v2, v3, shift, r)
            type(limb_t), intent(in) :: v1, v2, v3 
            integer, intent(in) :: shift
            type(limb_t), intent(inout) :: r
            integer(shortInt) i, maxLen, zeroInd1, zeroInd2, zeroInd3, zeroIndShift1, zeroIndShift2
            integer(longInt) :: overflow, st, borrow, val1, val2, val3

            if (checkInvalid(v1) .or. checkInvalid(v2) .or. checkInvalid(v3)) then 
                call setInvalid(r)
                return 
            end if

            r = initlimb()

            maxLen = max(v3%front + shift*2, v2%front+shift, v1%front)

            overflow = 0
            st = 0
            borrow = 0


            do i = 1, maxsz 
                if (i > maxLen .and. overflow == 0) then 
                    exit 
                end if 
                ! ! Indicators for whether or not the i value (shifted by shift) is outside of limb size
                ! zeroInd1 = max(min(v1%front+1-i, 1), 0)
                ! zeroInd2 = max(min(v2%front+shift+1-i, 1), 0)
                ! zeroInd3 = max(min(v3%front+shift*2+1-i, 1), 0)
                ! ! Indicators for whether or not the value of i (based on shift) can be included in the sum yet
                ! zeroIndShift1 = max(min(i-shift, 1), 0)
                ! zeroIndShift2 = max(min(i - (shift*2), 1), 0)

                ! ! Adds components based on sign
                ! ! Shifted based on shifts and multiplied by indicators
                ! ! If any index is 0, value is 0 (from initialisation)
                ! st = (v1%limbs(i*zeroInd1))*1_8 + (v2%limbs((i-shift)*zeroInd2*zeroIndShift1))*1_8 & 
                !     + (v3%limbs((i-shift*2)*zeroInd3*zeroIndShift2))*1_8 + overflow 


                if (i > v1%front) then 
                    val1 = 0 
                else 
                    val1 = v1%limbs(i)
                end if 

                if (i > v2%front+shift) then 
                    val2 = 0
                else if (i > shift) then 
                    val2 = v2%limbs(i-shift)
                else 
                    val2 = 0
                end if 

                if (i > v3%front+shift*2) then 
                    val3 = 0
                else if (i > shift*2) then 
                    val3 = v3%limbs(i-shift*2)
                else 
                    val3 = 0
                end if 

                st = val1*1_8 + val2*1_8 + val3*1_8 + overflow 


                ! Compute and remove overflow via 31-bit shifting
                overflow = shiftl(shiftr(st, 31), 31) 

                r%limbs(i) = st - overflow 
                ! Any overflow is by at most one bit
                overflow = shiftr(overflow, 31)

                r%front = i

            end do

            call adjustFront(r)


        end subroutine tripleShiftSumBy


        ! Potentially add multiple functions for other types of integer (currently only integer 4 allowed to prevent overflow)
        pure function multiplylimbmixedL(n,a) result(p)
            type(limb_t), intent(in) :: a 
            integer(longInt), intent(in) :: n 
            type(limb_t) :: p
            
            p = a * n


        end function multiplylimbmixedL


        pure function multiplylimbmixedR(a,n) result(p)
            type(limb_t), intent(in) :: a 
            integer(longInt), intent(in) :: n 
            type(limb_t) :: p, nl 

            nl = initlimb(n)

            call longMultiplication(a, nl, 1, a%front, 1, nl%front, p)

            p%sign = a%sign * nl%sign


        end function multiplylimbmixedR


        pure function dividelimbs(a,b) result(x)
            type(limb_t), intent(in) :: a,b 
            real(defReal) :: x, normA, normB
            integer(shortInt) onePosA, onePosB, expA, expB, expX


            call findFirst1(a, onePosA)
            call findFirst1(b, onePosB)

            call tosubnormalisedreal(a, normA)
            call tosubnormalisedreal(b, normB)

            expA = onePosA + (a%front-1)*31 + 1
            expB = onePosB + (b%front-1)*31 + 1


            expX = expA - expB 

            if (expX > 1023) then 
                x = 0 
                return 
            end if

            x = (normA/normB) * (2.0_real64 **expX) * (a%sign*b%sign*1.0_real64)


        end function dividelimbs




        pure function dividelimbs2(a,b) result(x)
            type(limb_t), intent(in) :: a,b 
            type(limb_t) :: at, bt
            real(defReal) :: ra, rb, x 
            real(defReal) :: tr8 
            integer(shortInt) tempfront 
            integer(shortInt) first1loc
            integer(shortInt) shift
            integer(shortInt) i
            integer(shortInt) sign

            if (checkInvalid(a) .or. checkInvalid(b) .or. limbiszero(b)) then 
                x = 0 
                return 
            end if


            ! First finds the position of the first one and records number of shifts made
            tr8 = 1

            sign = a%sign * b%sign 


            call tonormalisedreal(b, rb)
            rb = rb/(2_real64)
            


            tempfront = b%limbs(b%front)
            first1loc = 0
            


            call findFirst1(b, first1loc)

            

            shift = first1loc + (b%front - 1)*31

            !NOTE: issue with the numerator normalisation, likely with the floor() way of getting the front
            ! This was ceiling before, then was changed and it worked until it didnt (on 2654.0_real64 / 9988445522.0_real64)
            call tonormalisedrealby(a, shift, ra)

            ra = ra/(2_real64)




            x = (48.0_real64/17.0_real64) - (32.0_real64/17.0_real64) * rb 
            do i = 1, 6
                x = x + x * (1.0_real64 - (rb * x))
            end do

            x = ra * x * (sign*1_real64)



        end function dividelimbs2


        pure subroutine tonormalisedreal(a, r)
            type(limb_t), intent(in) :: a 
            real(defReal), intent(inout) :: r 
            real(defReal) :: tr8 

            integer(shortInt) tempfront 
            integer(shortInt) first1loc
            integer(shortInt) limit

            ! First finds the position of the first 1 and records number of shifts made
            tempfront = a%limbs(a%front)
            first1loc = 0
            tr8 = 1_real64
            r = 0
            
            call findFirst1(a, first1loc)

            tempfront = a%limbs(a%front)

            !limit = min(51, 31*(a%front -1) + first1loc)
            !limit = 31*(a%front -1) + first1loc

            call sumFractionalComponent(a, a%front, 31*(a%front -1) + first1loc, first1loc, 1, r)
            r = r + 1
        

        end subroutine tonormalisedreal



         pure subroutine tosubnormalisedreal(a, r)
            type(limb_t), intent(in) :: a 
            real(defReal), intent(inout) :: r 
            real(defReal) :: tr8 

            integer(shortInt) tempfront 
            integer(shortInt) first1loc
            integer(shortInt) limit

            ! First finds the position of the first 1 and records number of shifts made
            tempfront = a%limbs(a%front)
            first1loc = 0
            tr8 = 1_real64
            r = 0
            
            call findFirst1(a, first1loc)

            tempfront = a%limbs(a%front)

            !limit = min(51, 31*(a%front -1) + first1loc)
            !limit = 31*(a%front -1) + first1loc

            call sumFractionalComponent(a, a%front, 31*(a%front -1) + first1loc+1, first1loc, 0, r)
        

        end subroutine tosubnormalisedreal


        pure subroutine tonormalisedrealby(a, n, r)
            type(limb_t), intent(in) :: a 
            integer, intent(in) :: n
            real(defReal), intent(inout) :: r 
            real(defReal) :: tr8

            integer(shortInt) tempfront, front
            integer(shortInt) shift
            integer(shortInt) limit, tlimit, i,j
            integer(shortInt) counter


            tr8 = 1_real64
            r = 0_real64

            ! Uses the shift value given to find the position at which the decimal point will be

            
            ! NOTE the line here is equivalent to the code below: (but branchless)
            ! tempfront = a%limbs(max(ceiling(n*1.0/31.0), 1))
            ! ! Covers the fractional part
            ! if (ceiling(n*1.0/31.0) > a%front) then 
            !     tempfront = 0
            ! end if
            front = max(ceiling(n*1.0/31.0), 1) * max(min(a%front-ceiling(n*1.0/31.0)+1, 1), 0)
            tempfront = a%limbs(front)

    

        
            
            !limit = min(51, n)
            shift = n - (n/31)*31

            !call sumFractionalComponent(a, front, n, shift, r)

            tempfront = a%limbs(front)
            
            limit = min(51, n)
            limit = n


            r = 0_real64

            ! Sums the decimals in the limb where the decimal point is located
            do i=shift-1, 0, -1 
                r = r + ((0.5_real64**(shift-i)) * (mod(shiftr(tempfront,i),2)))
                limit = limit - 1
            end do



            tempfront = a%limbs(front)

        
            ! Deducts previously summed values from the limit and sums beyond
            counter = shift+1
            tlimit = limit 
            do i=1, (floor(limit/31.0)+1)

                if (tlimit <= 0) then 
                    exit 
                end if 
                tempfront = a%limbs((floor(n*1.0/31.0)+1-i) * max(min(a%front-(ceiling(n*1.0/31.0)+1-i)+1, 1), 0))
                !tempfront = a%limbs(max(ceiling(n*1.0/31.0)-i,0))

                do j=30, 0, -1
                    if (tlimit <= 0) then 
                        exit 
                    end if 
                

                    r = r + ((0.5_real64**(counter)) * ((mod(shiftr(tempfront, j), 2)))*1_real64)

                    counter = counter + 1
                    tlimit = tlimit - 1 
                end do
            
            end do
      


            call sumWholeComponent(a, front, n, shift, r)

    

        end subroutine tonormalisedrealby


     pure subroutine sumWholeComponent(a, front, limitBy, dpLoc, s) 
            type(limb_t), intent(in) :: a 
            integer, intent(in) :: limitBy
            integer(4), intent(in) :: front, dpLoc
            real(defReal), intent(inout) :: s 
            integer(longInt) :: tempfront
            integer(shortInt) i, j, counter

            tempfront = a%limbs(front)
            counter = 0

            do i=dpLoc, 30
                s = s + ((2_real64**(i-dpLoc)) * (mod(shiftr(tempfront,i), 2_8)))
                counter = counter + 1
            end do

            do i=max(2, (ceiling(limitBy*1.0/31.0)+1)), a%front
                
                tempfront = a%limbs(i)
                do j=0, 30
                    s = s + ((2_real64**(counter)) * (mod(shiftr(tempfront,j),2_8)))
                    counter = counter + 1
                end do 
            end do


        end subroutine sumWholeComponent 


     pure subroutine findFirst1(a, pos)
            type(limb_t), intent(in) :: a 
            integer, intent(inout) :: pos 
            integer(longInt) :: tempfront 
            integer(shortInt) i

            tempfront = a%limbs(a%front)
            pos = 0
            
            do i=0, 30
                if (tempfront == 1) then 
                    pos = i 
                    exit
                else 
                    tempfront = shiftr(tempfront, 1)
                end if 
            end do 

        end subroutine findFirst1


      pure subroutine sumFractionalComponent(a, firstFront, limitBy, dpLoc, normaliseTo, s)
            type(limb_t), intent(in) :: a 
            integer, intent(in) :: dpLoc, limitBy, normaliseTo
            integer(4), intent(in) :: firstFront
            real(defReal), intent(inout) :: s 
            integer(longInt) :: tempfront
            integer(shortInt) i, j, limit, tlimit, counter


            tempfront = a%limbs(firstFront)
            
            limit = min(51, limitBy)
            limit = limitBy


            s = 0

            ! Sums the decimals in the limb where the decimal point is located
            do i=dpLoc-1 + (1-normaliseTo), 0, -1 
                s = s + ((0.5_real64**(dpLoc-i)) * (mod(shiftr(tempfront,i),2_8)))
                limit = limit - 1
            end do

            tempfront = a%limbs(firstfront)

        
            ! Deducts previously summed values from the limit and sums beyond
            counter = dpLoc+1
            tlimit = limit
            do i=1, (ceiling(limit/31.0))

                if (tlimit <= 0) then 
                    exit 
                end if 
                ! Note for reference: this used to be a%front in place of firstFront, didnt cause issues, but logically incorrect
                tempfront = a%limbs(max(1,firstfront-i)) 

                do j=30, 0, -1
                    if (tlimit <= 0) then 
                        exit 
                    end if 

                    s = s + ((0.5_real64**(counter)) * (mod(shiftr(tempfront, j), 2_8)))
                    counter = counter + 1
                    tlimit = tlimit - 1 
                end do
            
            end do


        
        end subroutine sumFractionalComponent




     pure function dividelimbsmixedR(a,n) result(r)
            type(limb_t), intent(in) :: a
            type(limb_t) :: r 
            integer, intent(in) :: n 
            real(8) :: carry
            integer(shortInt) i
            real(8) :: r8 

            if (checkInvalid(a) .or. n == 0) then 
                call setInvalid(r)
                return 
            end if


            r8 = 1
            carry = 0
            r = a
            
            do i = a%front, 1, -1 
                carry = (carry*(int(2, 8)**31)) + ((a%limbs(i))/(n*r8))
                r%limbs(i) = floor(carry)
                carry = carry - floor(carry)
            end do

            r%limbs(1) = r%limbs(1) + ceiling(carry)
            
            call adjustFront(r)

            r%sign = sign(1,n) * a%sign
        
        end function dividelimbsmixedR


        pure function limbequality(a, b) result(r)
          type(limb_t), intent(in) :: a, b
          integer(shortInt)        :: i, maxfront, va, vb
          logical(defBool)         :: r

          if(checkInvalid(a) .or. checkInvalid(b)) then
            r = .false.
            return

          end if

          r = .false.
          maxfront = max(a%front, b%front)

          do i = 1, maxfront 
            if(a % front < i) then 
              va = 0

            else 
              va = a % limbs(i)

            end if 

            if(b % front < i) then 
              vb = 0

            else 
              vb = b % limbs(i)

            end if 

            if(va * a % sign < vb * b % sign .or. vb * b % sign < va * a % sign) return

          end do
          r = .true.
            
        end function limbequality


     pure function limbiszero(a) result (r)
            type(limb_t), intent(in) :: a
            logical :: r 

            if (checkInvalid(a)) then 
                r = .false.
                return 
            end if

            r = a%front == 0 .or. (a%front == 1 .and. a%limbs(a%front) == 0)

            
        end function limbiszero

     pure function limbgeq0(a) result(r)
            type(limb_t), intent(in) :: a
            logical :: r 

            if (checkInvalid(a)) then 
                r = .false.
                return 
            end if

            r = a%limbs(a%front)*a%sign >= 0


        end function limbgeq0


      pure function limbgneq0(a) result(r)
            type(limb_t), intent(in) :: a
            logical :: r 

            if (checkInvalid(a)) then 
                r = .false.
                return 
            end if

            r = a%limbs(a%front)*a%sign > 0


        end function limbgneq0


      pure function limbgneq(a,b) result(r)
            type(limb_t), intent(in) :: a,b
            logical :: r 

            if (checkInvalid(a) ) then 
                r = .false.
                return 
            end if

            r = absgneq(a,b)

            ! Case 1: magnitude of a > b with sign of a positive (a>b or a>-b)
            ! Case 2: magnitude of a < b with both signs negative
            r = (r .and. a%sign == 1) .or. ((.not. r) .and. b%sign == -1 .and. a%sign == -1) & 
                .or. ((.not. r) .and. a%sign == 1 .and. b%sign==-1)


        end function limbgneq


      pure function limbgeq(a,b) result(r)
            type(limb_t), intent(in) :: a,b
            logical :: r 
            integer(longInt) :: zeroIndA, zeroIndB
            integer(shortInt) i 
            integer(shortInt) maxfront 

            if (checkInvalid(a) .or. checkInvalid(b)) then 
                r = .false.
                return 
            end if

            r = .true.

        

            maxfront = max(a%front, b%front)

            do i=maxfront, 1, -1 
                zeroIndA = max(min(a%front+1-i, 1), 0)
                zeroIndB = max(min(b%front+1-i, 1), 0)
                if (a%limbs(i)*zeroIndA*a%sign > b%limbs(i)*zeroIndB*b%sign) then 
                    r = .true. 
                    exit
                else if (a%limbs(i)*zeroIndA*a%sign < b%limbs(i)*zeroIndB*b%sign) then
                    r = .false.
                    exit 
                end if 
            end do 


        end function limbgeq



        




        ! Shifts the given limb by n entries in the array, a shift corresponds to multiplying by 2*31
      pure function shiftbyn(a,n) result (s)
            type(limb_t), intent(in) :: a 
            integer, intent(in) :: n 
            type(limb_t) :: s 

            if (checkInvalid(a) .or. a%front+n > 1977) then 
                call setInvalid(s)
                return 
            end if

            if (a%front == 1 .and. a%limbs(1) == 0) then 
                s = initlimb()
                return 
            end if

            s = a 
            s%limbs = eoshift(s%limbs, (-1)*n, 0, 1)

            s%front = a%front+n
        end function shiftbyn


      pure subroutine assignlimbs(lout, lin)
            type(limb_t), intent(out) :: lout 
            type(limb_t), intent(in) :: lin 

            if (checkInvalid(lin)) then 
                call setInvalid(lout)
                return 
            end if

            lout%limbs = lin%limbs 
            lout%front = lin%front
            lout%sign = lin%sign

        end subroutine assignlimbs




      pure function abslimb(a) result(r)
            type(limb_t), intent(in) :: a 
            type(limb_t) :: r 

            if (checkInvalid(a)) then 
                call setInvalid(r)
                return 
            end if

            r = a 
            r%sign = abs(a%sign)
        
        end function abslimb


        subroutine printlimb(a)
            type(limb_t), intent(in) :: a 
            integer(shortInt) i

            if (checkInvalid(a)) then 
                return 
            end if

            do i=1, a%front
                print *, a%limbs(i)*a%sign
            end do

        end subroutine printlimb


        ! subroutine printFullLimbs(a)
        !     type(limb_t), intent(in) :: a 
        !     character(len=19770) :: num
        !     character(len=10) :: limb
        !     integer(shortInt) i
        !     !limb = ''

        !     do i=1, a%front 
        !         write(limb, '(i10)') a%limbs(i)
        !         num = num // limb
        !     end do

        !     print *, num

        ! end subroutine printFullLimbs







        ! Copies from the given index range of one limb into another to the front
      pure function copyrange(from, s, e) result (to)
            type(limb_t),intent(in) :: from
            integer, intent(in) :: s, e
            type(limb_t) :: to
            integer(shortInt) i 

            if (checkInvalid(from)) then 
                call setInvalid(to)
                return 
            end if

            to = initlimb()
            to%front = e-s+1
            to%sign = from%sign

            do i=s, e
                ! copy here
                to%limbs(i-s+1) = from%limbs(i)

            end do
                
        end function copyrange


end module limb_class


! program multtest
!     use limb_class 
!     !use ratint
!     use, intrinsic :: iso_fortran_env
!     use, intrinsic :: ieee_arithmetic
!     implicit none 
!     integer(shortInt) i
!     real(defReal) :: scalc, sactual, l3, l11, scalc2
!     type(limb_t) ::  v3, v4, v5, t1, t2, v1
!     real(defReal) ::n 
!     type(limb_t) :: d1, d2, d3, d4, d5, d6, dt1, dt2, dt3, v2
!     logical :: result
!    ! type(limb_t) :: sactual, scalc
!     real(defReal) :: va, vb, vc, eval
!     type(limb_t) :: l1, l2
    !type(ratint_t) :: ratint1


    ! v1 = initlimb(0_8)
    ! !call printlimb(v1)
    ! v2 = initlimb(59874_8)
    ! scalc2 = (v1 / v2)
    ! sactual = (0.0_real64 / 59874_8)
    ! print *, scalc2
    ! print *, sactual


!     l11 = 2654.0_real64 / 9988445522.0_real64
!     l1 = initlimb(2654_8)
!     l2 = initlimb(9988445522_8)
!     l3 = l1 / l2 
!     print *, '..'
!     print *, l3

    
!     t1 = initlimb() 
!     t1%front = 2 
!     t1%limbs(1) = 1698360743
!     t1%limbs(2) = 1168591
!     t2 = initlimb()
!     t2%front = 3
!     t2%limbs(1) = 0
!     t2%limbs(2) = 0
!     t2%limbs(3) = 2048
    
!     l3 = t1 / t2
!     print *, '???' 
!     print *, l3

!     v1 = initlimb(int(585297,8))
!     v2 = initlimb(int(692, 8))
!     scalc = (v1 / v2)
!     sactual = (585297.0_real64)/692
!     print *, '---'
!     print *, scalc 
!     print *, sactual
   ! result = (scalc == sactual)

    ! n =  1.0_real64 / (9223372036854775807.0_real64)
    ! v1 = initlimb(1_8)
    ! v2 = initlimb(9223372036854775807_8)

    ! call printlimb(v2)
    ! sactual = n
    ! print *, 'correct'
    ! print *, n
    ! scalc = v1 / v2 
    ! print *, 'calc'
    ! print *, scalc
    ! result = scalc == sactual

!     v1 = initlimbSize(128)
!     v2 = initlimbSize(128)
!     v3 = initlimb()

!     v1%limbs(65) = 1
!     v1%limbs(128) = 1 

!     v3 = v1 * v1
!     call printlimb(v3)

!     do i=1, 128
!         if (v3%limbs(i) /= 0) then 
!             print *,'oh'
!         end if 
!     end do

    

!     v1 = 2654.0_real64 / 9988445522.0_real64
!     v3 = initlimb()
!     v4 = initlimb() 
!     v3%front = 2 
!     v3%limbs(1) = 1698360743
!     v3%limbs(2) = 1168591
!     v4%front = 3 
!     v4%limbs(1) = 0 
!     v4%limbs(2) = 0
!     v4%limbs(3) = 2048
!     eval = v3 / v4
!     result = v1 == eval
    
!     print *, '----'
!     print *, eval 
!     print *, v1

!     v3 = initlimb(int(585297,8))
!     v4 = initlimb(int(692, 8))
!     scalc = (v3 / v4)
!     sactual = (585297.0_real64)/692
!     print *, scalc
!     print *, sactual
!     result = (scalc == sactual)



!     v1 = 1.0_real64 / 75.0_real64
!     print *, '----'
!     print *, v1

!     v3 = initlimb()
!     v3%front = 2 
!     v3%limbs(1) = 887626575
!     v3%limbs(2) = 3579139

!     v4 = initlimb()
!     v4%front = 2
!     v4%limbs(1) = 0 
!     v4%limbs(2) = 268435456
  
!     eval = v3 / v4

    !ratint1 = convert_ieee64(v1)
    !call printRatInt(ratint1)
    !eval = evaluate(ratint1)
    ! print *, 'evaluated:'
    ! print *, eval 

!      n =  -123456_8 * (-987654_8)
!       v1 = initlimb(-123456_8)
!       v2 = -987654_8
!       sactual = initlimb(n)
!       scalc = v1 * v2 
!       result = scalc == sactual
!       call printlimb(scalc)
!       call printlimb(sactual)

!     v1 = initlimb(123456_8)
!     v2 = 987654
!     n = 123456_8 * 987654
!     sactual = initlimb(n)
!     scalc = v1 * v2 
!     result = scalc == sactual
!     call printlimb(scalc)
!     print * ,result


!     v1 = initlimb(5468751218_8)
!     v2 = initlimb(1159874_8)
!     v3 = initlimb(5587985_8)
!     v4 = initlimb(22247_8)
!     v5 = initlimb(5546448787512_8)

!     d1 = ((((v1 * v2) * v3) * v4) * v5)
!    ! print *, '1'
!     d2 = (((v1 * (v2 * v3)) * v4) * v5)
!     !print *, '2'
!     d3 = ((v1 * (v2 * (v3 * v4))) * v5)
!     d4 = (v1 * (v2 * (v3 * (v4 * v5))))
!    ! print *, '4'
!     d5 = ((v1 * v2) * ((v3 * v4) * v5))
!    ! print *, '5'
!     d6 = (((v1 * v2) * v3) * (v4 * v5))
!    ! print *, '6'
!     result = (d1 == d2) .and. (d2 == d3) .and. (d3 == d4) .and. (d5 == d6)
!     print *, result




!    v1 = initlimb(0_8)
!    !call printlimb(v1)
!       v2 = initlimb(0_8)
!       v3 = initlimb(0_8)
!       v4 = initlimb(0_8)
!       v5 = initlimb(0_8)

!     d1 = ((((v1 * v2) * v3) * v4) * v5)
!     d2 = (((v1 * (v2 * v3)) * v4) * v5)
!     d3 = ((v1 * (v2 * (v3 * v4))) * v5)
!     d4 = (v1 * (v2 * (v3 * (v4 * v5))))
!     d5 = ((v1 * v2) * ((v3 * v4) * v5))
!     d6 = (((v1 * v2) * v3) * (v4 * v5))
    ! result = (d1 == d2) .and. (d2 == d3) .and. (d3 == d4) .and. (d5 == d6) 

!      print *, result

    ! !!print *, '1------'
    ! dt1 = (v3 * (v4 * v5))
    ! !!call printlimb(dt1)
    ! !!print *, '1------'
    ! dt2 = ((v3 * v4) * v5)
    ! !!call printlimb(dt2)

    ! print *, '2------'
    ! print *, 'v2'
    ! call printlimb(v2)
    ! print *, ',,,,'
    ! dt1 = (v2 * dt1)
    ! call printlimb(dt1)
    
    ! print *, '2.5-----'
    ! print *, 'v3'
    ! call printlimb(v3)
    ! print *, ',,,,'
    ! dt1 = v3 * dt1
    ! call printlimb(dt1)
    ! print *, v3%front


    ! print *, '2------'
    ! dt2 = (v1 * v2) * dt2
    ! call printlimb(dt2)

    ! print *, '3------'
    ! dt3 = v1 * v2 
    ! call printlimb(dt3)



    ! print *, 'd1'
    ! call printlimb(d1)

    ! print *, 'd2'
    ! call printlimb(d2)

    ! print *, 'd3'
    ! call printlimb(d3)

    ! print *, 'd4'
    ! call printlimb(d4)

    ! print *, 'd5'
    ! call printlimb(d5)

    ! print *, 'd6'
    ! call printlimb(d6)


    ! print *, 'PRODuCt1'

    ! dt1 = (v4 * v5)
    ! call printlimb(dt1)

    ! print *, '----'
    ! print *, 'PRODuCt2'
    ! dt2 = ((v1 * v2) * v3)
    ! call printlimb(dt2)

    ! print *, '----'
    ! print *, 'PRODuCt3'

    ! dt2 = dt2 * dt1

    ! call printlimb(dt2)







!  end program multtest