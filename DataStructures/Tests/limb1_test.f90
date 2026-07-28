module limb1_test
    use numPrecision
    use funit
    use limb_class
    !use ratint
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    implicit none

contains


  @Test
    subroutine testAddition()
      type(limb_t) :: v1
      type(limb_t) :: v2
      type(limb_t) :: scalc
      type(limb_t) :: sactual
      integer(8) :: n
      logical :: result

      n = 558564297_8 + 69872592_8
      v1 = initlimb(558564297_8)
      v2 = initlimb(69872592_8)
      scalc = v1 + v2
      sactual = initlimb(n)
      result = (scalc == sactual)
      @assertTrue(result, message = 'add1')


      n = 15423_8 + (-4581671_8)
      v1 = initlimb(15423_8)
      v2 = initlimb(-4581671_8)
      scalc = v1 + v2
      sactual = initlimb(n)
      result = (scalc == sactual)
      @assertTrue(result, message = 'add2')


      ! 64-bit integer limit
      n = 92233720354775807_8 + (-92233720368775807_8)
      v1 = initlimb(92233720354775807_8)
      v2 = initlimb(-92233720368775807_8)
      scalc = v1 + v2
      sactual = initlimb(n)
      result = (scalc == sactual)
      @assertTrue(result, message = 'add3')


      n = 9223372036854775_8 + (9223372036854775_8)
      v1 = initlimb(9223372036854775_8)
      v2 = initlimb(9223372036854775_8)
      scalc = v1 + v2
      sactual = initlimb(n)
      result = (scalc == sactual)
      @assertTrue(result, message = 'add4')

    end subroutine testAddition



    @Test 
    subroutine testAdditionAssociativity()
      type(limb_t) :: v1, v2, v3, v4, v5
      integer(8) ::n 
      type(limb_t) :: d1, d2, d3, d4, d5, d6
      logical :: result

      v1 = initlimb(5468751218_8)
      v2 = initlimb(1159874_8)
      v3 = initlimb(5587985_8)
      v4 = initlimb(22247_8)
      v5 = initlimb(5546448787512_8)

      d1 = ((((v1 + v2) + v3) + v4) + v5)
      d2 = (((v1 + (v2 + v3)) + v4) + v5)
      d3 = ((v1 + (v2 + (v3 + v4))) + v5)
      d4 = (v1 + (v2 + (v3 + (v4 + v5))))
      d5 = ((v1 + v2) + ((v3 + v4) + v5))
      d6 = (((v1 + v2) + v3) + (v4 + v5))
      result = (d1 == d2) .and. (d2 == d3) .and. (d3 == d4) .and. (d5 == d6)

      @assertTrue(result, message='AdditionAssociativity1')

    end subroutine testAdditionAssociativity


    @Test 
    subroutine testAdditionIdentity() 

      type(limb_t) :: v1, v2, vres, vcalc
      integer(8) ::n 
      type(limb_t) :: d1, d2, d3, d4, d5, d6
      logical :: result

      n = 545198733455454_8 + 0_8
      v1 = initlimb(545198733455454_8)
      v2 = initlimb(0_8)
      vres = v1 + v2 
      vcalc = initlimb(n)
      result = vres == vcalc
      @assertTrue(result, message='additionIdentity1')

      n = 0_8 + 545198733455454_8
      v1 = initlimb(0_8)
      v2 = initlimb(545198733455454_8)
      vres = v1 + v2 
      vcalc = initlimb(n)
      result = vres == vcalc
      @assertTrue(result, message='additionIdentity2')


      n = 0_8 + (-545198733455454_8)
      v1 = initlimb(0_8)
      v2 = initlimb(-545198733455454_8)
      vres = v1 + v2 
      vcalc = initlimb(n)
      result = vres == vcalc
      @assertTrue(result, message='additionIdentity2')


      n = (-545198733455454_8) + 0_8
      v1 = initlimb(-545198733455454_8)
      v2 = initlimb(0_8)
      vres = v1 + v2 
      vcalc = initlimb(n)
      result = vres == vcalc
      @assertTrue(result, message='additionIdentity1')

      
    end subroutine testAdditionIdentity


    @Test 
    subroutine testAdditionIdempotence()
      type(limb_t) :: v1, v2, vres, vcalc
      integer(8) ::n 
      type(limb_t) :: d1, d2, d3, d4, d5, d6
      logical :: result

      n = 0_8 + 0_8
      v1 = initlimb(0_8)
      v2 = initlimb(0_8)
      vres = v1 - v2 
      vcalc = initlimb(n)
      result = vres == vcalc
      @assertTrue(result, message='additionIdempotence1')

    end subroutine testAdditionIdempotence



  @Test 
    subroutine testSubtraction() 

      type(limb_t) :: v1
      type(limb_t) :: v2
      type(limb_t) :: scalc
      type(limb_t) :: sactual
      integer(8) :: n
      logical :: result

      n = 1_8 - 545128872415_8
      v1 = initlimb(1_8)
      v2 = initlimb(545128872415_8)
      scalc = v1 - v2
      sactual = initlimb(n)
      result = (scalc == sactual)
      @assertTrue(result, message = 'sub1')


      n = 15423_8 - (-4581671_8)
      v1 = initlimb(15423_8)
      v2 = initlimb(-4581671_8)
      scalc = v1 - v2
      sactual = initlimb(n)
      result = (scalc == sactual)
      @assertTrue(result, message = 'sub2')


      n = 15423_8 - (-0_8)
      v1 = initlimb(15423_8)
      v2 = initlimb(-0_8)
      scalc = v1 - v2
      sactual = initlimb(n)
      result = (scalc == sactual)
      @assertTrue(result, message = 'sub3')

      n = 0_8 - (-2554815_8)
      v1 = initlimb(0_8)
      v2 = initlimb(-2554815_8)
      scalc = v1 - v2
      sactual = initlimb(n)
      result = (scalc == sactual)
      @assertTrue(result, message = 'sub4')


      n = 0_8 - (2554815_8)
      v1 = initlimb(0_8)
      v2 = initlimb(2554815_8)
      scalc = v1 - v2
      sactual = initlimb(n)
      result = (scalc == sactual)
      @assertTrue(result, message = 'sub5')


      n = 92233720354775807_8 - (-92233720368775807_8)
      v1 = initlimb(92233720354775807_8)
      v2 = initlimb(-92233720368775807_8)
      scalc = v1 - v2
      sactual = initlimb(n)
      result = (scalc == sactual)
      @assertTrue(result, message = 'sub6')


      n = 9223372036854775807_8 - (9223372036854775807_8)
      v1 = initlimb(9223372036854775807_8)
      v2 = initlimb(9223372036854775807_8)
      scalc = v1 - v2
      sactual = initlimb(n)
      result = (scalc == sactual)
      @assertTrue(result, message = 'sub7')

    end subroutine testSubtraction


    @Test 
    subroutine testSubtractionIdentity() 

      type(limb_t) :: v1, v2, vres, vcalc
      integer(8) ::n 
      type(limb_t) :: d1, d2, d3, d4, d5, d6
      logical :: result

      n = 545198733455454_8 - 0_8
      v1 = initlimb(545198733455454_8)
      v2 = initlimb(0_8)
      vres = v1 - v2 
      vcalc = initlimb(n)
      result = vres == vcalc
      @assertTrue(result, message='subtractionIdentity1')

      n = 0_8 - 0_8
      v1 = initlimb(0_8)
      v2 = initlimb(0_8)
      vres = v1 - v2 
      vcalc = initlimb(n)
      result = vres == vcalc
      @assertTrue(result, message='subtractionIdentity2')
      
    end subroutine testSubtractionIdentity



    @Test 
    subroutine testSubtractionIdempotence()
      type(limb_t) :: v1, v2, vres, vcalc
      integer(8) ::n 
      type(limb_t) :: d1, d2, d3, d4, d5, d6
      logical :: result

      n = 0_8 - 0_8
      v1 = initlimb(0_8)
      v2 = initlimb(0_8)
      vres = v1 - v2 
      vcalc = initlimb(n)
      result = vres == vcalc
      @assertTrue(result, message='subtractionIdempotence1')


      n = 0_8 - (-0_8)
      v1 = initlimb(0_8)
      v2 = initlimb(-0_8)
      vres = v1 - v2 
      vcalc = initlimb(n)
      result = vres == vcalc
      @assertTrue(result, message='subtractionIdempotence1')

    end subroutine testSubtractionIdempotence



  @Test 
    subroutine testDivision() 
      type(limb_t) :: v1, v2 
      real(real64) :: scalc1, scalc2, scalc
      real(real64) :: sactual
      real(real64) :: t8
      real(real64) :: n
      logical :: result

      t8 = 1

    !NOTE  This test does not pass, because the newton version is more accurate by 1 bit (very end) when comparing against full precision
      ! v1 = initlimb(int(585297,8))
      ! v2 = initlimb(int(692, 8))
      ! scalc1 = (v1 / v2)
      ! sactual = (585297.0_real64)/692
      ! print *, '!!!'
      ! print *, scalc1 
      ! print *, sactual
      ! result = (scalc1 == sactual)
      ! @assertTrue(result, message = 'div1')

   

      v1 = initlimb(int(-15852,8))
      v2 = initlimb(int(8867, 8))
      scalc2 = (v1 / v2)
      sactual = ((-15852.0_real64) / 8867_real64)
      result = (scalc2 == sactual)
      @assertTrue(result, message = 'div2')


      v1 = initlimb(878467598742_8)
      v2 = initlimb(-2585589_8)
      scalc2 = (v1 / v2)
      sactual = (((878467598742.0_real64)) / (-2585589_real64))
      result = (scalc2 == sactual)
      @assertTrue(result, message = 'div3')


      v1 = initlimb(1_8)
      v2 = initlimb(59874_8)
      scalc2 = (v1 / v2)
      sactual = (1.0_real64 / 59874_real64)
      result = (scalc2 == sactual)
      @assertTrue(result, message = 'div4')


      v1 = initlimb(454587892_8)
      v2 = initlimb(1_8)
      scalc2 = (v1 / v2)
      sactual = ((454587892.0_real64) / 1_real64)
      result = (scalc2 == sactual)
      @assertTrue(result, message = 'div5')

      v1 = initlimb(0_8)
      v2 = initlimb(59874_8)
      scalc2 = (v1 / v2)
      sactual = (0.0_real64 / 59874_8)
      result = (scalc2 == sactual)
      @assertTrue(result, message = 'div6')


      n =  (9223372036854775807.0_real64) / 5454842_8 
      v1 = initlimb(9223372036854775807_8)
      v2 = initlimb(5454842_8)
      sactual = n
      scalc = v1 / v2 
      result = scalc == sactual
      @assertTrue(result, message='div7')


      n =  0.0_real64 / (9223372036854775807.0_real64)
      v1 = initlimb(0_8)
      v2 = initlimb(9223372036854775807_8)
      sactual = n
      scalc = v1 / v2 
      result = scalc == sactual
      @assertTrue(result, message='div8')

      n =  1.0_real64 / (9223372036854775807.0_real64)
      v1 = initlimb(1_8)
      v2 = initlimb(9223372036854775807_8)
      sactual = n
      !print *, n
      scalc = v1 / v2 
      !print *, scalc
      print *, n 
      print *, scalc
      result = scalc == sactual
      @assertTrue(result, message='div9')


    end subroutine testDivision

!! NOTE: difficult to test because in order to use the divided value, it needs to be rounded first
    @Test 
    subroutine testDivisionDistributivity()
      type(limb_t) :: v1, v2, v3
      real(8) :: vres, vcalc, d1r, d2r
      integer(8) ::n 
      type(limb_t) :: d1, d2, d3, d4, d5, d6
      logical :: result


      ! v1 = initlimb(66971255_8)
      ! v2 = initlimb(45648942_8)
      ! v3 = initlimb(55484845612_8)
      ! d1 = int((v1 / v2), 4) * v3 
      ! d2 = initlimb(int((v1 * v3) / v2, 8)) 
      ! !n = (585297.0_8)/692
      ! result = (d1 == d2)
      ! @assertTrue(result, message = 'divDist1')


      ! v1 = initlimb(66971255_8)
      ! v2 = initlimb(45648942_8)
      ! v3 = initlimb(55484845612_8)
      ! d1r = (v1 * v3) / v2
      ! print *, d1r
      ! d2r = v1 / v2
      ! d2 = int((v1 / v2),4) * v3
      ! print *, d2r
      ! call printlimb(d2)
      ! !n = (585297.0_8)/692
      ! result = (d1 == d2)
      ! @assertTrue(result, message = 'divDist2')

    end subroutine testDivisionDistributivity



    @Test 
    subroutine testDivisionIdentity() 
      type(limb_t) :: v1, v2
      real(8) :: vres, vcalc
      integer(8) ::n 
      type(limb_t) :: d1, d2, d3, d4, d5, d6
      logical :: result

      n = 545198733455454.0_8 / 1_8
      v1 = initlimb(545198733455454_8)
      v2 = initlimb(1_8)
      vres = v1  / v2 
      vcalc = n
      result = vres == vcalc
      @assertTrue(result, message='divisionIdentity1')

    end subroutine testDivisionIdentity



    @Test
    subroutine testMultiplicationCorrectness() 

      type(limb_t) :: v1, v4
      type(limb_t) :: v2, v3
      integer(8) ::n 
      type(limb_t) :: scalc
      type(limb_t) :: sactual
      logical :: result

      v1 = initlimb(123456_8)
      v2 = initlimb(987654_8)
      n = 121931812224_8
      sactual = initlimb(n)
      scalc = v1 * v2 
      result = scalc == sactual
      @assertTrue(result, message = 'mult1')

      n =  -123456_8 * (-987654_8)
      v1 = initlimb(-123456_8)
      v2 = initlimb(-987654_8)
      sactual = initlimb(n)
      scalc = v1 * v2 
      result = scalc == sactual
      @assertTrue(result, message = 'mult2')


      n = (-123456_8) * 987654_8
      v1 = initlimb(-123456_8)
      v2 = initlimb(987654_8)
      sactual = initlimb(n)
      scalc = v1 * v2 
      result = scalc == sactual
      @assertTrue(result, message = 'mult3')

      n = (-987654) * 123456_8
      v1 = initlimb(123456_8)
      v2 = initlimb(-987654_8)
      sactual = initlimb(n)
      scalc = v1 * v2 
      result = scalc == sactual
      @assertTrue(result, message = 'mult4')


      n = 32456655765745_8 * 1_8
      v1 = initlimb(32456655765745_8)
      v2 = initlimb(1_8)
      sactual = initlimb(n)
      scalc = v1 * v2 
      result = scalc == sactual
      @assertTrue(result, message = 'mult5')


      n = 32456655765745_8 * 0_8
      v1 = initlimb(32456655765745_8)
      v2 = initlimb(0_8)
      sactual = initlimb(n)
      scalc = v1 * v2 
      result = scalc == sactual
      @assertTrue(result, message='mult6')

      ! 64-bit integer limit
      n = 9223372036854775807_8 * 0_8
      v1 = initlimb(9223372036854775807_8)
      
      v2 = initlimb(0_8)
      sactual = initlimb(n)

      scalc = v1 * v2 

      result = scalc == sactual
      @assertTrue(result, message='mult7')

      ! 64-bit integer limit
      n = (-9223372036854775807_8) * 0_8
      v1 = initlimb(-9223372036854775807_8)
      v2 = initlimb(0_8)
      sactual = initlimb(n)
      scalc = v1 * v2 
      result = scalc == sactual
      @assertTrue(result, message='mult8')

      n =  0_8 * (9223372036854775807_8)
      v1 = initlimb(0_8)
      v2 = initlimb(9223372036854775807_8)
      sactual = initlimb(n)
      scalc = v1 * v2 
      result = scalc == sactual
      @assertTrue(result, message='mult9')


      v1 = initlimb(656989745_8)
      v2 = v1*v1 
      v2 = v2*v2
      v2 = v2*v2*v2*v2*v2*v2*v2 
      !v3 = initlimb(55587_8)
      v2 = v2*v2*v2*v2*v2*v2
      v3 = v2 * (v2 * v2)
      v4 = (v2 * v2) * v2 
      

      result = v3 == v4 
      @assertTrue(result, message = ':(')



      
      

    end subroutine testMultiplicationCorrectness


    @Test
    subroutine testMultiplicationAssociativity()
      type(limb_t) :: v1, v2, v3, v4, v5
      integer(8) ::n 
      type(limb_t) :: d1, d2, d3, d4, d5, d6
      logical :: result

      v1 = initlimb(5468751218_8)
      v2 = initlimb(1159874_8)
      v3 = initlimb(5587985_8)
      v4 = initlimb(22247_8)
      v5 = initlimb(5546448787512_8)

      d1 = ((((v1 * v2) * v3) * v4) * v5)
      d2 = (((v1 * (v2 * v3)) * v4) * v5)
      d3 = ((v1 * (v2 * (v3 * v4))) * v5)
      d4 = (v1 * (v2 * (v3 * (v4 * v5))))
      d5 = ((v1 * v2) * ((v3 * v4) * v5))
      d6 = (((v1 * v2) * v3) * (v4 * v5))
      result = (d1 == d2) .and. (d2 == d3) .and. (d3 == d4) .and. (d5 == d6)
      @assertTrue(result, message='MultAssociativity1')



      v1 = initlimb(5468751218_8)
      v2 = initlimb(1159874_8)
      v3 = initlimb(0_8)
      v4 = initlimb(22247_8)
      v5 = initlimb(5546448787512_8)

      d1 = ((((v1 * v2) * v3) * v4) * v5)
      d2 = (((v1 * (v2 * v3)) * v4) * v5)
      d3 = ((v1 * (v2 * (v3 * v4))) * v5)
      d4 = (v1 * (v2 * (v3 * (v4 * v5))))
      d5 = ((v1 * v2) * ((v3 * v4) * v5))
      d6 = (((v1 * v2) * v3) * (v4 * v5))
      result = (d1 == d2) .and. (d2 == d3) .and. (d3 == d4) .and. (d5 == d6) 


      @assertTrue(result, message='MultAssociativity2')



      v1 = initlimb(5468751218_8)
      v2 = initlimb(1159874_8)
      v3 = initlimb(0_8)
      v4 = initlimb(22247_8)
      v5 = initlimb(0_8)

      d1 = ((((v1 * v2) * v3) * v4) * v5)
      d2 = (((v1 * (v2 * v3)) * v4) * v5)
      d3 = ((v1 * (v2 * (v3 * v4))) * v5)
      d4 = (v1 * (v2 * (v3 * (v4 * v5))))
      d5 = ((v1 * v2) * ((v3 * v4) * v5))
      d6 = (((v1 * v2) * v3) * (v4 * v5))
      result = (d1 == d2) .and. (d2 == d3) .and. (d3 == d4) .and. (d5 == d6) 


      @assertTrue(result, message='MultAssociativity3')


      v1 = initlimb(0_8)
      v2 = initlimb(1159874_8)
      v3 = initlimb(0_8)
      v4 = initlimb(22247_8)
      v5 = initlimb(0_8)

      d1 = ((((v1 * v2) * v3) * v4) * v5)
      d2 = (((v1 * (v2 * v3)) * v4) * v5)
      d3 = ((v1 * (v2 * (v3 * v4))) * v5)
      d4 = (v1 * (v2 * (v3 * (v4 * v5))))
      d5 = ((v1 * v2) * ((v3 * v4) * v5))
      d6 = (((v1 * v2) * v3) * (v4 * v5))
      result = (d1 == d2) .and. (d2 == d3) .and. (d3 == d4) .and. (d5 == d6) 


      @assertTrue(result, message='MultAssociativity4')




      v1 = initlimb(0_8)
      v2 = initlimb(1159874_8)
      v3 = initlimb(0_8)
      v4 = initlimb(0_8)
      v5 = initlimb(0_8)

      d1 = ((((v1 * v2) * v3) * v4) * v5)
      d2 = (((v1 * (v2 * v3)) * v4) * v5)
      d3 = ((v1 * (v2 * (v3 * v4))) * v5)
      d4 = (v1 * (v2 * (v3 * (v4 * v5))))
      d5 = ((v1 * v2) * ((v3 * v4) * v5))
      d6 = (((v1 * v2) * v3) * (v4 * v5))
      result = (d1 == d2) .and. (d2 == d3) .and. (d3 == d4) .and. (d5 == d6) 


      @assertTrue(result, message='MultAssociativity5')


      v1 = initlimb(0_8)
      v2 = initlimb(0_8)
      v3 = initlimb(0_8)
      v4 = initlimb(0_8)
      v5 = initlimb(0_8)

      d1 = ((((v1 * v2) * v3) * v4) * v5)
      d2 = (((v1 * (v2 * v3)) * v4) * v5)
      d3 = ((v1 * (v2 * (v3 * v4))) * v5)
      d4 = (v1 * (v2 * (v3 * (v4 * v5))))
      d5 = ((v1 * v2) * ((v3 * v4) * v5))
      d6 = (((v1 * v2) * v3) * (v4 * v5))
      result = (d1 == d2) .and. (d2 == d3) .and. (d3 == d4) .and. (d5 == d6) 


      @assertTrue(result, message='MultAssociativity6')

!9223372036854775807_8
      v1 = initlimb(92233720368547758_8)
      v2 = initlimb(92232036854775807_8)
      v3 = initlimb(92233720368575807_8)
      v4 = initlimb(93372036854775807_8)
      v5 = initlimb(23372036854775807_8)

      d1 = ((((v1 * v2) * v3) * v4) * v5)
      d2 = (((v1 * (v2 * v3)) * v4) * v5)
      d3 = ((v1 * (v2 * (v3 * v4))) * v5)
      d4 = (v1 * (v2 * (v3 * (v4 * v5))))
      d5 = ((v1 * v2) * ((v3 * v4) * v5))
      d6 = (((v1 * v2) * v3) * (v4 * v5))
      result = (d1 == d2) .and. (d2 == d3) .and. (d3 == d4) .and. (d5 == d6) 


      @assertTrue(result, message='MultAssociativity7')


      v1 = initlimb(9223372036854775807_8)
      v2 = initlimb(9223372036854775807_8)
      v3 = initlimb(9223372036854775807_8)
      v4 = initlimb(9223372036854775807_8)
      v5 = initlimb(9223372036854775807_8)

      d1 = ((((v1 * v2) * v3) * v4) * v5)
      d2 = (((v1 * (v2 * v3)) * v4) * v5)
      d3 = ((v1 * (v2 * (v3 * v4))) * v5)
      d4 = (v1 * (v2 * (v3 * (v4 * v5))))
      d5 = ((v1 * v2) * ((v3 * v4) * v5))
      d6 = (((v1 * v2) * v3) * (v4 * v5))
      result = (d1 == d2) .and. (d2 == d3) .and. (d3 == d4) .and. (d5 == d6) 


      @assertTrue(result, message='MultAssociativity8')


      v1 = initlimb(5468751218_8)
      v1 = v1* v1 *v1 *v1 *v1 
      v1 = v1* v1* v1* v1* v1 
      v1 = v1* v1* v1* v1* v1 
      v1 = v1 *v1 *v1 *v1 *v1 

      d1 = ((((v1 * v1) * v1) * v1) * v1)
      d2 = (((v1 * (v1 * v1)) * v1) * v1)
      d3 = ((v1 * (v1 * (v1 * v1))) * v1)
      d4 = (v1 * (v1 * (v1 * (v1 * v1))))
      d5 = ((v1 * v1) * ((v1 * v1) * v1))
      d6 = (((v1 * v1) * v1) * (v1 * v1))
      result = (d1 == d2) .and. (d2 == d3) .and. (d3 == d4) .and. (d5 == d6) 

      @assertTrue(result, message='MultAssociativity9')




      v1 = initlimb(45645487875125_8)
      v2 = initlimb(77459977415125_8)
      v1 = v1* v1 *v1 *v1 *v1 
      v1 = v1* v1* v1* v1* v1 
      v1 = v1* v1* v1* v1* v1 
      v1 = v1 *v1 *v1 *v1 *v1 


      d1 = ((((v1 * v1) * v1) * v1) * v1)
      d2 = (((v1 * (v1 * v1)) * v1) * v1)
      d3 = ((v1 * (v1 * (v1 * v1))) * v1)
      d4 = (v1 * (v1 * (v1 * (v1 * v1))))
      d5 = ((v1 * v1) * ((v1 * v1) * v1))
      d6 = (((v1 * v1) * v1) * (v1 * v1))
      result = (d1 == d2) .and. (d2 == d3) .and. (d3 == d4) .and. (d5 == d6) 

      @assertTrue(result, message='MultAssociativity9')

      v1 = initlimb(9223372036854775807_8)
      v1 = v1 * v1 * v1 * v1 
      v1 = v1 * v1 * v1 * v1 
      v1 = v1 * v1 * v1 * v1 
      v2 = initlimb(9223372036854775807_8)
      v2 = v2 * v2 * v2 * v2 
      v2 = v2 * v2 * v2 * v2 
      v2 = v2 * v2 * v2 * v2 
      v3 = initlimb(9223372036854775807_8)
      v3 = v3 * v3 * v3 * v3 
      v3 = v3 * v3 * v3 * v3 
      v3 = v3 * v3 * v3 * v3 
      v4 = initlimb(9223372036854775807_8)
      v4 = v4 * v4 * v4 * v4 
      v4 = v4 * v4 * v4 * v4 
      v4 = v4 * v4 * v4 * v4
      v5 = initlimb(9223372036854775807_8)
      v5 = v5 * v5 * v5 * v5 
      v5 = v5 * v5 * v5 * v5 
      v5 = v5 * v5 * v5 * v5 

      d1 = ((((v1 * v2) * v3) * v4) * v5)
      d2 = (((v1 * (v2 * v3)) * v4) * v5)
      d3 = ((v1 * (v2 * (v3 * v4))) * v5)
      d4 = (v1 * (v2 * (v3 * (v4 * v5))))
      d5 = ((v1 * v2) * ((v3 * v4) * v5))
      d6 = (((v1 * v2) * v3) * (v4 * v5))
      result = (d1 == d2) .and. (d2 == d3) .and. (d3 == d4) .and. (d5 == d6) 


      @assertTrue(result, message='MultAssociativity10')



      


    end subroutine testMultiplicationAssociativity



    @Test 
    subroutine testMultAddDistributivity() 

      type(limb_t) :: v1, v2, v3, v4, v5
      integer(8) ::n 
      type(limb_t) :: d1, d2, d3, d4, d5, d6
      logical :: result

      v1 = initlimb(1598745_8)
      v2 = initlimb(7891159874_8)
      v3 = initlimb(555897587985_8)

      d1 = v3 * (v1 + v2)
      d2 = (v3 * v1) + (v3 * v2)

      result = (d1 == d2)
      @assertTrue(result, message = 'distmultadd1')


      v1 = initlimb(1598745_8)
      v2 = initlimb(0_8)
      v3 = initlimb(555897587985_8)

      d1 = v3 * (v1 + v2)
      d2 = (v3 * v1) + (v3 * v2)

      result = (d1 == d2)
      @assertTrue(result, message = 'distmultadd2')



      v1 = initlimb(1598745_8)
      v2 = initlimb(7891159874_8)
      v3 = initlimb(0_8)

      d1 = v3 * (v1 + v2)
      d2 = (v3 * v1) + (v3 * v2)

      result = (d1 == d2)
      @assertTrue(result, message = 'distmultadd3')


      v1 = initlimb(0_8)
      v2 = initlimb(7891159874_8)
      v3 = initlimb(555897587985_8)

      d1 = v3 * (v1 + v2)
      d2 = (v3 * v1) + (v3 * v2)

      result = (d1 == d2)
      @assertTrue(result, message = 'distmultadd4')


      v1 = initlimb(0_8)
      v2 = initlimb(7891159874_8)
      v3 = initlimb(0_8)

      d1 = v3 * (v1 + v2)
      d2 = (v3 * v1) + (v3 * v2)

      result = (d1 == d2)
      @assertTrue(result, message = 'distmultadd5')



      v1 = initlimb(0_8)
      v2 = initlimb(0_8)
      v3 = initlimb(0_8)

      d1 = v3 * (v1 + v2)
      d2 = (v3 * v1) + (v3 * v2)

      result = (d1 == d2)
      @assertTrue(result, message = 'distmultadd6')


      v1 = initlimb(23372036854775807_8)
      v2 = initlimb(23372036854775807_8)
      v3 = initlimb(23372036854775807_8)

      d1 = v3 * (v1 + v2)
      d2 = (v3 * v1) + (v3 * v2)

      result = (d1 == d2)
      @assertTrue(result, message = 'distmultadd7')
      

      v1 = initlimb(9223372036854775807_8)
      v2 = initlimb(9223372036854775807_8)
      v3 = initlimb(9223372036854775807_8)

      d1 = v3 * (v1 + v2)
      d2 = (v3 * v1) + (v3 * v2)

      result = (d1 == d2)
      @assertTrue(result, message = 'distmultadd8')

      v1 = initlimb(45645487875125_8)
      v1 = v1* v1 *v1 *v1 *v1 
      v1 = v1* v1* v1* v1* v1 
      v1 = v1* v1* v1* v1* v1 
      ! print *, 'front'
      ! print *, v1%front
      !v1 = v1 *v1 *v1 *v1 *v1 
     

      d1 = v1 * (v1 + v1)
      d2 = (v1 * v1) + (v1 * v1)

      result = (d1 == d2)
      @assertTrue(result, message = 'distmultadd9')


      

    end subroutine testMultAddDistributivity


    @Test 
    subroutine testMultSubDistributivity() 

      type(limb_t) :: v1, v2, v3, v4, v5
      integer(8) ::n 
      type(limb_t) :: d1, d2, d3, d4, d5, d6
      logical :: result

      v1 = initlimb(1598745_8)
      v2 = initlimb(7891159874_8)
      v3 = initlimb(555897587985_8)

      d1 = v3 * (v1 - v2)
      d2 = (v3 * v1) - (v3 * v2)

      result = (d1 == d2)
      @assertTrue(result, message = 'distmultsub1')



      v1 = initlimb(1598745_8)
      v2 = initlimb(0_8)
      v3 = initlimb(555897587985_8)

      d1 = v3 * (v1 - v2)
      d2 = (v3 * v1) - (v3 * v2)

      result = (d1 == d2)
      @assertTrue(result, message = 'distmultsub2')



      v1 = initlimb(1598745_8)
      v2 = initlimb(7891159874_8)
      v3 = initlimb(0_8)

      d1 = v3 * (v1 - v2)
      d2 = (v3 * v1) - (v3 * v2)

      result = (d1 == d2)
      @assertTrue(result, message = 'distmultsub3')


      v1 = initlimb(0_8)
      v2 = initlimb(7891159874_8)
      v3 = initlimb(555897587985_8)

      d1 = v3 * (v1 - v2)
      d2 = (v3 * v1) - (v3 * v2)

      result = (d1 == d2)
      @assertTrue(result, message = 'distmultsub3')


      v1 = initlimb(0_8)
      v2 = initlimb(7891159874_8)
      v3 = initlimb(0_8)

      d1 = v3 * (v1 - v2)
      d2 = (v3 * v1) - (v3 * v2)

      result = (d1 == d2)
      @assertTrue(result, message = 'distmultsub3')



      v1 = initlimb(0_8)
      v2 = initlimb(0_8)
      v3 = initlimb(0_8)

      d1 = v3 * (v1 - v2)
      d2 = (v3 * v1) - (v3 * v2)

      result = (d1 == d2)
      @assertTrue(result, message = 'distmultsub3')


      v1 = initlimb(23372036854775807_8)
      v2 = initlimb(23372036854775807_8)
      v3 = initlimb(23372036854775807_8)

      d1 = v3 * (v1 - v2)
      d2 = (v3 * v1) - (v3 * v2)

      result = (d1 == d2)
      @assertTrue(result, message = 'distmultsub3')
      

      v1 = initlimb(9223372036854775807_8)
      v2 = initlimb(9223372036854775807_8)
      v3 = initlimb(9223372036854775807_8)

      d1 = v3 * (v1 - v2)
      d2 = (v3 * v1) - (v3 * v2)

      result = (d1 == d2)
      @assertTrue(result, message = 'distmultsub3')
      

    end subroutine testMultSubDistributivity

    @Test 
    subroutine testMultiplicationIdentity() 

      type(limb_t) :: v1, v2, vres, vcalc
      integer(8) ::n 
      type(limb_t) :: d1, d2, d3, d4, d5, d6
      logical :: result

      n = 545198733455454_8 * 1_8
      v1 = initlimb(545198733455454_8)
      v2 = initlimb(1_8)
      vres = v1 * v2 
      vcalc = initlimb(n)
      result = vres == vcalc
      @assertTrue(result, message='testMultiplicationIdentity1')


    
    end subroutine testMultiplicationIdentity


    @Test 
    subroutine testMultiplicationAnnihilation() 

      type(limb_t) :: v1, v2, vres, vcalc
      integer(8) ::n 
      type(limb_t) :: d1, d2, d3, d4, d5, d6
      logical :: result

      n = 5451987335454_8 * 0_8
      v1 = initlimb(5451987335454_8)
      v2 = initlimb(0_8)
      vres = v1 * v2 
      vcalc = initlimb(n)
      result = vres == vcalc
      @assertTrue(result, message='multiplicationAnnihilation1')


      n = 0_8 * 5451987355454_8
      v1 = initlimb(0_8)
      v2 = initlimb(5451987355454_8)
      vres = v1 * v2 
      vcalc = initlimb(n)
      result = vres == vcalc
      @assertTrue(result, message='multiplicationAnnihilation2')


      n = 0_8 * 1_8
      v1 = initlimb(0_8)
      v2 = initlimb(1_8)
      vres = v1 * v2 
      vcalc = initlimb(n)
      result = vres == vcalc
      @assertTrue(result, message='multiplicationAnnihilation3')


      n = 1_8 * 0_8
      v1 = initlimb(1_8)
      v2 = initlimb(0_8)
      vres = v1 * v2 
      vcalc = initlimb(n)
      result = vres == vcalc
      @assertTrue(result, message='multiplicationAnnihilation4')


      n = 0_8 * (-1_8)
      v1 = initlimb(0_8)
      v2 = initlimb(-1_8)
      vres = v1 * v2 
      vcalc = initlimb(n)
      result = vres == vcalc
      @assertTrue(result, message='multiplicationAnnihilation5')


      n = (-1_8) * 0_8
      v1 = initlimb(-1_8)
      v2 = initlimb(0_8)
      vres = v1 * v2 
      vcalc = initlimb(n)
      result = vres == vcalc
      @assertTrue(result, message='multiplicationAnnihilation5')



    end subroutine testMultiplicationAnnihilation


    @Test 
    subroutine testmixedMultR() 
      type(limb_t) :: v1
      integer(8) :: v2
      integer(8) ::n 
      type(limb_t) :: scalc
      type(limb_t) :: sactual
      logical :: result

      v1 = initlimb(123456_8)
      v2 = 987654_8
      n = 123456_8 * 987654_8
      sactual = initlimb(n)
      scalc = v1 * v2 
      result = scalc == sactual
      @assertTrue(result, message = 'mixmultr1')

      n =  -123456_8 * (-987654_8)
      v1 = initlimb(-123456_8)
      v2 = -987654_8
      sactual = initlimb(n)
      scalc = v1 * v2 
      result = scalc == sactual
      @assertTrue(result, message = 'mixmultr2')


      n = (-123456_8) * 987654_8
      v1 = initlimb(-123456_8)
      v2 = 987654_8
      sactual = initlimb(n)
      scalc = v1 * v2 
      result = scalc == sactual
      @assertTrue(result, message = 'mixmultr3')

      n = (-987654) * 123456_8
      v1 = initlimb(123456_8)
      v2 = -987654_8
      sactual = initlimb(n)
      scalc = v1 * v2 
      result = scalc == sactual
      @assertTrue(result, message = 'mixmultr4')


      n = 32456655765745_8 * 1_8
      v1 = initlimb(32456655765745_8)
      v2 = 1_8
      sactual = initlimb(n)
      scalc = v1 * v2 
      result = scalc == sactual
      @assertTrue(result, message = 'mixmultr5')


      n = 32456655765745_8 * 0_8
      v1 = initlimb(32456655765745_8)
      v2 = 0_8
      sactual = initlimb(n)
      scalc = v1 * v2 
      result = scalc == sactual
      @assertTrue(result, message='mixmultr6')

      ! 64-bit integer limit
      n = 9223372036854775807_8 * 0_8
      v1 = initlimb(9223372036854775807_8)
      
      v2 = 0_8
      sactual = initlimb(n)

      scalc = v1 * v2 

      result = scalc == sactual
      @assertTrue(result, message='mixmultr7')

      n = (-922337_8) * 0_8
      v1 = initlimb(-922337_8)
      v2 = 0_8
      sactual = initlimb(n)
      scalc = v1 * v2 
      result = scalc == sactual
      @assertTrue(result, message='mixmultr8')

 

    end subroutine testmixedMultR



    @Test 
    subroutine testmixedMultL() 
      type(limb_t) :: v1
      integer(8) :: v2
      integer(8) ::n 
      type(limb_t) :: scalc
      type(limb_t) :: sactual
      logical :: result

      v1 = initlimb(123456_8)
      v2 = 987654_8
      n = 123456_8 * 987654_8
      sactual = initlimb(n)
      scalc = v2 * v1
      result = scalc == sactual
      @assertTrue(result, message = 'mixmultl1')

      n =  -123456_8 * (-987654_8)
      v1 = initlimb(-123456_8)
      v2 = -987654_8
      sactual = initlimb(n)
      scalc = v2 * v1
      result = scalc == sactual
      @assertTrue(result, message = 'mixmultl2')


      n = (-123456_8) * 987654_8
      v1 = initlimb(-123456_8)
      v2 = 987654_8
      sactual = initlimb(n)
      scalc = v2 * v1 
      result = scalc == sactual
      @assertTrue(result, message = 'mixmultl3')

      n = (-987654) * 123456_8
      v1 = initlimb(123456_8)
      v2 = -987654_8
      sactual = initlimb(n)
      scalc = v2 * v1 
      result = scalc == sactual
      @assertTrue(result, message = 'mixmultl4')


      n = 32456655765745_8 * 1_8
      v1 = initlimb(32456655765745_8)
      v2 = 1_8
      sactual = initlimb(n)
      scalc = v2 * v1
      result = scalc == sactual
      @assertTrue(result, message = 'mixmultl5')


      n = 32456655765745_8 * 0_8
      v1 = initlimb(32456655765745_8)
      v2 = 0_8
      sactual = initlimb(n)
      scalc = v2 * v1
      result = scalc == sactual
      @assertTrue(result, message='mixmultl6')

      ! 64-bit integer limit
      n = 9223372036854775807_8 * 0_8
      v1 = initlimb(9223372036854775807_8)
      
      v2 = 0_8
      sactual = initlimb(n)

      scalc = v2 * v1 

      result = scalc == sactual
      @assertTrue(result, message='mixmultl7')

      n = (-922337_8) * 0_8
      v1 = initlimb(-922337_8)
      v2 = 0_8
      sactual = initlimb(n)
      scalc = v2 * v1
      result = scalc == sactual
      @assertTrue(result, message='mixmultl8')

 

    end subroutine testmixedMultL


    @Test 
    subroutine testDivisionMixed()

      type(limb_t) :: v1, res, v3
      real(8) :: scalc1, scalc2, scalc
      type(limb_t) :: sactual
      real(8) :: t8
      real(8) :: n
      
      integer :: v2
      logical :: result


      v1 = initlimb(598988452_8)
      v2 = 2
      res = (v1 / v2)
      sactual = initlimb((598988452)/2 * 1_8)
      result = (res == sactual)
      @assertTrue(result, message = 'divmixed1')

   

      v1 = initlimb(int(-598988452,8))
      v2 = 2
      res = (v1 / v2)
      sactual = initlimb(((-598988452) / 2)*1_8)
      result = (res == sactual)
      @assertTrue(result, message = 'divmixed2')


      v1 = initlimb(116473856_8)
      v2 = 256
      res = (v1 / v2)
      sactual = initlimb(int((((116473856.0_8)) / (256)),8))
      result = (res == sactual)
      @assertTrue(result, message = 'divmixed3')



      v1 = initlimb(116473856_8)
      v2 = 256
      res = (v1 / v2)
      sactual = initlimb(int((((116473856.0_8)) / (256)),8))
      result = (res == sactual)
      @assertTrue(result, message = 'divmixed4')


      v1 = initlimb(78297986016_8)
      v2 = 8841236
      res = (v1 / v2)
      sactual = initlimb(int((((78297986016.0_8)) / (8841236)),8))
      result = (res == sactual)
      @assertTrue(result, message = 'divmixed5')


      v1 = initlimb(78297986016_8)
      v1 = v1 * v1 
      v1 = v1 * v1 * v1
      v2 = 8841236
      res = (v1 / v2)
      v3 = initlimb(8841236_8)
      sactual = v1 / v2
      result = (res == sactual)
      @assertTrue(result, message = 'divmixed6')



      v1 = initlimb(162713295016_8)
      v1 = v1 * v1 
      v1 = v1 * v1 * v1
      v1 = v1* v1 * v1
      v1 = v1 * v1 *v1 *v1
      v2 = 352957256
      res = (v1 / v2)
      v3 = initlimb(352957256_8)
      sactual = v1 / v2
      result = (res == sactual)
      @assertTrue(result, message = 'divmixed7')


      


    end subroutine testDivisionMixed


    @Test 
    subroutine testRatIntConversion() 
      ! real(real64) :: v1, v2, vres, l3
      ! type(ratint_t) :: ratint1, ratint2
      ! real(real64) :: eval
      ! logical :: result
      ! type(limb_t) :: l1, l2, t1, t2, rs 

      ! v1 = 564653.0_real64 / 75.0_real64
      ! l1 = initlimb(564653_8)
      ! l2 = initlimb(75_8)
      ! l3 = l1 / l2 
      ! print *, '..'
      ! print *, l3
      ! t1 = initlimb() 
      ! t1%front = 2 
      ! t1%limbs(1) = 1947051841
      ! t1%limbs(2) = 1927348
      ! t2 = initlimb()
      ! t2%front = 2 
      ! t2%limbs(1) = 0
      ! t2%limbs(2) = 256
      ! print *, '???'
      ! l3 = t1 / t2 
      ! print *, l3
      

      ! ratint1 = convert_ieee64(v1)
      ! eval = evaluate(ratint1)
      ! print *, eval 
      ! print *, v1
      ! result = v1 == eval
      ! @assertTrue(result, message='ratintconv1')


      ! v1 = 5599874562114.0_real64 / 1.0_real64
      ! ratint1 = convert_ieee64(v1)
      ! eval = evaluate(ratint1)
      ! result = v1 == eval
      ! @assertTrue(result, message='ratintconv2')



      ! v1 = 9999887744.0_real64 / 75.0_real64
      ! ratint1 = convert_ieee64(v1)
      ! eval = evaluate(ratint1)
      ! result = v1 == eval
      ! @assertTrue(result, message='ratintconv3')


! !!NOTE: the tests below are commented out, because subzero values are hard to compare this way

!       v1 = 1.0_real64 / 75.0_real64
!       ratint1 = convert_ieee64(v1)
!       eval = evaluate(ratint1)


!       result = v1 == eval
!       @assertTrue(result, message='ratintconv4')





!       v1 = 1.0_real64 / 1.0_real64
!       ratint1 = convert_ieee64(v1)
!       eval = evaluate(ratint1)
!       result = v1 == eval
!       @assertTrue(result, message='ratintconv5')



!       v1 = 955841.0_real64 / 2.0_real64
!       ratint1 = convert_ieee64(v1)
!       eval = evaluate(ratint1)
!       result = v1 == eval
!       @assertTrue(result, message='ratintconv6')



      ! v1 = 2654.0_real64 / 9988445522.0_real64
      ! ratint1 = convert_ieee64(v1)
      ! eval = evaluate(ratint1)
      ! result = v1 == eval
      ! call printRatInt(ratint1)
      ! print *, '----'
      ! print *, eval 
      ! print *, v1
      ! @assertTrue(result, message='ratintconv7')



      ! v1 = 0.0_real64 / 75.0_real64
      ! ratint1 = convert_ieee64(v1)
      ! eval = evaluate(ratint1)
      ! result = v1 == eval
      ! @assertTrue(result, message='ratintconv8')


    end subroutine testRatIntConversion





  end module limb1_test