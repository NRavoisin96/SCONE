module thermalScatteringData_iTest

  use aceNeutronDatabase_class, only : aceNeutronDatabase
  use aceNeutronNuclide_class,  only : aceNeutronNuclide, aceNeutronNuclide_CptrCast
  use CENeutron_class,          only : CENeutron
  use ceNeutronCache_mod,       only : nuclideCache
  use ceNeutronNuclide_inter,   only : ceNeutronNuclide, ceNeutronNuclide_CptrCast
  use CEParticleState_class,    only : buildCEParticleStatePayload
  use dictionary_class,         only : dictionary
  use dictParser_func,          only : charToDict
  use funit
  use materialMenu_mod,         only : mm_init => init
  use neutronXSPackages_class,  only : neutronMicroXSs
  use nuclearDatabase_inter,    only : nuclearDatabase
  use numPrecision

  implicit none

  ! Material definitions
  character(*), parameter :: MAT_INPUT_STR =       &
  & "water {                                       &
  &       moder {1001.03 (h-h2o.49); }             &
  &       composition {                            &
  &         type rawAtomicDensities;               &
  &         nuclides {                             &
  &           1001.03  2.0E-3;                     &
  &           8016.03  1.0E-3;                     &
  &         }                                      &
  &       }                                        &
  &  }                                             &
  &  graphite {                                    &
  &          moder {6012.06  (grph30.46);}         &
  &          composition {                         &
  &            type rawAtomicDensities;            &
  &            nuclides {                          &
  &              6012.06 2.0E-3;                   &
  &            }                                   &
  &          }                                     &
  & }                                              &
  & waterMix {                                     &
  &       temp 500;                                &
  &       moder {1001.03 (h-h2o.50 h-h2o.49); }    &
  &       composition {                            &
  &         type rawAtomicDensities;               &
  &         nuclides {                             &
  &           1001.03  2.0E-3;                     &
  &           8016.03  1.0E-3;                     &
  &         }                                      &
  &       }                                        &
  & }  "

  ! CE Neutron Database specification
  character(*), parameter :: ACE_INPUT_STR = &
  & "aceLibrary ./IntegrationTestFiles/testLib; "

  ! Variables.
  type(aceNeutronDatabase), target :: data

contains
@After
  subroutine cleanUp()

    call data % kill()

  end subroutine cleanUp

  !!
  !! Test the use of thermal scattering libraries
  !!
@Test
  subroutine test_thermalScatteringData()
    class(aceNeutronNuclide), pointer :: ACENuc, C12, H1, H1_2, O16
    class(ceNeutronNuclide), pointer  :: nuc
    class(nuclearDatabase), pointer   :: ptr
    integer(shortInt)                 :: C12_Idx, H1_Idx, H1_2_Idx, i, O16_Idx
    real(defReal)                     :: val
    real(defReal), dimension(2)       :: eBounds, kTBounds
    type(buildCEParticleStatePayload) :: payload
    type(CENeutron)                   :: p
    type(dictionary)                  :: dataDict, matDict
    type(neutronMicroXSs)             :: microXSs
    real(defReal), parameter          :: TOL = 1.0e-6_defReal

    ! Prepare dictionaries
    call charToDict(matDict, MAT_INPUT_STR)
    call charToDict(dataDict, ACE_INPUT_STR)

    ! Build material menu
    call mm_init(matDict)

    ! Initialise data
    ptr => data
    call data % init(dataDict, ptr, silent = .true.)
    call data % activate(([1, 2, 3]), silent = .true.)

    !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !! Perform tests
    !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Get nuclides
    do i = 1, 4
      ACENuc => aceNeutronNuclide_CptrCast(data % getNuclide(i))
      select case(trim(adjustl(ACENuc % ZAID)))
        case('8016.03c')
          O16_Idx = i
        
        case('1001.03c')
          if (ACENuc % stochasticMixing) then
            H1_2_Idx = i

          else
            H1_Idx = i

          end if

        case('6012.06c')
          C12_Idx = i

      end select

    end do
    C12 => aceNeutronNuclide_CptrCast(data % getNuclide(C12_Idx))
    H1 => aceNeutronNuclide_CptrCast(data % getNuclide(H1_Idx))
    H1_2 => aceNeutronNuclide_CptrCast(data % getNuclide(H1_2_Idx))
    O16 => aceNeutronNuclide_CptrCast(data % getNuclide(O16_Idx))

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test scattering tables

    @assertTrue(H1 % hasThData)
    @assertFalse(O16 % hasThData)

    @assertFalse(H1 % thData(1) % hasElastic)
    @assertFalse(H1 % stochasticMixing)
    @assertEqual(size(H1 % thData), 1)

    @assertTrue(C12 % hasThData)
    @assertTrue(C12 % thData(1) % hasElastic)
    @assertTrue(C12 % thData(1) % isCoherent)
    @assertFalse(C12 % stochasticMixing)

    @assertTrue(H1_2 % hasThData)
    @assertTrue(H1_2 % stochasticMixing)
    @assertEqual(size(H1_2 % thData), 2)
    @assertFalse(H1_2 % thData(2) % hasElastic)

    !<><><><><><><><><><><><><><><><><><><><>
    ! Test energy bounds
    eBounds = H1 % thData(1) % getEbounds('inelastic')

    @assertEqual(1.000E-11_defReal, eBounds(1), TOL)
    @assertEqual(1.000E-5_defReal,  eBounds(2), TOL)

    @assertEqual(O16 % SabInel(1), ZERO)
    @assertEqual(O16 % SabInel(2), ZERO)

    eBounds = C12 % thData(1) % getEbounds('elastic')

    @assertEqual(1.000E-11_defReal, eBounds(1), TOL)
    @assertEqual(4.9000E-06,  eBounds(2), TOL)

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test temperature bounds of libraries

    kTbounds = H1_2 % getSabTBounds()
    @assertEqual(4.0812E-8, kTbounds(1), TOL)
    @assertEqual(4.3087E-8, kTbounds(2), TOL)
    
    kTbounds = C12 % getSabTBounds()
    @assertEqual(8.6173E-8, kTbounds(1), TOL)
    @assertEqual(8.6173E-8, kTbounds(2), TOL)

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test sampling from tables

    val = H1 % thData(1) % getInelXS(1.8E-6_defReal)
    @assertEqual(21.018654322_defReal, val, TOL)

    val = H1 % thData(1) % getElXS(1.8E-6_defReal)
    @assertEqual(ZERO, val, TOL)
    
    val = H1_2 % thData(1) % getInelXS(1.8E-6_defReal)
    @assertEqual(21.018654322_defReal, val, TOL)
    
    val = H1_2 % thData(2) % getInelXS(1.8E-6_defReal)
    @assertEqual(21.024875613_defReal, val, TOL)
    
    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test Getting material XSs
    ! water
    payload % uGlobal = [ONE, ZERO, ZERO]
    call p % init(payload)

    ! Total XS of water
    call p % setEnergy(1.8e-6_defReal)
    @assertEqual(ONE, data % getTotalMatXS(p, 1) / 0.0459700882_defReal, TOL)

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test getting XSs
    ! H-1
    nuc => ceNeutronNuclide_CptrCast(data % getNuclide(H1_Idx))
    nuclideCache(H1_Idx) % E_tot = ONE

    call nuc % getMicroXSs(1.8E-6_defReal, ZERO, microXSs, p % getRNGPtr())

    @assertEqual(ONE, 21.05810233858_defReal / microXSs % total, TOL)
    @assertEqual(ONE, 21.01865432_defReal / microXSs % inelasticScatter, TOL)
    @assertEqual(ONE, 3.94480160E-002_defReal / microXSs % capture, TOL)
    @assertEqual(ZERO, microXSs % elasticScatter, TOL)
    @assertEqual(ZERO, microXSs % fission, TOL)
    @assertEqual(ZERO, microXSs % nuFission, TOL)

  end subroutine test_thermalScatteringData

end module thermalScatteringData_iTest