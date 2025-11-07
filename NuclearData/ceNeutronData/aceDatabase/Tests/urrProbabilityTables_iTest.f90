module urrProbabilityTables_iTest

  use aceNeutronDatabase_class, only : aceNeutronDatabase
  use aceNeutronNuclide_class,  only : aceNeutronNuclide, aceNeutronNuclide_CptrCast
  use ceNeutronCache_mod,       only : zaidCache, nuclideCache
  use ceNeutronNuclide_inter,   only : ceNeutronNuclide, ceNeutronNuclide_CptrCast
  use dictionary_class,         only : dictionary
  use dictParser_func,          only : charToDict
  use funit
  use materialMenu_mod,         only : mm_init => init
  use neutronXSPackages_class,  only : neutronMicroXSs
  use nuclearDatabase_inter,    only : nuclearDatabase
  use numPrecision

  implicit none

  ! Material definitions
  character(*), parameter :: MAT_INPUT_STR = &
  & " uo2  {                   &
  &        composition {       &
  &        92235.03 1.0E-3;    &
  &        8016.03  2.0E-3;    &
  &        }                   &
  &      }"

  ! CE Neutron Database specification
  character(*), parameter :: ACE_INPUT_STR = &
  & "aceLibrary ./IntegrationTestFiles/testLib; ures 1 ; majorant 1; "

  ! Variables.
  type(aceNeutronDatabase), target :: data

contains
@After
  subroutine cleanUp()

    call data % kill()

  end subroutine cleanUp

  !!
  !! Test the use of probability tables
  !!
@Test
  subroutine test_urrProbabilityTables()
    class(aceNeutronNuclide), pointer :: ACENuc, O16, U235
    class(ceNeutronNuclide), pointer  :: nuc
    class(nuclearDatabase), pointer   :: ptr
    integer(shortInt)                 :: i, O16_Idx, U235_Idx
    real(defReal), dimension(2)       :: eBounds
    real(defReal), dimension(3)       :: val
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
    call data % activate([1], silent = .true.)

    !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !! Perform tests
    !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Get nuclides
    do i = 1, 2
      ACENuc => aceNeutronNuclide_CptrCast(data % getNuclide(i))
      select case(trim(adjustl(ACENuc % ZAID)))
      case('92235.03c')
        U235_Idx = i

      case('8016.03c')
        O16_Idx = i

      end select

    end do
    U235 => aceNeutronNuclide_CptrCast(data % getNuclide(U235_Idx))
    O16 => aceNeutronNuclide_CptrCast(data % getNuclide(O16_Idx))

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test probability tables

    @assertTrue(U235 % hasProbTab)
    @assertFalse(O16 % hasProbTab)

    !<><><><><><><><><><><><><><><><><><><><>
    ! Test energy bounds

    eBounds = U235 % probTab % getEbounds()

    @assertEqual(2.25E-3_defReal, eBounds(1), TOL)
    @assertEqual(2.5E-2_defReal,  eBounds(2), TOL)

    @assertEqual(O16 % urrE(1), ZERO)
    @assertEqual(O16 % urrE(2), ZERO)

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test sampling from tables

    call U235 % probTab % sampleXSs(9.1E-3_defReal, 0.347_defReal, val)

    @assertEqual(0.98499622_defReal, val(1), TOL)
    @assertEqual(0.83939802_defReal, val(2), TOL)
    @assertEqual(0.8515398_defReal, val(3), TOL)

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test getting XSs

    ! U-235
    nuc  => ceNeutronNuclide_CptrCast(data % getNuclide(U235_Idx))
    zaidCache(U235_Idx) % E  = 9.1E-3_defReal
    zaidCache(U235_Idx) % xi = 0.347_defReal
    nuclideCache(U235_Idx) % E_tot = ONE

    call nuc % getMicroXSs(9.1E-3_defReal, ZERO, microXSs)

    @assertEqual(ONE, 15.317184903738868_defReal/ microXSs % total, TOL)
    @assertEqual(ONE, 11.662135262310867_defReal/ microXSs % elasticScatter, TOL)
    @assertEqual(ONE, 0.5743300000E-5_defReal   / microXSs % inelasticScatter, TOL)
    @assertEqual(ONE, 0.999051523404001_defReal / microXSs % capture, TOL)
    @assertEqual(ONE, 2.655992374724002_defReal / microXSs % fission, TOL)
    @assertEqual(ONE, 6.462838469906821_defReal / microXSs % nuFission, TOL)

  end subroutine test_urrProbabilityTables

end module urrProbabilityTables_iTest