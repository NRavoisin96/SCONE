module ceNeutronMaterial_class

  use CENeutron_class,         only : castCENeutronPtr, CENeutron
  use ceNeutronCache_mod,      only : materialCache, nuclideCache
  use ceNeutronDatabase_inter, only : ceNeutronDatabase
  use ceNeutronNuclide_inter,  only : ceNeutronNuclide, ceNeutronNuclide_CptrCast
  use coordList_class,         only : coordList
  use genericProcedures,       only : fatalError, numToChar
  use materialHandle_inter,    only : materialHandle
  use neutronMaterial_inter,   only : neutronMaterial
  use neutronXsPackages_class, only : neutronMacroXSs
  use numPrecision
  use RNG_class,               only : RNG
  use scalarField_inter,       only : getScalarFieldValue, scalarField
  use scatteringKernels_func,  only : relativeEnergy_constXS, dopplerCorrectionFactor
  use transportObject_inter,   only : transportObject
  use universalVariables

  implicit none
  private

  ! Parameters.
  integer(shortInt), parameter :: FISSION = 1, SCATTER = 2, SCATTER_WITH_FISSION = 3

  !!
  !! Public Pointer Cast
  !!
  public ceNeutronMaterial_CptrCast, ceNeutronMaterial_TptrCast

  !!
  !! An abstract class that represent all CE Neutron Material Data
  !!
  !! Exist mainly in order to decouple caching logic from the database implementation
  !! so there is no need to repeat it in every database type. Thus it will be easier to
  !! mantain and optimise.
  !!
  !! Note that a material without any composition is not allowed.
  !! Makes no assumption about the range of nucIdx. Allows for -ve values
  !!
  !! Interface:
  !!   materialHandle Interface
  !!   getMacroXSs -> return package of macroscopic XSs directly from Energy and RNG
  !!   set         -> Set data related to material by keyword association
  !!   setComposition -> Set composition of material from densities and nucIdxs
  !!   sampleNuclide  -> sample collision nuclide
  !!   sampleFission  -> sample collision nuclide given that fission reaction has happened
  !!   sampleScatter  -> sample collision nuclide given that Scattering has happened
  !!   sampleScatterWithFission -> sample collision nuclide given ther scatter or fission has
  !!     happened
  !!
  type, public, extends(neutronMaterial) :: ceNeutronMaterial
    character(nameLen)                           :: name = ''
    integer(shortInt)                            :: matIdx = 0
    class(ceNeutronDatabase), pointer            :: data => null()
    real(defReal), dimension(:), allocatable     :: dens
    integer(shortInt), dimension(:), allocatable :: nuclides
    logical(defBool)                             :: fissile = .false., hasTMS = .false.
    real(defReal)                                :: eLowerURR = ZERO, eUpperSab = ZERO, kT = ZERO
  contains
    ! Superclass procedures
    procedure                  :: isFissile
    procedure                  :: kill
    generic                    :: getMacroXSs => getMacroXSs_byE
    procedure, non_overridable :: getMacroXSs_byE
    procedure                  :: getMacroXSs_byP
    procedure, non_overridable :: sampleFission
    procedure, non_overridable :: sampleNuclide
    procedure, private         :: sampleReactionNuclide
    procedure, non_overridable :: sampleScatter
    procedure, non_overridable :: sampleScatterWithFission
    procedure, non_overridable :: set
    procedure, non_overridable :: setComposition
    procedure                  :: useTMS
  end type ceNeutronMaterial

contains
  !!
  !! Cast materialHandle pointer to ceNeutronMaterial pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class materialHandle
  !!
  !! Result:
  !!   Null is source is not of ceNeutronMaterial
  !!   Pointer to source if source is ceNeutronMaterial class
  !!
  pure function ceNeutronMaterial_CptrCast(source) result(ptr)
    class(materialHandle), pointer, intent(in) :: source
    class(ceNeutronMaterial), pointer          :: ptr

    select type(source)
      class is(ceNeutronMaterial)
        ptr => source

      class default
        ptr => null()

    end select

  end function ceNeutronMaterial_CptrCast

  !!
  !! Cast materialHandle pointer to ceNeutronMaterial pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class materialHandle
  !!
  !! Result:
  !!   Null is source is not of ceNeutronMaterial
  !!   Pointer to source if source is ceNeutronMaterial class
  !!
  pure function ceNeutronMaterial_TptrCast(source) result(ptr)
    class(materialHandle), pointer, intent(in) :: source
    type(ceNeutronMaterial), pointer            :: ptr

    select type(source)
      type is(ceNeutronMaterial)
        ptr => source

      class default
        ptr => null()
        
    end select

  end function ceNeutronMaterial_TptrCast

  !!
  !! Return Macroscopic XSs for the material
  !!
  !! Args:
  !!   densityFactor [in] -> local density scaling factor. Negative values should be ignored.
  !!   E [in]             -> Requested energy [MeV]
  !!   kT [in]            -> local temperature [MeV]. Negative values should be ignored.
  !!   rand [inout]       -> Random Number Generator
  !!   xss [out]          -> Cross section package to store the data
  !!
  !! Errors:
  !!   fatalError if E is out-of-bounds for the stored data
  !!
  subroutine getMacroXSs_byE(self, densityFactor, E, kT, rand, xss)
    class(ceNeutronMaterial), intent(in) :: self
    real(defReal), intent(in)            :: densityFactor, E, kT
    type(RNG), intent(inout), optional   :: rand
    type(neutronMacroXSs), intent(out)   :: xss

    ! Check Cache and update if needed
    associate(matCache => materialCache(self % matIdx))
      if (any([matCache % densityFactor_tot, matCache % E_tail, matCache % E_tot, matCache % kT_tot] /= &
              [densityFactor, E, E, kT])) call self % data % updateMacroXSs(densityFactor, E, kT, self % matIdx, rand)
      xss = matCache % xss

    end associate

  end subroutine getMacroXSs_byE

  !!
  !! Return Macroscopic XSs for the material given particle
  !!
  !! See neutronMaterial_inter for details
  !!
  subroutine getMacroXSs_byP(self, object, xss)
    class(ceNeutronMaterial), intent(in) :: self
    class(transportObject), intent(in)   :: object
    type(neutronMacroXSs), intent(out)   :: xss
    class(CENeutron), pointer            :: CENeutronPtr
    type(coordList), pointer             :: coordListPtr

    CENeutronPtr => castCENeutronPtr(object)
    coordListPtr => CENeutronPtr % getCoordsPtr()
    call self % getMacroXSs(getScalarFieldValue(nameDensity, ONE, coordListPtr), CENeutronPtr % getEnergy(), &
                            getScalarFieldValue(nameTemperature, self % kT, coordListPtr, kBoltzmann_MeV), &
                            CENeutronPtr % getRNGPtr(), xss)

  end subroutine getMacroXSs_byP

  !!
  !! Return .true. if material is fissile
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   .true. if fissile, .false. otherwise
  !!
  !! Errors:
  !!   None
  !!
  elemental function isFissile(self) result(isIt)
    class(ceNeutronMaterial), intent(in) :: self
    logical(defBool)                     :: isIt

    isIt = self % fissile

  end function isFissile

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(ceNeutronMaterial), intent(inout) :: self

    self % matIdx = 0
    self % kT = ZERO
    self % data => null()
    if (allocated(self % dens)) deallocate(self % dens)
    if (allocated(self % nuclides)) deallocate(self % nuclides)
    self % fissile = .false.
    self % hasTMS = .false.
    self % eUpperSab = ZERO
    self % eLowerURR = ZERO

  end subroutine kill

  !!
  !! Sample fission nuclide given that a fission neutron was produced
  !!
  !! Basically samples from P(nucIdx| fission neutron produced in material)
  !! Useful when generating fission sites
  !!
  !! As such it uses nu*sigma_f
  !!
  !! For a non-fissile material return nucIdx <= 0 !
  !!
  !! Args:
  !!   densityFactor [in] -> Local density scaling factor [-].
  !!   E [in]             -> incident energy [MeV]
  !!   kT [in]            -> Boltzmann constant * local temperature [MeV].
  !!   rand [inout]       -> random number generator
  !!
  !! Result:
  !!   nucIdx of the sampled nuclide for collision.
  !!
  !! Errors:
  !!   fatalError if sampling fails for some reason (E.G. random number > 1)
  !!   fatalError if E is out-of-bounds of the present data
  !!   Returns nucIdx <= if material is not fissile
  !!
  function sampleFission(self, densityFactor, E, kT, rand) result(nucIdx)
    class(ceNeutronMaterial), intent(in) :: self
    real(defReal), intent(in)            :: densityFactor, E, kT
    type(RNG), intent(inout)             :: rand
    integer(shortInt)                    :: nucIdx

    ! Short-cut for nonFissile material
    nucIdx = 0
    if (.not. self % fissile) return
    call self % sampleReactionNuclide(FISSION, densityFactor, E, kT, rand, nucIdx)

  end function sampleFission

  !!
  !! Sample collision nuclide at energy E
  !!
  !! This function randomly determines the exact nuclide for a collision
  !! It uses nuclide total XSs to determine nuclide
  !!
  !! Args:
  !!   E [in]       -> incident energy [MeV]
  !!   rand [inout] -> random number generator
  !!   nucIdx [out] -> sampled nuclide index
  !!   eOut [out]   -> relative energy between neutron and target (may be /= E in case of TMS)
  !!
  !! Errors:
  !!   fatalError if sampling fails for some reason (E.G. random number > 1)
  !!   fatalError if E is out-of-bounds of the present data
  !!
  subroutine sampleNuclide(self, densityFactor, E, kT, rand, nucIdx, eOut)
    class(ceNeutronMaterial), intent(in) :: self
    real(defReal), intent(in)            :: densityFactor, E, kT
    type(RNG), intent(inout)            :: rand
    integer(shortInt), intent(out)       :: nucIdx
    real(defReal), intent(out)           :: eOut
    class(ceNeutronNuclide), pointer     :: nuc
    integer(shortInt)                    :: i
    real(defReal)                        :: P_acc, eMin, eMax, eRel, &
                                            trackMatXS, totNucXS, dens, randomNumber
    character(*), parameter              :: HERE = 'sampleNuclide (ceNeutronMaterial_class.f90)'

    ! Get material tracking XS
    associate(matCache => materialCache(self % matIdx))
      if (any([matCache % densityFactor_track, matCache % E_track, matCache % kT_track] /= [densityFactor, E, kT])) &
      call self % data % updateTrackMatXS(densityFactor, E, kT, self % matIdx, rand)
      call rand % generate(trackMatXS, mult = matCache % trackXS)

      ! Loop over nuclides
      do i = 1, size(self % nuclides)
        nucIdx = self % nuclides(i)
        dens = self % dens(i) * densityFactor

        associate (nucCache => nuclideCache(nucIdx))
          ! Retrieve nuclide XS from cache
          if (self % useTMS(E)) then
            ! If the material is using TMS, the nuclide temperature majorant is needed
            ! The check for the right values stored in cache happens inside the subroutine
            call self % data % updateTotalTempNucXS(E, kT, nucIdx)
            totNucXS = nucCache % tempMajXS * nucCache % doppCorr

          else
            ! Update nuclide cache if needed
            if (nucCache % E_tot /= E .or. matCache % kT_tot /= kT) call self % data % updateTotalNucXS(E, nucIdx, kT, rand)
            totNucXS = nucCache % xss % total

          end if
          trackMatXS = trackMatXS - totNucXS * dens

          ! Nuclide temporarily accepted: check TMS condition
          if (trackMatXS < ZERO) then
            ! Save energy to be used to sample reaction
            eOut = E

            if (self % useTMS(E)) then
              ! If the material is using TMS, retrieve nuclide and nuclide information
              nuc => ceNeutronNuclide_CptrCast(self % data % getNuclide(nucIdx))
              if (.not. associated(nuc)) call fatalError(HERE, 'Failed to retrieve CE Neutron Nuclide')

              ! Sample relative energy
              eRel = relativeEnergy_constXS(E, nuc % getMass(), nucCache % deltakT, rand)

              ! Call through system minimum and maximum energies
              call self % data % energyBounds(eMin, eMax)

              ! Ensure relative energy is within energy bounds
              eRel = min(max(eRel, eMin), eMax)

              ! Get relative energy nuclide cross section
              call self % data % updateTotalNucXS(eRel, nucIdx, kT, rand)

              ! Calculate acceptance probability using ratio of relative energy xs to temperature majorant
              P_acc = nucCache % xss % total * nucCache % doppCorr / totNucXS

              ! Accept or reject the sampled nuclide
              call rand % generate(randomNumber)
              if (P_acc <= randomNumber) nucIdx = REJECTED

              ! Overwrite energy to be used to sample reaction
              eOut = eRel

            end if

            ! Exit function, return the sampled nucIdx
            return

          end if

        end associate

      end do

    end associate

    ! Print error message as the inversion failed
    call fatalError(HERE,'Nuclide sampling loop failed to terminate')

  end subroutine sampleNuclide

  !!
  !!
  !!
  subroutine sampleReactionNuclide(self, what, densityFactor, E, kT, rand, nuclideIdx)
    class(ceNeutronMaterial), intent(in) :: self
    integer(shortInt), intent(in)        :: what
    real(defReal), intent(in)            :: densityFactor, E, kT
    type(RNG), intent(inout)             :: rand
    integer(shortInt), intent(out)       :: nuclideIdx
    class(ceNeutronNuclide), pointer     :: nuc
    integer(shortInt)                    :: i
    real(defReal)                        :: factor, mult, subtractedAmount, xs
    character(*), parameter              :: HERE = 'sampleReactionNuclide (ceNeutronMaterial_class.f90)'

    associate(matCache => materialCache(self % matIdx))
      ! Update cache (without checking the energy) to get the correct results with TMS. The relative energy flag 
      ! cached is cleaned to make sure cross sections are updated.
      matCache % E_rel = ZERO
      call self % data % updateMacroXSs(densityFactor, E, kT, self % matIdx, rand)

      ! Compute multiplication factor and generate macroscopic cross section.
      select case(what)
        case(FISSION)
          mult = matCache % xss % nuFission

        case(SCATTER)
          mult = matCache % xss % elasticScatter + matCache % xss % inelasticScatter

        case(SCATTER_WITH_FISSION)
          mult = matCache % xss % elasticScatter + matCache % xss % inelasticScatter + matCache % xss % fission

        case default
          call fatalError(HERE, 'Invalid reaction channel.')

      end select
      call rand % generate(xs, mult = mult)

      ! Loop over all nuclides in the material.
      do i = 1, size(self % nuclides)
        ! The nuclide cache should be at the right energy after updating the material
        ! In the case of TMS where the macro xss are at a relative energy, the nuclide
        ! xss to be used are at the relative energy just sampled.
        nuclideIdx = self % nuclides(i)
        associate(nucCache => nuclideCache(nuclideIdx))
          select case(what)
            case(FISSION)
              subtractedAmount = nucCache % xss % nuFission

            case(SCATTER)
              subtractedAmount = nucCache % xss % elasticScatter + nucCache % xss % inelasticScatter

            case(SCATTER_WITH_FISSION)
              subtractedAmount = nucCache % xss % elasticScatter + nucCache % xss % inelasticScatter + nucCache % xss % fission

            case default
              call fatalError(HERE, 'Invalid reaction channel.')

          end select

          ! Compute Doppler correction factor then update cross section.
          factor = ONE
          if (self % useTMS(E)) then
            nuc => ceNeutronNuclide_CptrCast(self % data % getNuclide(nuclideIdx))
            if (.not. associated(nuc)) call fatalError(HERE, 'Failed to retrieve CE neutron Nuclide.')
            factor = dopplerCorrectionFactor(E, nuc % getMass(), kT - nuc % getkT())

          end if
          xs = xs - subtractedAmount * factor * self % dens(i) * densityFactor
          if (xs < ZERO) return

        end associate

      end do

      ! Print error message if the inversion failed.
      call fatalError(HERE, 'Nuclide sampling loop failed to terminate.')

    end associate

  end subroutine sampleReactionNuclide

  !!
  !! Sample collision nuclide given that any scattering has happened
  !! Treats fission as a capture!
  !!
  !! For a pure-absorbing material return nucIdx <= 0 !
  !!
  !! Args:
  !!   densityFactor [in] -> Local density scaling factor [-].
  !!   E [in]             -> incident energy [MeV]
  !!   kT [in]            -> Boltzmann constant * local temperature [MeV].
  !!   rand [inout] -> random number generator
  !!
  !! Result:
  !!   nucIdx of the sampled nuclide for collision.
  !!
  !! Errors:
  !!   fatalError if sampling fails for some reason (E.G. random number > 1)
  !!   fatalError if E is out-of-bounds of the present data
  !!   Returns nucIdx <= if material is a pure-absorber (with fission as absorbtion)
  !!
  function sampleScatter(self, densityFactor, E, kT, rand) result(nucIdx)
    class(ceNeutronMaterial), intent(in) :: self
    real(defReal), intent(in)            :: densityFactor, E, kT
    type(RNG), intent(inout)             :: rand
    integer(shortInt)                    :: nucIdx

    call self % sampleReactionNuclide(SCATTER, densityFactor, E, kT, rand, nucIdx)

  end function sampleScatter

  !!
  !! Sample collision nuclide given that any scattering or fission has happened
  !! Treats fission as a scattering!
  !!
  !! For a pure-capture material return nucIdx <= 0 !
  !!
  !! Args:
  !!   densityFactor [in] -> Local density scaling factor [-].
  !!   E [in]             -> incident energy [MeV]
  !!   kT [in]            -> Boltzmann constant * local temperature [MeV].
  !!   rand [inout] -> random number generator
  !!
  !! Result:
  !!   nucIdx of the sampled nuclide for collision.
  !!
  !! Errors:
  !!   fatalError if sampling fails for some reason (E.G. random number > 1)
  !!   fatalError if E is out-of-bounds of the present data
  !!   Returns nucIdx <= if material is a pure-capture (with fission as scattering)
  !!
  function sampleScatterWithFission(self, densityFactor, E, kT, rand) result(nucIdx)
    class(ceNeutronMaterial), intent(in) :: self
    real(defReal), intent(in)            :: densityFactor, E, kT
    type(RNG), intent(inout)             :: rand
    integer(shortInt)                    :: nucIdx

    call self % sampleReactionNuclide(SCATTER_WITH_FISSION, densityFactor, E, kT, rand, nucIdx)

  end function sampleScatterWithFission

  !!
  !! Set matIdx, pointer to a database and fissile flag
  !!
  !! All arguments are optional. Use with keyword association e.g.
  !!   call mat % set(matIdx = 7)
  !!
  !! Use this procedure ONLY during build. NEVER during transport.
  !! IT IS NOT THREAD SAFE!
  !!
  !! NOTE: eUpperSab and eLowerURR are fed by the aceNeutronDatabase, and they are the
  !! strictest (respectively highest and lowest) energy limits among all nuclides
  !! in the material composition.
  !!
  !! Args:
  !!   name [in]      -> material name
  !!   matIdx [in]    -> material index
  !!   database [in]  -> pointer to a database that updates XSs on the ceNeutronCache
  !!   fissile [in]   -> flag indicating whether fission data is present
  !!   hasTMS [in]    -> flag indicating whether TMS is on
  !!   temp [in]      -> TMS material temperature
  !!   eUpperSab [in] -> upper energy of S(a,b) range in the material
  !!   eLowerURR [in] -> lower energy of ures range in the material
  !!
  subroutine set(self, name, matIdx, database, fissile, hasTMS, temp, eUpperSab, eLowerURR)
    class(ceNeutronMaterial), intent(inout)                 :: self
    character(nameLen), intent(in), optional                :: name
    integer(shortInt), intent(in), optional                 :: matIdx
    class(ceNeutronDatabase), pointer, optional, intent(in) :: database
    logical(defBool), intent(in), optional                  :: fissile
    logical(defBool), intent(in), optional                  :: hasTMS
    real(defReal), intent(in), optional                     :: temp
    real(defReal), intent(in), optional                     :: eUpperSab
    real(defReal), intent(in), optional                     :: eLowerURR
    character(*), parameter :: Here = 'set (ceNeutronMaterial_class.f90)'

    if (present(name)) self % name = name
    if (present(database)) self % data => database
    if (present(fissile)) self % fissile = fissile
    if (present(matIdx)) self % matIdx = matIdx
    if (present(hasTMS)) self % hasTMS = hasTMS
    if (present(temp)) self % kT = temp * kBoltzmann_MeV

    if (present(eUpperSab)) then
      if (eUpperSab < ZERO) call fatalError (Here, 'Upper Sab energy limit of material '&
                                             &//numToChar(matIdx)//' is negative.')
      self % eUpperSab  = eUpperSab
    end if

    if (present(eLowerURR)) then
      if (eLowerURR < ZERO) call fatalError (Here, 'Lower URR energy limit of material '&
                                             &//numToChar(matIdx)//' is negative.')
      self % eLowerURR  = eLowerURR
    end if

    ! Check to make sure URR energy and S(a,b) energy do not overlap
    ! Possible in principle, but would be quite strange...
    if (present(eLowerURR) .and. present(eUpperSab)) then
      if (eUpperSab > eLowerUrr) call fatalError(Here,self % name//&
              ' has an overlap in URR and S(alpha, beta) energy ranges. Dodgy data?')
    end if

  end subroutine set

  !!
  !! Set composition of the material in terms of nucIdx and atomic density
  !!
  !! Use this procedure ONLY during build. NEVER during transport.
  !! IT IS NOT THREAD SAFE!
  !!
  !! Args:
  !!   dens    [in] -> array of atomic densities [1/barn/cm] of nuclides
  !!   nucIdxs [in] -> correpsonding array with nucIdxs
  !!
  !! Errors:
  !!   FatalError if arrays have different size
  !!   FatalError if dens contains -ve values
  !!   FatalError if dens has size of 0 -> no composition
  !!
  subroutine setComposition(self, dens, nucIdxs)
    class(ceNeutronMaterial), intent(inout)     :: self
    real(defReal), dimension(:), intent(in)     :: dens
    integer(shortInt), dimension(:), intent(in) :: nucIdxs
    character(*), parameter :: Here = 'setComposition (ceNeutronMaterial_class.f90)'

    ! Check input
    if (size(dens) /= size(nucIdxs)) call fatalError(Here,'Different sizes of density and nuclide vector')
    if (any(dens < ZERO)) call fatalError(Here,'-ve nuclide densities are present')
    if (size(dens) == 0)  call fatalError(Here,'Empty composition is not allowed')

    ! Clean any current content
    if (allocated(self % dens))     deallocate(self % dens)
    if (allocated(self % nuclides)) deallocate(self % nuclides)

    ! Load values
    self % dens     = dens
    self % nuclides = nucIdxs

  end subroutine setComposition

  !!
  !! Return .true. if TMS is on in the material and the provided energy is not
  !! within the ures or S(a,b) range of any nuclide in the material.
  !!
  !! Args:
  !!   E [in] -> test energy
  !!
  !! Result:
  !!   .true. if conditions for using TMS are satisfied, .false. otherwise
  !!
  !! Errors:
  !!   None
  !!
  elemental function useTMS(self, E) result(shouldIt)
    class(ceNeutronMaterial), intent(in) :: self
    real(defReal), intent(in)            :: E
    logical(defBool)                     :: shouldIt

    shouldIt = self % hasTMS .and. E < self % eLowerURR .and. self % eUpperSab < E

  end function useTMS

end module ceNeutronMaterial_class