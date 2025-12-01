module ceNeutronDatabase_inter

  use CENeutron_class,       only : castCENeutronPtr, CENeutron
  use ceNeutronCache_mod,    only : materialCache, majorantCache, trackingCache
  use charMap_class,         only : charMap
  use coordList_class,       only : coordList
  use genericProcedures,     only : fatalError
  use intMap_class,          only : intMap
  use materialHandle_inter,  only : materialHandle
  use nuclearDatabase_inter, only : nuclearDatabase
  use nuclideHandle_inter,   only : nuclideHandle
  use numPrecision
  use reactionHandle_inter,  only : reactionHandle
  use RNG_class,             only : RNG
  use scalarField_inter,     only : getScalarFieldValue
  use transportObject_inter, only : transportObject
  use universalVariables,    only : kBoltzmann_MeV, MAJORANT_XS, MATERIAL_XS, nameDensity, nameTemperature, TRACKING_XS

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public ceNeutronDatabase_CptrCast

  !!
  !! An abstract base class for all nulcear databases that support CE Neutron
  !!
  !! Its primary goal is to contain CE Neutron caching logic so there is
  !! no need to reproduce it in each database implementation.
  !!
  !! It is also used by material and nuclide handles for CE Neutron data to order an
  !! update of XSs on the cache
  !!
  !! Public Members:
  !!   mapDBRCnuc -> map to link indexes of DBRC nuclides with their corresponding 0K
  !!
  !! Interface:
  !!   nuclearDatabase Interface
  !!   energyBounds       -> return maximum and minimum energy
  !!   updateTotalMatXS   -> update Total Material XS on CE Neutron Cache
  !!   updateMajorantXS   -> update Majorant XS on CE Neutron Cache
  !!   updateMacroXSs     -> update Macroscopic XSs for a selected material
  !!   updateTotalXS      -> update Total XS for a selected nuclide
  !!   updateMicroXSs     -> update Microscopic XSs for a selected nuclide
  !!   getScattMicroMajXS -> returns nuclide elastic scattering temperature majorant (for DBRC)
  !!   updateTotalTempNucXS -> returns nuclide total temperature majorant (for TMS)
  !!
  type, public, abstract, extends(nuclearDatabase) :: ceNeutronDatabase
    type(intMap) :: mapDBRCnuc
  contains
    procedure(energyBounds), deferred         :: energyBounds
    procedure                                 :: getMajorantXS
    procedure(getMaterial_kT), deferred       :: getMaterial_kT
    procedure(getScattMicroMajXS), deferred   :: getScattMicroMajXS
    procedure                                 :: getTotalMatXS
    procedure                                 :: getTrackingXS
    procedure                                 :: getTrackMatXS
    procedure(updateMacroXSs), deferred       :: updateMacroXSs
    procedure(updateMajorantXS), deferred     :: updateMajorantXS
    procedure(updateMicroXSs), deferred       :: updateMicroXSs
    procedure(updateTotalMatXS), deferred     :: updateTotalMatXS
    procedure(updateTotalNucXS), deferred     :: updateTotalNucXS
    procedure(updateTotalTempNucXS), deferred :: updateTotalTempNucXS
    procedure(updateTrackMatXS), deferred     :: updateTrackMatXS
  end type ceNeutronDatabase

  abstract interface
    !!
    !! Return energy bounds for data in the database
    !!
    !! eMin and eMax are minimun and maximumum energy such that data
    !! for ALL nuclides if avalible
    !!
    !! Args:
    !!   eMin [out] -> minimum value of energy [MeV]
    !!   eMax [out] -> maximum value of energy [MeV]
    !!
    !! Errors:
    !!   None
    !!
    subroutine energyBounds(self, eMin, eMax)
      import :: ceNeutronDatabase, defReal
      class(ceNeutronDatabase), intent(in) :: self
      real(defReal), intent(out)           :: eMin
      real(defReal), intent(out)           :: eMax
    end subroutine energyBounds

    !!
    !!
    !!
    function getMaterial_kT(self, matIdx) result(kT)
      import                               :: ceNeutronDatabase, defReal, shortInt
      class(ceNeutronDatabase), intent(in) :: self
      integer(shortInt), intent(in)        :: matIdx
      real(defReal)                        :: kT
    end function getMaterial_kT

    !!
    !! Function to get the elastic scattering majorant cross section in a nuclide
    !! over a certain energy range, defined as a function of a given temperature
    !!
    !! NOTE: This function is called by the collision operator to apply DBRC; nucIdx
    !!       should correspond to a nuclide with temperature 0K, while kT is the
    !!       temperature of the target nuclide the neutron is colliding with
    !!
    !! Args:
    !!   A  [in]   -> Nuclide atomic weight ratio
    !!   kT [in]   -> Thermal energy of nuclide [MeV]
    !!   E  [in]   -> Energy of neutron incident to target for which majorant needs to be found
    !!   maj [out] -> Majorant cross section
    !!
    function getScattMicroMajXS(self, E, kT, A, nucIdx) result(maj)
      import :: ceNeutronDatabase, defReal, shortInt
      class(ceNeutronDatabase), intent(in) :: self
      real(defReal), intent(in)            :: E
      real(defReal), intent(in)            :: kT
      real(defReal), intent(in)            :: A
      integer(shortInt), intent(in)        :: nucIdx
      real(defReal)                        :: maj
    end function getScattMicroMajXS

    !!
    !! Make sure that the macroscopic XSs for the material with matIdx are set
    !! to energy E in ceNeutronCache
    !!
    !! ANY CHANGE in ceNeutronCache is POSSIBLE
    !!   E.G. Extra materials may be set to energy E as well
    !!
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   matIdx [in]  -> material index that needs to be updated
    !!   rand [inout] -> random number generator
    !!
    subroutine updateMacroXSs(self, densityFactor, E, kT, matIdx, rand)
      import :: ceNeutronDatabase, defReal, shortInt, RNG
      class(ceNeutronDatabase), intent(in) :: self
      real(defReal), intent(in)            :: densityFactor, E, kT
      integer(shortInt), intent(in)        :: matIdx
      type(RNG), optional, intent(inout)  :: rand
    end subroutine updateMacroXSs

    !!
    !! Make sure that the majorant of ALL Active materials is at energy E
    !! in ceNeutronCache
    !!
    !! ANY CHANGE in ceNeutronCache is POSSIBLE
    !!   E.G. All material XSs may be updated to energy E
    !!
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   rand [inout] -> random number generator
    !!
    subroutine updateMajorantXS(self, E, rand)
      import :: ceNeutronDatabase, defReal, RNG
      class(ceNeutronDatabase), intent(in) :: self
      real(defReal), intent(in)            :: E
      type(RNG), optional, intent(inout)  :: rand
    end subroutine updateMajorantXS

    !!
    !! Make sure that the microscopic XSs for the nuclide with nucIdx are set
    !! to energy E in ceNeutronCache
    !!
    !! ANY CHANGE in ceNeutronCache is POSSIBLE
    !!   E.G. Extra nuclides may be set to energy E as well
    !!
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   nucIdx [in]  -> material index that needs to be updated
    !!   kT [in]      -> thermal energy of material [MeV]
    !!   rand [inout] -> random number generator
    !!
    subroutine updateMicroXSs(self, E, nucIdx, kT, rand)
      import :: ceNeutronDatabase, defReal, shortInt, RNG
      class(ceNeutronDatabase), intent(in) :: self
      real(defReal), intent(in)            :: E
      integer(shortInt), intent(in)        :: nucIdx
      real(defReal), intent(in)            :: kT
      type(RNG), optional, intent(inout)  :: rand
    end subroutine updateMicroXSs

    !!
    !! Make sure that totalXS of material with matIdx is at energy E
    !! in ceNeutronCache
    !!
    !! ANY CHANGE in ceNeutronCache is POSSIBLE
    !!   E.G. All material XSs may be updated to energy E
    !!
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   matIdx [in]  -> material index that needs to be updated
    !!   rand [inout] -> random number generator
    !!
    subroutine updateTotalMatXS(self, densityFactor, E, kT, matIdx, rand)
      import :: ceNeutronDatabase, defReal, shortInt, RNG
      class(ceNeutronDatabase), intent(in) :: self
      real(defReal), intent(in)            :: densityFactor, E, kT
      integer(shortInt), intent(in)        :: matIdx
      type(RNG), optional, intent(inout)  :: rand
    end subroutine updateTotalMatXS

    !!
    !! Make sure that totalXS of nuclide with nucIdx is at energy E
    !! in ceNeutronCache
    !!
    !! ANY CHANGE in ceNeutronCache is POSSIBLE
    !!   E.G. All nuclide XSs may be updated to energy E
    !!
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   nucIdx [in]  -> nuclide index that needs to be updated
    !!   kT [in]      -> thermal energy of material [MeV]
    !!   rand [inout] -> random number generator
    !!
    subroutine updateTotalNucXS(self, E, nucIdx, kT, rand)
      import :: ceNeutronDatabase, defReal, shortInt, RNG
      class(ceNeutronDatabase), intent(in) :: self
      real(defReal), intent(in)            :: E
      integer(shortInt), intent(in)        :: nucIdx
      real(defReal), intent(in)            :: kT
      type(RNG), intent(inout), optional  :: rand
    end subroutine updateTotalNucXS

    !!
    !! Subroutine to retrieve the nuclide total majorant cross section over a range
    !! of relative energies a nuclide can see given the material temperature. The
    !! energy range is calculated based on the nuclide and material kT, on the atomic
    !! mass of the nuclide, and the incident neutron energy
    !!
    !! The nuclide cache is updated with the 'temperature' majorant value, incident
    !! neutron energy, deltakT and Doppler correction factor
    !!
    !! Args:
    !!   E [in]      -> required energy [MeV]
    !!   kT [in]     -> thermal energy of TMS material [MeV]
    !!   nucIdx [in] -> material index that needs to be updated
    !!
    !! Errors:
    !!   FatalError if material kT is smaller than the nuclide kT
    !!
    subroutine updateTotalTempNucXS(self, E, kT, nucIdx)
      import :: ceNeutronDatabase, defReal, shortInt
      class(ceNeutronDatabase), intent(in) :: self
      real(defReal), intent(in)            :: E
      real(defReal), intent(in)            :: kT
      integer(shortInt), intent(in)        :: nucIdx
    end subroutine updateTotalTempNucXS

    !!
    !! Make sure that trackXS of material with matIdx is at energy E = E_track
    !! in ceNeutronCache
    !!
    !! The tracking xs corresponds to the material total cross section unless TMS
    !! is used. In that case, this is the material temperature majorant xs.
    !!
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   matIdx [in]  -> material index that needs to be updated
    !!   rand [inout] -> random number generator
    !!
    subroutine updateTrackMatXS(self, densityFactor, E, kT, matIdx, rand)
      import :: ceNeutronDatabase, defReal, shortInt, RNG
      class(ceNeutronDatabase), intent(in) :: self
      real(defReal), intent(in)            :: densityFactor, E, kT
      integer(shortInt), intent(in)        :: matIdx
      type(RNG), optional, intent(inout)  :: rand
    end subroutine updateTrackMatXS

  end interface

contains
  !!
  !! Cast nuclearDatabase pointer to ceNeutronDatabase pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class nuclearDatabase
  !!
  !! Result:
  !!   Null is source is not of ceNuclearDatabase class
  !!   Target points to source if source is ceNuclearDatabase class
  !!
  pure function ceNeutronDatabase_CptrCast(source) result(ptr)
    class(nuclearDatabase), pointer, intent(in) :: source
    class(ceNeutronDatabase), pointer           :: ptr

    select type(source)
      class is(ceNeutronDatabase)
        ptr => source

      class default
        ptr => null()

    end select

  end function ceNeutronDatabase_CptrCast

  !!
  !! Return Majorant XS
  !!
  !! See nuclearDatabase_inter for details
  !!
  !! Error:
  !!   fatalError if particle is not CE Neutron
  !!
  function getMajorantXS(self, object) result(xs)
    class(ceNeutronDatabase), intent(inout) :: self
    class(transportObject), intent(in)      :: object
    real(defReal)                           :: energy, xs
    type(CENeutron), pointer                :: CENeutronPtr

    ! Check dynamic type of the particle
    CENeutronPtr => castCENeutronPtr(object)
    energy = CENeutronPtr % getEnergy()
    associate(majCache => majorantCache(1))
      ! Check Cache and update if needed
      if (majorantCache(1) % E /= energy) call self % updateMajorantXS(energy, CENeutronPtr % getRNGPtr())

      ! Return Cross-Section
      xs = majorantCache(1) % xs

    end associate

  end function getMajorantXS

  !!
  !! Return Total XS for matIdx
  !!
  !! See nuclearDatabase_inter for details!
  !!
  !! Error:
  !!   fatalError if particle is not CE Neutron
  !!
  function getTotalMatXS(self, object, matIdx) result(xs)
    class(ceNeutronDatabase), intent(inout) :: self
    class(transportObject), intent(in)      :: object
    integer(shortInt), intent(in)           :: matIdx
    real(defReal)                           :: densityFactor, energy, kT, xs
    type(CENeutron), pointer                :: CENeutronPtr
    type(coordList), pointer                :: coordListPtr

    ! Check dynamic type of the particle.
    CENeutronPtr => castCENeutronPtr(object)
    coordListPtr => CENeutronPtr % getCoordsPtr()

    associate(matCache => materialCache(matIdx))
      ! Check Cache and update if needed
      densityFactor = getScalarFieldValue(nameDensity, ONE, coordListPtr)
      energy = CENeutronPtr % getEnergy()
      kT = getScalarFieldValue(nameTemperature, self % getMaterial_kT(matIdx), coordListPtr, kBoltzmann_MeV)
      if (any([matCache % densityFactor_tot, matCache % E_tot, matCache % kT_tot] /= [densityFactor, energy, kT])) &
      call self % updateTotalMatXS(densityFactor, energy, kT, matIdx, CENeutronPtr % getRNGPtr())

      ! Return Cross-Section
      xs = matCache % xss % total

    end associate

  end function getTotalMatXS

  !!
  !! Return tracking XS requested
  !!
  !! See nuclearDatabase_inter for details!
  !!
  !! Error:
  !!   fatalError if particle is not CE Neutron
  !!
  function getTrackingXS(self, object, matIdx, what) result(xs)
    class(ceNeutronDatabase), intent(inout) :: self
    class(transportObject), intent(in)      :: object
    integer(shortInt), intent(in)           :: matIdx, what
    class(CENeutron), pointer               :: CENeutronPtr
    real(defReal)                           :: energy, xs
    character(*), parameter                 :: HERE = 'getTrackingXS (ceNeutronDatabase_inter.f90)'

    ! Check dynamic type of physical particle then process request.
    CENeutronPtr => castCENeutronPtr(object)
    energy = CENeutronPtr % getEnergy()
    select case(what)
      case (MATERIAL_XS)
        xs = self % getTrackMatXS(CENeutronPtr, matIdx)

      case (MAJORANT_XS)
        xs = self % getMajorantXS(CENeutronPtr)

      case (TRACKING_XS)
        ! READ ONLY - read from previously updated cache
        if (energy == trackingCache(1) % E) then
          xs = trackingCache(1) % xs
          return

        else
          call fatalError(HERE, 'Failed to update cache during tracking.')

        end if

      case default
        call fatalError(HERE, 'Neither material xs nor majorant xs was asked')

    end select

    ! Update Cache
    trackingCache(1) % E = energy
    trackingCache(1) % xs = xs

  end function getTrackingXS

  !!
  !! Return tracking XS for matIdx
  !!
  !! This is the regular material total cross section unless TMS is used.
  !! If TMS is used, this is the material temperature majorant cross section.
  !!
  !! See nuclearDatabase_inter for details!
  !!
  !! Error:
  !!   fatalError if particle is not CE Neutron
  !!
  function getTrackMatXS(self, object, matIdx) result(xs)
    class(ceNeutronDatabase), intent(inout) :: self
    class(transportObject), intent(in)      :: object
    integer(shortInt), intent(in)           :: matIdx
    real(defReal)                           :: densityFactor, energy, kT, xs
    type(CENeutron), pointer                :: CENeutronPtr
    type(coordList), pointer                :: coordListPtr

    ! Check dynamic type of the particle
    CENeutronPtr => castCENeutronPtr(object)
    coordListPtr => CENeutronPtr % getCoordsPtr()

    ! Check Cache and update if needed
    associate(matCache => materialCache(matIdx))
      densityFactor = getScalarFieldValue(nameDensity, ONE, coordListPtr)
      energy = CENeutronPtr % getEnergy()
      kT = getScalarFieldValue(nameTemperature, self % getMaterial_kT(matIdx), coordListPtr, kBoltzmann_MeV)
      if (any([matCache % densityFactor_track, matCache % E_track, matCache % kT_track] /= [densityFactor, energy, kT])) &
      call self % updateTrackMatXS(densityFactor, energy, kT, matIdx, CENeutronPtr % getRNGPtr())

      ! Return Cross-Section
      xs = matCache % trackXS

    end associate

  end function getTrackMatXS

end module ceNeutronDatabase_inter