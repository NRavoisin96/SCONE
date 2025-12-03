!!
!! Material Menu is a module (Singleton) that contains global definitions of different materials
!!
!! It exists to make it easier for all databases to refer to the same materials by the same
!! name and index. This is necessary to avoid confusion resulting from different materials with the
!! same name or index in different databases.
!!
!! Public Members:
!!   materialDefs -> array of material definitions of type materialItem
!!   nameMap      -> Map that maps material name to matIdx
!!   colourMap    -> Map that maps matIdx to 24bit colour (to use for visualisation)
!!
!! Interface:
!!   init      -> Load material definitions from a dictionary
!!   kill      -> Return to uninitialised state
!!   display   -> Display information about all defined materials to console
!!   getMatPtr -> Return pointer to a detailed material information (materialItem)
!!   nMat      -> Return number of materials
!!   matName   -> Return material Name given Index
!!   matIdx    -> Return material Index given Name
!!
module materialMenu_mod

  use atomicDensitiesCalculator_inter,       only : atomicDensitiesCalculator
  use atomicDensitiesCalculatorFactory_func, only : new_atomicDensitiesCalculator
  use charMap_class,                         only : charMap
  use colours_func,                          only : rgb24bit
  use dictionary_class,                      only : dictionary
  use errors_mod,                            only : fatalError
  use genericProcedures,                     only : charToInt, numToChar
  use intMap_class,                          only : intMap
  use nuclideInfo_class,                     only : buildNuclideInfoPayload, nuclideInfo
  use numPrecision
  use physicalPropertyLaw_inter,             only : physicalPropertyLaw
  use physicalPropertyLawFactory_func,       only : new_physicalPropertyLaw
  use universalVariables,                    only : NOT_FOUND, OUTSIDE_MAT, UNDEF_MAT, VOID_MAT

  implicit none
  private

  !!
  !! This is a type which collects all information about a single material definition
  !!
  !! Public Members:
  !!   name      -> name of material
  !!   matIdx    -> material index of the material
  !!   T         -> material temperature [K]
  !!   dens      -> vector of densities [1/barn/cm]
  !!   nuclides  -> associated vector of nuclide types
  !!   extraInfo -> dictionary with extra keywords
  !!
  !! Interface:
  !!   init -> build material item from dictionary
  !!   kill -> return to uninitialised state
  !!
  !! Sample Input Dictionary:
  !!
  !!   matDef {
  !!     temp 273;
  !!     #moder {1001.03 (h-h2o.43);}#
  !!     #tms 1;#
  !!     composition {
  !!       1001.03  5.028E-02;
  !!       8016.03  2.505E-02;
  !!       5010.03  2.0E-005;
  !!     }
  !!     xsFile /home/uberMoffTarkin/XS/mat1.xs;
  !!     #rgb (255 0 0); # // RGB colour to be used in visualisation
  !!   }
  !!
  !! Sample with stochastic mixing:
  !!   matDef {
  !!     temp 300;
  !!     moder {1001.03 (h-h2o.43 h-h2o.53);}
  !!     composition {
  !!       1001.03  5.028E-02;
  !!       8016.03  2.505E-02;
  !!       5010.03  2.0E-005;
  !!     }
  !!   }
  !!
  !! NOTE: the moder dictionary is optional, necessary only if S(a,b) thermal scattering
  !!       data are used. If some nuclides are included in moder but not in composition,
  !!       an error is raised.
  !!       Including two entries in moder will invoke stochastic mixing, i.e.,
  !!       stochastic interpolation between the two data libraries.
  !!
  type, public :: materialItem
    character(:), allocatable                     :: name
    class(atomicDensitiesCalculator), allocatable :: densitiesCalculator
    class(physicalPropertyLaw), allocatable       :: densityLaw, thermalConductivityLaw
    integer(shortInt)                             :: matIdx = 0
    logical(defBool)                              :: hasTMS = .false.
    real(defReal)                                 :: density = ZERO, inverseDensity = ZERO, T = ZERO
    type(dictionary)                              :: extraInfo
    type(nuclideInfo), dimension(:), allocatable  :: nuclides
  contains
    procedure :: getAtomicDensities
    procedure :: getAtomicDensity
    procedure :: init => init_materialItem
    procedure :: initNuclide
    procedure :: kill => kill_materialItem
    procedure :: display => display_materialItem
  end type materialItem

  !! Parameters
  integer(shortInt), parameter :: OUTSIDE_COLOUR = int(z'ffffff', shortInt), UNDEFINED_COLOUR = int(z'00ff00', shortInt), &
                                  VOID_COLOUR = int(z'000000', shortInt)


  !! MODULE COMPONENTS
  type(materialItem), dimension(:), allocatable, target, public :: materialDefs
  type(charMap), target, public                                 :: nameMap
  type(intMap), public                                          :: colourMap

  public :: display, getMatPtr, init, kill, matIdx, matName, nMat

contains

  !!
  !! Initialises materialMenu from a dictionary with material definitions
  !!
  !! Args:
  !!   dict [in] -> dictionary with material definitions
  !!
  !! Errors:
  !!   None from Here
  !!
  subroutine init(dict)
    class(dictionary), intent(in)                :: dict
    character(nameLen), dimension(:), allocatable :: matNames
    integer(shortInt)                           :: i
    character(nameLen)                          :: temp

    ! Clean whatever may be already present
    call kill()

    ! Load all material names
    call dict % keys(matNames,'dict')

    ! Allocate space
    allocate(materialDefs(size(matNames)))

    ! Load definitions
    do i = 1, size(matNames)
      call materialDefs(i) % init(matNames(i), i, dict % getDictPtr(matNames(i)))
      call nameMap % add(matNames(i), i)

    end do

    ! Add special Material keywords to the dictionary
    temp = 'void'
    call nameMap % add(temp, VOID_MAT)
    temp = 'outside'
    call nameMap % add(temp, OUTSIDE_MAT)

    !! Load colours for the special materials
    call colourMap % add(VOID_MAT, VOID_COLOUR)
    call colourMap % add(OUTSIDE_MAT, OUTSIDE_COLOUR)
    call colourMap % add(UNDEF_MAT, UNDEFINED_COLOUR)

  end subroutine init


  !!
  !! Returns material Menu to an uninitialised state
  !!
  subroutine kill()
    integer(shortInt) :: i

    call nameMap % kill()
    call colourMap % kill()
    if (allocated(materialDefs)) then
      do i = 1, size(materialDefs)
        call materialDefs(i) % kill()

      end do
      deallocate(materialDefs)

    end if

  end subroutine kill

  !!
  !! Print material definition information to the console
  !!
  !! Args:
  !!   None
  !! Errors:
  !!   None
  !!
  subroutine display()
    integer(shortInt) :: i

    print '(A60)', repeat('<>',30)
    print '(A)', "^^ MATERIAL DEFINITIONS ^^"

    do i = 1,size(materialDefs)
      call materialDefs(i) % display()
      ! Print separation line
      print '(A)', " ><((((*>  +  <*))))><"
    end do

    print '(A60)', repeat('<>',30)

  end subroutine display

  !!
  !! Return Material Name given index
  !!
  !! Args:
  !!   idx [in] -> Material Index
  !!
  !! Result:
  !!   nameLen long character with material name
  !!
  !! Error:
  !!   If idx is -ve or larger then number of defined materials
  !!   Empty string '' is returned as its name
  !!
  function matName(idx) result(name)
    integer(shortInt), intent(in) :: idx
    character(nameLen)            :: name

    if (idx <= 0 .or. nMat() < idx) then
      name = ''

    else
      name = materialDefs(idx) % name
    end if

  end function matName

  !!
  !! Return material index Given Name
  !!
  !! Args:
  !!   name [in] -> material name
  !!
  !! Result:
  !!   matIdx corresponding to name
  !!
  !! Error:
  !!   If name does not correspond to any defined material NOT_FOUND is returned
  !!
  function matIdx(name) result(idx)
    character(*), intent(in) :: name
    integer(shortInt)        :: idx

    idx = nameMap % getOrDefault(name, NOT_FOUND)

  end function matIdx

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! TYPE PROCEDURES
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !!
  !!
  !!
  pure function getAtomicDensities(self) result(atomicDensities)
    class(materialItem), intent(in)          :: self
    real(defReal), dimension(:), allocatable :: atomicDensities

    if (allocated(self % nuclides)) then
      allocate(atomicDensities(size(self % nuclides)))
      atomicDensities = self % nuclides % getDensity()

    else
      allocate(atomicDensities(0))

    end if

  end function getAtomicDensities

  !!
  !!
  !!
  elemental function getAtomicDensity(self, idx) result(atomicDensity)
    class(materialItem), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    real(defReal)                   :: atomicDensity

    atomicDensity = self % nuclides(idx) % getDensity()

  end function getAtomicDensity

  !!
  !! Initialise material definition from a dictionary and name
  !!
  !! Args:
  !!   name [in] -> character with material name
  !!   idx  [in] -> material index
  !!   dict [in] -> dictionary with material definition
  !!
  !! Errors:
  !!   FatalError if dictionary does not contain valid material definition.
  !!
  subroutine init_materialItem(self, name, idx, dict)
    class(materialItem), intent(inout)            :: self
    character(nameLen), intent(in)                :: name
    integer(shortInt), intent(in)                 :: idx
    type(dictionary), intent(in)                  :: dict
    character(nameLen), dimension(:), allocatable :: keys, moderKeys
    integer(shortInt)                             :: i, nNuclides, nSab, foundModer
    integer(shortInt), dimension(:), allocatable  :: temp
    real(defReal), dimension(:), allocatable      :: densitiesCalculatorInputArray
    type(dictionary), pointer                     :: compDict, moderDict, nuclidesDict
    character(*), parameter                       :: HERE = 'init_materialItem (materialMenu_mod.f90)'

    ! Return to initial state
    call self % kill()

    ! Load easy components properties
    self % name = trim(name)
    self % matIdx = idx

    ! Check TMS flag and read temperature
    call dict % getOrDefault(self % hasTMS, 'tms', .false.)
    call dict % getOrDefault(self % T, 'temp', ZERO)
    if (self % T < ZERO) &
    call fatalError(HERE, 'Temperature of material: '//self % name//' is negative: '//numToChar(self % T)//'K.')

    ! Build density law and compute density.
    if (dict % isPresent('densityLaw')) then
      call new_physicalPropertyLaw(dict % getDictPtr('densityLaw'), self % densityLaw)
      self % density = self % densityLaw % computeProperty(self % T)
      self % inverseDensity = ONE / self % density

    end if

    ! Get composition dictionary and load composition.
    if (.not. dict % isPresent('composition')) &
    call fatalError(HERE, "Missing 'composition' subdictionary for material: "//self % name//'.')
    compDict => dict % getDictPtr('composition')

    ! Allocate space for nuclide information.
    nNuclides = 0
    if (compDict % isPresent('nuclides')) then
      nuclidesDict => compDict % getDictPtr('nuclides')
      call nuclidesDict % keys(keys)
      nNuclides = size(keys)
      allocate(densitiesCalculatorInputArray(nNuclides), self % nuclides(nNuclides))
      call new_atomicDensitiesCalculator(compDict, self % densitiesCalculator)

    end if

    ! Check if S(a,b) files are specified.
    nSab = 0
    if (dict % isPresent('moder')) then
      moderDict => dict % getDictPtr('moder')
      call moderDict % keys(moderKeys)
      nSab = size(moderKeys)
      
    end if

    ! Build nuclides.
    foundModer = 0
    if (0 < nNuclides) then
      do i = 1, nNuclides
        ! Retrieve input density for current nuclide and initialise current nuclide.
        call nuclidesDict % get(densitiesCalculatorInputArray(i), keys(i))
        call self % initNuclide(keys(i), nSab, moderDict, foundModer, self % nuclides(i))

      end do
      call self % densitiesCalculator % computeAtomicDensities(self % density, densitiesCalculatorInputArray, self % nuclides)

    end if

    ! Make sure if a moderator is provided the nuclide is present in the composition
    if (foundModer /= nSab) then
      print *, moderKeys
      call fatalError(HERE, 'Nuclides requested for S(alpha, beta) are not present in composition. '// &
              numToChar(nSab)//' nuclides requested but '//numToChar(foundModer)//' nuclides found.')
    end if

    ! Add colour info if present
    if (dict % isPresent('rgb')) then
      call dict % get(temp, 'rgb')
      if (size(temp) /= 3) call fatalError(HERE, "'rgb' keyword must have 3 values.")
      call colourMap % add(idx, rgb24bit(temp(1), temp(2), temp(3)))

    end if

    ! Save dictionary
    self % extraInfo = dict

    ! TODO: Remove composition subdictionary from extraInfo
    !       Or rather do not copy it in the first place

  end subroutine init_materialItem

  !!
  !!
  !!
  subroutine initNuclide(self, key, nSab, moderDict, nModerators, nuclide)
    class(materialItem), intent(in)               :: self
    character(nameLen), intent(in)                :: key
    integer(shortInt), intent(in)                 :: nSab
    type(dictionary), intent(in)                  :: moderDict
    integer(shortInt), intent(inout)              :: nModerators
    type(nuclideInfo), intent(out)                :: nuclide
    character(nameLen), dimension(:), allocatable :: filenames
    integer(shortInt)                             :: dot, nFiles
    logical(defBool)                              :: flag
    type(buildNuclideInfoPayload)                 :: payload
    character(*), parameter                       :: HERE = 'initNuclide (materialMenu_mod.f90)'

    ! First ensure nuclide definition is valid.
    if (.not. isNucDefinition(key)) call fatalError(HERE, 'Input is not ZZZAAA.TT formated definition: '//trim(key)//'.')

    ! Assemble payload. Find location of the dot and catch leading zeros in ZA id.
    dot = scan(key, '.')
    if (key(1:1) == '0') call fatalError(HERE, 'ZA id begins with 0.')
    
    payload % atomicNumber = charToInt(key(1:dot - 4), error = flag)
    payload % massNumber = charToInt(key(dot - 3:dot - 1), error = flag)
    payload % evaluationNumber = charToInt(key(dot + 1:len_trim(key)), error = flag)
    if (flag) call fatalError(HERE, 'Failed to convert: '//trim(key)//' into nuclide information.')

    ! Now check if S(alpha, beta) is on and required for that nuclide.
    if (0 < nSab .and. moderDict % isPresent(key)) then
      payload % hasSab = .true.
      nModerators = nModerators + 1

      ! Check for stochastic mixing - this will depend on the size of the array of files produced.
      call moderDict % get(filenames, key)
      nFiles = size(filenames)
      select case(nFiles)
        case(1, 2)
          payload % sabFiles = filenames
          if (nFiles == 2) payload % sabMix = .true.

        case default
          print *, filenames
          call fatalError(HERE, 'Unexpectedly long moder contents. Should be 1 or 2 entries.')

      end select

    end if
    call nuclide % init(payload)

    contains
      !!
      !! Helper function to identify nuclide definition string
      !!
      !! Nuclide definition string has a following format:
      !!   ZZZAAA.TT
      !!
      !! ZZZ -> Up to 3 Digits   [0-9] that specify Atomic Number (minimum 1)
      !! AAA -> EXACTLY 3 Digits [0-9] that specify Mass Number
      !! TT  -> EXACTLY 2 Digits [0-9] that specify Evaluation Number
      !!
      !! Must also be left-adjusted and padded only with spaces.
      !!
      !! Args:
      !!   str [in] -> character string that may or may not contain Nuclide definition
      !!
      !! Result:
      !!   True if str matches the Nuclide Definition format. False otherwise
      !!
      function isNucDefinition(str) result(isIt)
        character(nameLen), intent(in) :: str
        character(:), allocatable      :: subStr
        logical(defBool)               :: isIt
        integer(shortInt)              :: L
        character(*), parameter        :: SET = '0123456789'

        ! Initialise isIt = .false.
        isIt = .false.

        ! Save trim length of the string
        L = len_trim(str)

        ! Check that length is as expected
        if (L < 7 .or. 9 < L) return

        ! Verify that the location of the dot is consistent from number of digits before and after dot.
        subStr = str(1:L)
        isIt = all([verify(subStr, SET), verify(subStr, SET, back = .true.)] == scan(subStr, '.'))

      end function isNucDefinition

  end subroutine initNuclide

  !!
  !! Return material Item to uninitialised state
  !!
  subroutine kill_materialItem(self)
    class(materialItem), intent(inout) :: self
    integer(shortInt)                  :: i

    ! Return static components to default
    self % matIdx = 0
    self % density = ZERO
    self % inverseDensity = ZERO
    self % T = ZERO

    ! Deallocate allocatable components
    if (allocated(self % name)) deallocate(self % name)
    if (allocated(self % densitiesCalculator)) deallocate(self % densitiesCalculator)
    if (allocated(self % densityLaw)) then
      call self % densityLaw % kill()
      deallocate(self % densityLaw)

    end if
    if (allocated(self % nuclides)) then
      do i = 1, size(self % nuclides)
        call self % nuclides(i) % kill()

      end do
      deallocate(self % nuclides)

    end if
    call self % extraInfo % kill()

  end subroutine kill_materialItem

  !!
  !! Prints the definition of material to the console
  !! Uses up to 60 columns
  !!
  !! Args:
  !!   None
  !! Errors:
  !!   None
  !!
  subroutine display_materialItem(self)
    class(materialItem), intent(in) :: self
    integer(shortInt)               :: i

    print '(A)', 'Material: '// trim(self % name) //' with index: ' // numToChar(self % matIdx)
    print '(A)', 'Temperature [K]: '//numToChar(self % T)
    print '(A)', 'Nuclide Composition:'
    print '(3A13, A20)', 'Atomic #', 'Mass #', 'Evaluation #', 'Density [1/barn/cm]'

    do i = 1, size(self % nuclides)
      call self % nuclides(i) % display('(3I13, ES20.10)')

    end do

  end subroutine display_materialItem

  !!
  !! Get pointer to a material definition under matIdx
  !!
  !! Args:
  !!   idx [in] -> Index of the material
  !!
  !! Result:
  !!   Pointer to a materialItem with the definition
  !!
  !! Errors:
  !!   FatalError if idx does not correspond to any defined material
  !!   FatalError if material definitions were not loaded
  !!
  function getMatPtr(idx) result(ptr)
    integer(shortInt), intent(in) :: idx
    type(materialItem), pointer   :: ptr
    character(*), parameter       :: HERE = 'getMatPtr (materialMenu_mod.f90)'

    ! Check if materialMenu is initialised
    if (.not. allocated(materialDefs)) call fatalError(HERE, "Material definitions were not loaded.")

    ! Verify matIdx
    if (idx < 1 .or. nMat() < idx) call fatalError(HERE, "matIdx: "//numToChar(idx)//&
                                                   " does not correspond to any material.")

    ! Attach pointer.
    ptr => materialDefs(idx)

  end function getMatPtr


  !!
  !! Return number of materials
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Number of defined materials
  !!
  !! Errors:
  !!   Return 0 if materialMenu was not yet loaded
  !!
  elemental function nMat() result(N)
    integer(shortInt) :: N

    N = 0
    if (allocated(materialDefs)) N = size(materialDefs)

  end function nMat

end module materialMenu_mod