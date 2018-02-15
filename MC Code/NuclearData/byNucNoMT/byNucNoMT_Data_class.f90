module byNucNoMT_Data_class

  use numPrecision
  use genericProcedures, only : fatalError, openToRead, removeDuplicates, linFind, &
                                findDuplicates, arrayConcat
  use aceNoMT_class,     only : aceNoMT

  implicit none
  private

  type, public :: byNucNoMT_Data
    !private
    ! Material Data
    character(matNameLen),dimension(:),pointer :: matNames    => null()
    integer(shortInt),dimension(:),pointer     :: matNumNuc   => null()
    integer(shortInt),dimension(:,:),pointer   :: matNucIdx   => null()
    character(ZZidLen),dimension(:,:),pointer  :: matNucNames => null()
    real(defReal),dimension(:,:),pointer       :: matNucDens  => null()
    real(defReal),dimension(:),pointer         :: matTemp     => null()
    ! Isotope Data
    character(zzIdLen),dimension(:),pointer    :: nucNames    => null()
    type(aceNoMT),dimension(:), pointer        :: nucXsData   => null()

  contains
    procedure :: readFrom
    procedure :: print

    procedure,private :: createMatArrays
    procedure,private :: readMaterials
    procedure,private :: createNuclideList
    procedure,private :: assignNucIndices
    procedure,private :: readNuclides

  end type byNucNoMT_Data

contains

  !!
  !! Reads materials and isotopes from provided material and ACE library input files
  !!
  subroutine readFrom(self,matInput,nuclideLib)
    class(byNucNoMT_Data), intent(inout) :: self
    character(*), intent(in)             :: matInput
    character(*), intent(in)             :: nuclideLib

    call self % readMaterials(matInput)

    call self % createNuclideList()

    call self % assignNucIndices()
    call self % readNuclides(nuclideLib)

  end subroutine readFrom

  !!
  !! Read material data from the  input file
  !!
  subroutine readMaterials(self, inputFile)
    class(byNucNoMT_Data), intent(inout) :: self
    character(*), intent(in)         :: inputFile

    integer(shortInt),parameter      :: input=66
    character(99),parameter          :: here='readFrom in byNucNoMT_Data_class.f03'
    character(99)                    :: readMsg

    character(matNameLen)            :: matName
    character(ZZidLen)               :: ZZid
    integer(shortInt)                :: numMat, numNuc, maxNuc
    real(defReal)                    :: temp, numDen
    integer(shortInt)                :: i, j, readStat

    call openToRead(Input,InputFile)

    call readMaxNumNuc(Input,maxNuc)

    ! Read Number of Materials in input file
    read(unit = Input, fmt=* , iostat = readStat, iomsg = readMsg) numMat
    if (readStat > 0) call fatalError(Here, readMsg)

    call self % createMatArrays(numMat,maxNuc)

    ! Read Materials
    do i=1,numMat

      read(unit = Input,    &
           fmt=*,           &
           iostat=readStat, &
           iomsg = readMsg) matName, &
                            temp,    &
                            numNuc

      if (readStat > 0) call fatalError(here, readMsg)
      if (numNuc <= 0)  call fatalError(here,'Number of Material Isotopes is 0 or -ve')

      self % matNames(i)  = matName
      self % matTemp(i)   = temp
      self % matNumNuc(i) = numNuc

      do j=1,numNuc

        read(unit = Input,      &
             fmt=*,             &
             iostat = readStat, &
             iomsg = readMsg)   ZZid ,&
                                numDen

        if (readStat > 0) call fatalError(here, readMsg)
        if (numDen <= 0 ) call fatalError(here,'Numerical Density is -ve')

        self % matNucNames(i,j) = ZZid
        self % matNucDens(i,j)  = numDen
      end do
    end do

    ! Give error if there are more enteries in the input
    read (unit = Input,    &
          fmt =*,          &
          iostat = readStat) ZZid, numDen
    if (readStat /= endOfFile) call fatalError(here,'End of file was not reached, check the ' // &
                                                    'number of isotopes in the last material')

    close(Input)

  end subroutine readMaterials


  !!
  !! Remove repetitions from material definitions and create a list of all used nuclide data
  !! cards without any repetitions
  !!
  subroutine createNuclideList(self)
    class(byNucNoMT_Data),intent(inout)           :: self
    integer(shortInt)                             :: maxNucNames
    character(zzIdLen),dimension(:),allocatable   :: withRepetition
    character(zzIdLen),dimension(:),allocatable   :: noRepetition
    integer(shortInt)                             :: i,j

    maxNucNames=sum(self % matNumNuc)

    ! Crate array of isotope names with repetitions
    allocate(withRepetition(maxNucNames))
    j=1
    do i = 1,size(self % matNames)
      ! Load material names for material i
      withRepetition(j:j+self % matNumNuc(i)) = self % matNucNames(i,1:self % matNumNuc(i))
      j = j + self % matNumNuc(i)
    end do

    noRepetition = removeDuplicates(withRepetition)
    allocate(self % nucNames( size(noRepetition) ))
    self % nucNames = noRepetition

  end subroutine createNuclideList


  !!
  !! Assign nuclide indexes to material nuclides on the unified nuclide list.
  !!
  subroutine assignNucIndices(self)
    class(byNucNoMT_Data),intent(inout)  :: self
    integer(shortInt)                    :: i,j

    do i = 1,size(self % matNames)
      do j = 1,self % matNumNuc(i)
        self % matNucIdx(i,j) = linFind(self % nucNames, self % matNucNames(i,j))
        if (self % matNucIdx(i,j) == -1 ) then
          call fatalError('assignIsoIndices (byNucNoMT_class.f90)', &
                          'Isotope ' // self % matNucNames(i,j) //' was not found')
        end if
      end do
    end do
  end subroutine assignNucIndices

  !!
  !! Read library of ACE nuclide cards and read nuclide data
  !!
  subroutine readNuclides(self,libraryPath)
    class(byNucNoMT_Data), intent(inout)         :: self
    character(*),intent(in)                      :: libraryPath
    integer(shortInt),parameter                  :: library=78
    character(99)                                :: readMsg
    character(zzIdLen),dimension(:),allocatable  :: zzIDs
    integer(shortInt),dimension(:),allocatable   :: startLine
    character(pathLen),dimension(:),allocatable  :: nucPath
    integer(shortInt)                            :: i, j, readStat
    integer(shortInt)                            :: libLen
    character(zzIDLen)                           :: currentZZId, &
                                                    pastZZId
    call openToRead(library,libraryPath)

    ! Find length of isotope library
    ! There is a discrepancy in behaviour on Linux and Windows-Cygwin system.
    ! On Windows after readeing last line one addoitional read is made in which last line is read
    ! again but with iostat=-1. On lunux iostat = -1 is returned when READING LAST LINE THE FIRST TIME.
    ! Thus if code that runs porperly on linux is reused on windows last line is read twice!. Explicit
    ! check whether that happens is required to obtain correct behaviour on both systems.
    libLen=0
    do
      read(unit = library, fmt=*,iostat=readStat,iomsg = readMsg) currentZZId
      if(readStat == -1 .and. trim(currentZZId)==trim(pastZZId)) exit ! See character equality to indicate double read of last line
      libLen = libLen + 1
      pastZZId = currentZZId
    end do
    rewind(library)

    ! Allocate and read library
    allocate(zzIDs(libLen))
    allocate(startLine(libLen))
    allocate(nucPath(libLen))

    do i=1,libLen
      read(library,"(A10, I12, A100)" ) zzIds(i), startLine(i), nucPath(i)
    end do

    ! Check library for repeted zzIDs
    if (size(zzIds) /= size(removeDuplicates(zzIds))) then
      call fatalError('readIsotopes (byNucNoMT_Data_class.f90)', &
                      'Duplicate zzIds found in ACE library file: ' //&
                      arrayConcat(findDuplicates(zzIds)))
    end if

    ! Allocate Memory for isotopic data
    allocate(self % nucXsData(size(self % nucNames)))

    ! Read Isotope Data
    do i=1,size(self % nucNames)
      ! **** This Search need to be modernised ****!
      j = linFind(zzIds,self % nucNames(i))
      if (j == -1) then
        call fatalError('readIsotopes (byNucNoMT_Data_class.f90)', &
                        'Isotope ' // self % nucNames(i) //' was not found')
      end if

      print *, "Reading : ", zzIds(j), startLine(j), nucPath(j)
      call self % nucXsData(i) % init(nucPath(j),startLine(j))

    end do
    print *, size(self % nucXSData)


    close(library)

  end subroutine readNuclides

  !!
  !! Allocate space for all arrays that contain data about materials
  !!
  subroutine createMatArrays(self,numMat,maxNuc)
    class(byNucNoMT_Data), intent(inout)  :: self
    integer(shortInt),intent(in)          :: numMat, maxNuc

    allocate(self % matNames(numMat))
    allocate(self % matNumNuc(numMat))
    allocate(self % matNucIdx(numMat,maxNuc))
    allocate(self % matNucNames(numMat,maxNuc))
    allocate(self % matNucDens(numMat,maxNuc))
    allocate(self % matTemp(numMat))

  end subroutine createMatArrays


  !!
  !! Read input file and indentify maximum number of nuclides in the problem
  !!
  subroutine readMaxNumNuc(Input,maxNuc)
    integer(shortInt),intent(in)      :: Input
    integer(shortInt),intent(out)     :: maxNuc
    character(99),parameter           :: Here='readMaxNumNuc (byNucNoMT_Data_class.f90)'
    character(99)                     :: readMsg
    integer(shortInt)                 :: numMat, numNuc, readStat, i, j
    character(3)                      :: dummyChar
    real(defReal)                     :: dummyReal

    rewind(Input)

    read(unit = Input, fmt=* , iostat = readStat, iomsg = readMsg) numMat
    if (readStat > 0) call fatalError(Here, readMsg)

    maxNuc = -1
    do i=1,numMat
      read(unit = Input, fmt=*, iostat=readStat, iomsg = readMsg) dummyChar, &
                                                                  dummyReal,    &
                                                                  numNuc
      if (readStat > 0) call fatalError(Here, readMsg)
      if (numNuc <= 0)  call fatalError(Here,'Number of Material Nuclides is 0 or -ve')
      maxNuc=max(maxNuc,numNuc)
      ! Skip lines to get to next material header
      do j=1,numNuc
        read(Input,*) dummyChar, dummyReal
      end do
    end do

    rewind(Input)

  end subroutine readMaxNumNuc


  !!
  !! Print data about all materials to the console
  !!
  subroutine print(self)
    class(byNucNoMT_Data), intent(in) :: self
    character(99)                     :: format
    integer(shortInt)                 :: i

    print '(a)', 'Material Names:'
    print '(a)', self % matNames

    print '(a)', 'Nuclide Numbers:'
    print '(i5)', self % matNumNuc

    print '(a)', 'Materials Temperatures:'
    print '(f10.3)', self % matTemp

    print '(a)', 'Nuclide Names:'
    do i=1,size(self % matNucNames,1)
      print '(9999a)', self % matNucNames(i,:)
    end do

    print '(a)', 'Nuclide Indexes:'
    do i=1,size(self % matNucNames,1)
      print '(9999I5)', self % matNucIdx(i,:)
    end do

    print '(a)', 'Nuclide Densities:'
    do i=1,size(self % matNucNames,1)
      print '(9999es15.5)', self % matNucDens(i,:)
    end do

    print '(a)', 'Nuclide Names Array'
    print '(a)', self % nucNames

  end subroutine print
    
end module byNucNoMT_Data_class