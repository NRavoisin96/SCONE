module cellZoneShelf_class
  
  use numPrecision
  use universalVariables, only : targetNotFound
  use genericProcedures,  only : openToRead, linFind
  use element_class,      only : element
  use cellZone_class,     only : cellZone
  
  implicit none
  private
  
  !!
  !! Storage space for cell zones in a given OpenFOAM mesh.
  !!
  !! Public members:
  !!   shelf          -> Array to store cell zones.
  !!
  !! Interface:
  !!   init           -> Initialises from cellZones files.
  !!   findCellZone   -> Returns the cell zone containing a given element.
  !!   getIdxFromName -> Returns the cell zone corresponding to a given name.
  !!   kill        -> Returns to an uninitialised state.
  !!
  type, public :: cellZoneShelf
    private
    type(cellZone), dimension(:), allocatable, public :: shelf
  contains
    procedure                                         :: findCellZone
    procedure                                         :: getIdxFromName
    procedure                                         :: init
    procedure                                         :: kill
  end type cellZoneShelf

contains
  
  !! Function 'findCellZone'
  !!
  !! Basic description:
  !!   Returns the index of the cell zone to which a given element belongs.
  !!
  !! Arguments:
  !!   elementIdx [in] -> Index of the element.
  !!
  !! Result:
  !!   cellZoneIdx -> Index of the cell zone to which the element belongs.
  !!
  elemental function findCellZone(self, elementIdx) result(cellZoneIdx)
    class(cellZoneShelf), intent(in) :: self
    integer(shortInt), intent(in)    :: elementIdx
    integer(shortInt)                :: cellZoneIdx
    
    do cellZoneIdx = 2, size(self % shelf) + 1
      if (self % shelf(cellZoneIdx - 1) % containsElement(elementIdx)) return
    end do
  end function findCellZone
  
  !! Function 'getIdxFromName'
  !!
  !! Basic description:
  !!   Returns the index of the cell zone whose name matches targetName. Returns 'targetNotFound' if
  !!   no cell zone name matches targetName.
  !!
  !! Arguments:
  !!   targetName [in] -> Name of the cell zone.
  !!
  !! Result:
  !!   cellZoneIdx     -> Index of the cell zone whose name matches targetName.
  !!
  elemental function getIdxFromName(self, targetName) result(cellZoneIdx)
    class(cellZoneShelf), intent(in)                  :: self
    character(*), intent(in)                          :: targetName
    integer(shortInt)                                 :: cellZoneIdx, i
    character(nameLen), dimension(size(self % shelf)) :: names
    
    do i = 1, size(self % shelf)
      names(i) = self % shelf(i) % getName()

    end do
    cellZoneIdx = linFind(names, targetName)
  end function getIdxFromName
  
  !! Subroutine 'init'
  !!
  !! Basic description:
  !!   Initialises the shelf from the cellZones file.
  !!
  !! Notes:
  !!   Due to the structure of the cellZones file the subroutine first scans the file a first time
  !!   to obtain the name of each cell zone and then rewinds it to read the indices of the elements
  !!   contained in each zone.
  !!
  !! Arguments:
  !!   folderPath [in] -> Path of the folder containing the OpenFOAM mesh files.
  !!
  subroutine init(self, folderPath)
    class(cellZoneShelf), intent(inout)          :: self
    character(*), intent(in)                     :: folderPath
    integer(shortInt)                            :: i, j, unit = 10, nElements, nCellZones
    integer(shortInt), allocatable, dimension(:) :: elementIndices
    character(100)                               :: string
    character(:), allocatable                    :: name
    logical(defBool)                             :: singleLine

    ! Open the 'cellZones' data file and read it. Skip lines until a blank line is encountered.
    call openToRead(unit, folderPath//'cellZones')
    read(unit, "(a)") string
    do while (len_trim(string) > 0)
      read(unit, "(a)") string
    end do
    
    ! Read the current line and copy the number of cell zones into the variable 'nCellZones'.
    read(unit, "(a)") string
    read(string, *) nCellZones
    
    ! Allocate memory and loop through all cell zones.
    allocate(self % shelf(nCellZones))
    do i = 1, nCellZones
      ! Initialise singleLine = .false. and skip lines until a '{' is encountered.
      singleLine = .false.
      do while (index(string, "{") == 0)
        read(unit, "(a)") string

      end do

      ! Go back to the previous line and read the name of the current cell zone.
      backspace(unit)
      read(unit, "(a)") string
      name = trim(string)
      
      ! Skip lines until the word 'cellLabels' is encoutered.
      do while (index(string, "cellLabels") == 0)
        read(unit, "(a)") string

      end do

      ! If the line contains a ')' then all the elements contained in the current cell zones are
      ! written on the same line.
      if (index(string, ")") > 0) singleLine = .true.

      ! Read the number of elements in the current cell zone depending on whether they are all
      ! written a single line or not.
      if (singleLine) then
        read(string(index(string, ">") + 2:index(string, "(") - 1), *) nElements

      else
        ! Skip one line and retrieve the number of elements in the cell zone..
        read(unit, "(a)") string
        read(string, *) nElements

      end if

      ! Allocate memory.
      allocate(elementIndices(nElements))

      ! Now read the indices of the elements in the current cell zone. Again this depends on
      ! whether they are all written on a single line or not.
      if (singleLine) then
        ! Retrieve the element indices.
        read(string(index(string, "(") + 1:index(string, ")") - 1), *) elementIndices

      else
        ! Skip two lines and retrieve the index of each element in the cell zone line by line..
        do j = 1, 2
          read(unit, "(a)") string

        end do
        do j = 1, nElements
          read(string, *) elementIndices(j)
          read(unit, "(a)") string

        end do

      end if

      ! Set the start and end elements in the cell zone (note that we encrement by one since Fortran
      ! starts indexing at one instead of zero) and reset elementIndices array.
      call self % shelf(i) % init(name, elementIndices(1) + 1, elementIndices(nElements) + 1)
      deallocate(elementIndices)

    end do
    
    ! Close the 'cellZones' file.
    close(unit)
  end subroutine init
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(cellZoneShelf), intent(inout) :: self
    
    if (allocated(self % shelf)) deallocate(self % shelf)

  end subroutine kill
  
end module cellZoneShelf_class