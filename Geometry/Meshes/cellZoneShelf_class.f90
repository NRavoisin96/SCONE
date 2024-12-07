module cellZoneShelf_class
  
  use numPrecision
  use universalVariables, only : targetNotFound
  use genericProcedures,  only : openToRead, linFind
  use element_class,      only : element
  use cellZone_class,     only : cellZone
  
  implicit none
  private
  
  !! Storage space for element zones in a given unstructured mesh.
  !!
  !! Public members:
  !!   shelf          -> Array to store element zones.
  !!
  !! Interface:
  !!   findCellZone   -> Returns the cell zone containing a given element.
  !!   getIdxFromName -> Returns the cell zone corresponding to a given name.
  !!   kill           -> Returns to an uninitialised state.
  !!
  type, public :: cellZoneShelf
    private
    type(cellZone), dimension(:), allocatable, public :: shelf
  contains
    procedure                                         :: findCellZone
    procedure                                         :: getIdxFromName
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