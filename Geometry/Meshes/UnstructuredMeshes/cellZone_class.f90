module cellZone_class
  
  use numPrecision
  
  implicit none
  private
  
  !!
  !! Group of elements of an OpenFOAM mesh. Used to assign localIds during mesh importation.
  !!
  type, public :: cellZone
    private
    integer(shortInt)         :: startElementIdx = 0, endElementIdx = 0
    character(:), allocatable :: name
  contains
    procedure                 :: containsElement
    procedure                 :: init
    procedure                 :: getName
    procedure                 :: kill
  end type cellZone

contains
  
  !! Function 'containsElement'
  !!
  !! Basic description:
  !!   Check whether a given element belongs to the cell zone.
  !!
  !! Arguments:
  !!   elementIdx [in] -> Index of the element.
  !!
  !! Result:
  !!   doesIt          -> A logical which is true if the cell zone contains the element.
  !!
  elemental function containsElement(self, elementIdx) result(doesIt)
    class(cellZone), intent(in)   :: self
    integer(shortInt), intent(in) :: elementIdx
    logical(defBool)              :: doesIt
    
    doesIt = self % startElementIdx <= elementIdx .and. elementIdx <= self % endElementIdx
  end function containsElement
  
  !! Subroutine 'init'
  !!
  !! Basic description:
  !!   Sets the first and last element indices in the cell zone.
  !!
  !! Arguments:
  !!   name [in]            -> Name of the cell zone.
  !!   startElementIdx [in] -> Index of the first element in the cell zone.
  !!   endElementIdx [in]   -> Index of the last element in the cell zone.
  !!
  elemental subroutine init(self, name, startElementIdx, endElementIdx)
    class(cellZone), intent(inout) :: self
    character(*), intent(in)       :: name
    integer(shortInt), intent(in)  :: startElementIdx, endElementIdx
    
    self % name = name
    self % startElementIdx = startElementIdx
    self % endElementIdx = endElementIdx
  end subroutine init
  
  !! Function 'getName'
  !!
  !! Basic description:
  !!   Returns the name of the cell zone.
  !!
  !! Result:
  !!   name -> A character string corresponding to the name of the cell zone.
  !!
  elemental function getName(self) result(name)
    class(cellZone), intent(in) :: self
    character(len(self % name)) :: name
    
    name = self % name
  end function getName
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns the cell zone to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(cellZone), intent(inout) :: self
    
    self % startElementIdx = 0
    self % endElementIdx = 0
    if (allocated(self % name)) deallocate(self % name)
  end subroutine kill
end module cellZone_class