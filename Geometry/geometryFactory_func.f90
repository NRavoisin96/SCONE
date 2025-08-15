!!
!! Module to build new geometries and add them to the geometry registry
!!
module geometryFactory_func

  use numPrecision
  use fieldFactory_func,  only : new_field
  use genericProcedures,  only : fatalError
  use dictionary_class,   only : dictionary
  use charMap_class,      only : charMap

  ! Geometry
  use geometryMesh_class, only : geometryMesh
  use geometryStd_class,  only : geometryStd
  use geometryReg_mod,    only : gr_addGeom => addGeom
  use geometry_inter,     only : geometry

  ! Material interface
  use materialMenu_mod,   only : mm_nameMap => nameMap

  implicit none
  private


  !! Parameters
  character(nameLen), dimension(*), parameter :: AVAILABLE_GEOMETRIES = ['geometryMesh', &
                                                                         'geometryStd ']

  ! Public interface
  public :: new_geometry

contains

  !!
  !! Initialises a geometry from dictionary. This is then allocated inside
  !! the geometry registry
  !!
  !! Args:
  !!   dict [in]  -> Dictionary with geometry definition
  !!   name [in]  -> Name of the geometry for the geometry registry
  !!   silent [in] -> Optional. Set to .true. to surpress console messeges. Default .false.
  !!     Note that errors will still be printed with silent=.true.
  !!
  !! Errors:
  !!   fatalError is type of geometry is unknown
  !!
  subroutine new_geometry(dict, name, silent)
    class(dictionary), intent(in)                 :: dict
    character(nameLen), intent(in)                :: name
    logical(defBool), optional, intent(in)        :: silent
    class(dictionary), pointer                    :: fieldsDict
    class(geometry), allocatable                  :: geom
    character(nameLen)                            :: type
    character(nameLen), dimension(:), allocatable :: fieldNames
    integer(shortInt)                             :: i
    logical(defBool)                              :: silent_l
    character(*), parameter                       :: Here = 'new_geometry (geometryFactory_func.f90)'

    ! Get silent flag
    silent_l = .false.
    if (present(silent)) silent_l = silent

    ! Get type
    call dict % get(type, 'type')

    ! Allocate to right type
    select case (type)
      case('geometryMesh')
        allocate(geometryMesh :: geom)

      case('geometryStd')
        allocate(geometryStd :: geom)

      case default
        print '(A)', 'AVAILABLE GEOMETRIES'
        print '(A)', AVAILABLE_GEOMETRIES
        call fatalError(Here, trim(type)// ' is not valid geometry. See list above.')

    end select

    ! Initialise geometry
    call geom % init(dict, mm_nameMap, silent_l)

    ! Call geometry registry to add geometry
    call gr_addGeom(geom, name)

    ! Build fields associated with the current geometry.
    if (dict % isPresent('fields')) then
      fieldsDict => dict % getDictPtr('fields')
      call fieldsDict % keys(fieldNames, 'dict')
      do i = 1, size(fieldNames)
        call new_field(fieldsDict % getDictPtr(fieldNames(i)), fieldNames(i))

      end do

    end if

  end subroutine new_geometry

end module geometryFactory_func
