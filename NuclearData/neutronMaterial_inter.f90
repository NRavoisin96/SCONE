module neutronMaterial_inter

  use materialHandle_inter,    only : materialHandle
  use neutronXsPackages_class, only : neutronMacroXSs
  use numPrecision
  use transportObject_inter,   only : transportObject

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: neutronMaterial_CptrCast

  !!
  !! Abstract interface far all neutron Materials (CE and MG)
  !!
  !! It was created to expose access to some key information in the context of
  !! tallying where one is not interested whether MG or CE data is used
  !!
  !! Interface:
  !!   materialHandle interface
  !!   isFissle    -> Return true if material is fissile
  !!   getMacroXSs -> Return Macroscopic XSs given particle with energy data
  !!
  type, public, abstract, extends(materialHandle) :: neutronMaterial
    private
  contains
    generic                              :: getMacroXSs => getMacroXSs_byP
    procedure(getMacroXSs_byP), deferred :: getMacroXSs_byP
    procedure(isFissile), deferred       :: isFissile
  end type neutronMaterial

  abstract interface
    !!
    !! Return Macroscopic XSs for the material given particle
    !!
    !! Args:
    !!   xss [out]    -> Cross section package to store the data
    !!   p [in]       -> Particle that provides energy or energy group
    !!
    !! Errors:
    !!   fatalError if energy value/group is outside bounds
    !!   fatalError if MG particle is given to CE data and vice versa
    !!
    subroutine getMacroXSs_byP(self, object, xss)
      import                             :: neutronMacroXSs, neutronMaterial, transportObject
      class(neutronMaterial), intent(in) :: self
      class(transportObject), intent(in) :: object
      type(neutronMacroXSs), intent(out) :: xss
    end subroutine getMacroXSs_byP

    !!
    !! Return .true. if the MG material is fissile
    !!
    !! Args:
    !!   None
    !!
    !! Errors:
    !!   None
    !!
    elemental function isFissile(self) result(isIt)
      import                             :: defBool, neutronMaterial
      class(neutronMaterial), intent(in) :: self
      logical(defBool)                   :: isIt
    end function isFissile

  end interface

contains


  !!
  !! Cast materialHandle pointer to neutronMaterial pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class materialHandle
  !!
  !! Result:
  !!   Null is source is not of neutronMaterial
  !!   Pointer to source if source is neutronMaterial class
  !!
  pure function neutronMaterial_CptrCast(source) result(ptr)
    class(materialHandle), pointer, intent(in) :: source
    class(neutronMaterial), pointer            :: ptr

    select type(source)
      class is(neutronMaterial)
        ptr => source

      class default
        ptr => null()
    end select

  end function neutronMaterial_CptrCast

end module neutronMaterial_inter