add_subdirectory(CompoundSurfaces)
add_subdirectory(QuadSurfaces)

add_sources(./surface_inter.f90
            ./surfaceFactory_func.f90
            ./surfaceShelf_class.f90
           )

add_unit_tests( 
                ./Tests/surfaceShelf_test.f90
              )
