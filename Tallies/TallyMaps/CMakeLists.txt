# Add Source Files to the global list
add_subdirectory(Maps1D)
add_sources(
            ./tallyMap_inter.f90
            ./tallyMapFactory_func.f90
            ./multiMap_class.f90
            ./cylindricalMap_class.f90            
           )

add_unit_tests(
               ./Tests/multiMap_test.f90
               ./Tests/cylindricalMap_test.f90
              )
