# Add Source Files to a global List
add_subdirectory(BoundingBoxes)
add_subdirectory(CartesianGrid)
add_subdirectory(IntegerMaps)
add_subdirectory(Trees)

add_sources( 
            ./charMap_class.f90
            ./dictionary_class.f90
            ./dictParser_func.f90
            ./dynArray_class.f90
            ./limb_class.f90
            ./ratint_class.f90
            ./stack_class.f90
           )

# Add Unit Tests to a global List
add_unit_tests(
               ./Tests/charMap_test.f90
               ./Tests/dictionary_test.f90
               ./Tests/dictParser_test.f90
               ./Tests/dynArray_test.f90
               ./Tests/limb1_test.f90
              )

add_integration_tests( 
                      ./Tests/dictParser_iTest.f90
                     )
