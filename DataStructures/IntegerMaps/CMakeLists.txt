# Add Source Files to a global List
add_sources(
            ./longIntMap_class.f90
            ./intMap_class.f90
           )

add_unit_tests(
               ./Tests/intMap_test.f90
               ./Tests/longIntMap_test.f90
              )