# Add Source Files to a global List
add_sources(
            ./objectKDTree_class.f90
            ./vertexKDTree_class.f90
           )

# Add Unit Tests to a global List
add_unit_tests( 
               ./Tests/vertexKDTree_test.f90
              )