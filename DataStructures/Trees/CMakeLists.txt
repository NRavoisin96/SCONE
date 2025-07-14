# Add Source Files to a global List
add_subdirectory(Nodes)

add_sources(
            ./kdTree_class.f90
            ./tree_inter.f90
            ./octree_class.f90
           )

add_unit_tests(
               ./Tests/kdTree_test.f90
              )