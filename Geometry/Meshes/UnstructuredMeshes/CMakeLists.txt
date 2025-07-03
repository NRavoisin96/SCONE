# Add Source Files to the global list
add_subdirectory(AccelerationStructures)
add_subdirectory(TopologicalObjects)
add_sources(
            ./OpenFOAMMesh_class.f90
            ./unstructuredMesh_inter.f90
           )

add_integration_tests(
                      ./Tests/OpenFOAMMesh_iTest.f90
                      ./Tests/OpenFOAMMesh_tetrahedral_iTest.f90
                      ./Tests/OpenFOAMMesh_triangulated_iTest.f90
                     )