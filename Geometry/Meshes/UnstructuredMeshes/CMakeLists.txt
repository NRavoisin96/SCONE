# Add Source Files to the global list
add_subdirectory(AccelerationStructures)
add_subdirectory(Elements)
add_subdirectory(Faces)
add_sources(
            ./edge_class.f90
            ./edgeShelf_class.f90
            ./OpenFOAMMesh_class.f90
            ./unstructuredMesh_inter.f90
            ./vertex_class.f90
            ./vertexShelf_class.f90
           )

add_integration_tests(
                      ./Tests/OpenFOAMMesh_iTest.f90
                      ./Tests/OpenFOAMMesh_tetrahedral_iTest.f90
                      ./Tests/OpenFOAMMesh_triangulated_iTest.f90
                     )