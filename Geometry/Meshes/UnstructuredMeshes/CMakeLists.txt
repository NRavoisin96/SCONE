# Add Source Files to the global list
add_subdirectory(TriangulatedMeshes)
add_sources(
            ./edge_class.f90
            ./edgeShelf_class.f90
            ./element_class.f90
            ./elementShelf_class.f90
            ./face_class.f90
            ./faceShelf_class.f90
            ./OpenFOAMFunctions.f90
            ./OpenFOAMMesh_class.f90
            ./unstructuredMesh_inter.f90
            ./vertex_class.f90
            ./vertexShelf_class.f90
           )

add_integration_tests(
                      ./Tests/OpenFOAMMesh_iTest.f90
                     )