# Add Source Files to the global list
add_sources(
            ./pyramid_class.f90
            ./pyramidShelf_class.f90
            ./tetrahedron_class.f90
            ./tetrahedronShelf_class.f90
            ./triangle_class.f90
            ./triangleShelf_class.f90
            ./triOpenFOAMMesh_class.f90
            ./triUnstructuredMesh_inter.f90
           )

add_integration_tests(
                      ./Tests/triOpenFOAMMesh_iTest.f90
                     )