# Add Source Files to the global list
add_sources(
            ./cellZone_class.f90
            ./cellZoneShelf_class.f90
            ./edge_class.f90
            ./edgeShelf_class.f90
            ./element_class.f90
            ./elementShelf_class.f90
            ./face_class.f90
            ./faceShelf_class.f90
            ./OpenFOAMMesh_class.f90
            ./pyramid_class.f90
            ./pyramidShelf_class.f90
            ./tetrahedron_class.f90
            ./tetrahedronShelf_class.f90
            ./triangle_class.f90
            ./triangleShelf_class.f90
            ./vertex_class.f90
            ./vertexShelf_class.f90
            )

add_integration_tests(
                      ./Tests/meshShelf_iTest.f90
                      ./Tests/OpenFOAMMesh_iTest.f90
                      )
