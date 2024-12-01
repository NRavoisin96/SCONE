# Add Source Files to the global list
add_subdirectory(UnstructuredMeshes)

add_sources(./meshFactory_func.f90
            ./mesh_inter.f90
            ./meshShelf_class.f90)

add_integration_tests(
                      ./Tests/meshShelf_iTest.f90
                     )