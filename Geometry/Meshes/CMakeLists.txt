add_subdirectory(OpenFOAMMesh)
# Add Source Files to the global list 
add_sources(./meshFactory_func.f90
            ./mesh_inter.f90
            ./meshShelf_class.f90)
