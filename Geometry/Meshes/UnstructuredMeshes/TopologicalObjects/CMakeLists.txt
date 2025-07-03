# Add Source Files to the global list
add_subdirectory(Factories)
add_subdirectory(Shelves)
add_sources(
            ./edge_class.f90
            ./element_class.f90
            ./face_class.f90
            ./topologicalObject_inter.f90
            ./vertex_class.f90
           )