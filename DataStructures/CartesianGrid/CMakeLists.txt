# Add Source Files to a global List
add_subdirectory(CartesianGridSingle)
add_subdirectory(CartesianGridMulti)

add_sources( 
            ./cartesianGenericProcedures.f90
            ./cartesianInitProcedures.f90
            ./ragged3dMatrix_class.f90
           )