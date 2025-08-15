# Add Source Files to the global list
add_subdirectory(NeutronCECollisionProcessors)
add_sources( 
            ./collisionProcessor_inter.f90
            ./collisionProcessorFactory_func.f90
            ./neutronMGstd_class.f90
            ./neutronMGimp_class.f90
           )