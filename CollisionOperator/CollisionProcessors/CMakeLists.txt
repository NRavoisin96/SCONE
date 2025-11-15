# Add Source Files to the global list
add_subdirectory(NeutronCECollisionProcessors)
add_subdirectory(NeutronMGCollisionProcessors)

add_sources( 
            ./collisionProcessor_inter.f90
            ./collisionProcessorFactory_func.f90
            ./populationComber_class.f90
           )