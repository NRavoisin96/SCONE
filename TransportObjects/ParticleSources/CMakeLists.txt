# Add Source Files to the global list
add_subdirectory(ConfigurableSources)
add_subdirectory(DistributedParticleSources)

add_sources( 
             ./source_inter.f90
	           ./sourceFactory_func.f90
            )