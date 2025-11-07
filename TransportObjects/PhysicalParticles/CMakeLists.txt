# Add Source Files to the global list
add_subdirectory(ContinuousEnergyParticles)
add_subdirectory(MultiGroupParticles)

add_sources(
            ./physicalParticle_inter.f90
            ./physicalParticleFactory_func.f90
            ./physicalParticleState_class.f90
            ./testPhysicalParticle_class.f90
           )