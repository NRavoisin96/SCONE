# Add Source Files to the global list
add_subdirectory(ParticleSources)
add_subdirectory(PhysicalParticles)
add_subdirectory(VirtualParticles)

add_sources( 
            ./particle_class.f90
            ./particleDungeon_class.f90
            ./testTransportObject_class.f90
            ./transportObject_inter.f90
            ./transportObjectState_class.f90
           )

add_unit_tests( 
               ./Tests/particle_test.f90
               ./Tests/particleDungeon_test.f90
              )
