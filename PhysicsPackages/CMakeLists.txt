# Add Source Files to the global list
add_subdirectory(ParticlePhysicsPackages)
add_subdirectory(RandomWalkerPhysicsPackages)
add_subdirectory(RayPhysicsPackages)
add_sources(
            ./physicsPackage_inter.f90
            ./physicsPackageFactory_func.f90
           )
