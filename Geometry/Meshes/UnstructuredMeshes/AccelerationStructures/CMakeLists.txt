# Add Source Files to the global list
add_subdirectory(AuxiliaryStructuredCartesianGridAccelerationStructures)
add_sources(
            ./accelerationStructure_inter.f90
            ./accelerationStructureFactory_func.f90
            ./octreeAcceleration_class.f90
            #./analysisDistribution.f90
           )