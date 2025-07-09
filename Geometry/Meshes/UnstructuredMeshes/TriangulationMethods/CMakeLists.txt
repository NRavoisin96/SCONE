# Add Source Files to the global list
add_sources(
            ./centroidTriangulationMethod_class.f90
            ./DompierreTriangulationMethod_class.f90
            ./noTriangulation_class.f90
            ./triangulationFactory_func.f90
            ./triangulationMethod_inter.f90
           )

add_integration_tests(
                      ./Tests/centroidTriangulationMethod_iTest.f90
                      ./Tests/DompierreTriangulationMethod_iTest.f90
                     )