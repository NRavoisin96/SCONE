add_subdirectory(ThermalScatteringCEReactions)

add_sources(
            ./elasticNeutronScatter_class.f90
            ./fissionCE_class.f90
            ./neutronScatter_class.f90
            ./pureCapture_class.f90
            ./uncorrelatedReactionCE_inter.f90
           )

add_unit_tests(
               ./Tests/pureCapture_test.f90
              )

add_integration_tests(
                      ./Tests/elasticScattering_iTest.f90
                      ./Tests/fissionCE_iTest.f90
                      ./Tests/neutronScattering_iTest.f90
                     )