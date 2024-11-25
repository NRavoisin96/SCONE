add_sources(
            ./cellUniverse_class.f90
            ./latUniverse_class.f90
            ./meshUniverse_class.f90
            ./pinUniverse_class.f90
            ./rootUniverse_class.f90
            ./uniFills_class.f90
            ./universe_inter.f90
            ./universeFactory_func.f90
            ./universeShelf_class.f90
           )

add_unit_tests(
               ./Tests/universe_test.f90
               ./Tests/cellUniverse_test.f90
               ./Tests/pinUniverse_test.f90
               ./Tests/latUniverse_test.f90
               ./Tests/rootUniverse_test.f90
               ./Tests/uniFills_test.f90
              )

add_integration_tests(
                      ./Tests/meshUniverse_iTest.f90
                      ./Tests/universeShelf_iTest.f90
                     )