add_sources(
            ./fissionMG_class.f90
            ./multiScatterMG_class.f90
            ./multiScatterP1MG_class.f90
            ./reactionMG_inter.f90
           )

add_unit_tests(
               ./Tests/fissionMG_test.f90
               ./Tests/multiScatterMG_test.f90
               ./Tests/multiScatterP1MG_test.f90
              )