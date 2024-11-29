add_sources(
            ./box_class.f90
            ./compoundSurface_inter.f90
            ./squareCylinder_class.f90
            ./truncCone_class.f90
            ./truncCylinder_class.f90
           )
           
add_unit_tests(
               ./Tests/box_test.f90
               ./Tests/squareCylinder_test.f90
               ./Tests/truncCone_test.f90
               ./Tests/truncCylinder_test.f90
              )
