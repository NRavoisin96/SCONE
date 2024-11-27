add_sources(
            ./quadSurface_inter.f90
            ./sphere_class.f90
            ./cylinder_class.f90
            ./aPlane_class.f90
            ./plane_class.f90
           )
            
add_unit_tests(
               ./Tests/sphere_test.f90
               ./Tests/cylinder_test.f90
               ./Tests/aPlane_test.f90
               ./Tests/plane_test.f90
              )
