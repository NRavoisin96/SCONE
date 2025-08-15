add_subdirectory(PiecewiseConstantScalarFields)

add_sources( 
            ./mixedScalarField_class.f90
            ./scalarField_inter.f90
            ./uniformScalarField_class.f90
           )

add_unit_tests( ./Tests/uniformScalarField_test.f90)
