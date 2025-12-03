add_subdirectory(ceNeutronData)
add_subdirectory(DataDecks)
add_subdirectory(emissionENDF)
add_subdirectory(mgNeutronData)
add_subdirectory(NuclearDataStructures)
add_subdirectory(NuclideCompositions)
add_subdirectory(PhysicalPropertyLaws)
add_subdirectory(Reactions)
add_subdirectory(testNeutronData)
add_subdirectory(xsPackages)

add_sources(
            ./materialHandle_inter.f90
            ./materialMenu_mod.f90
            ./neutronMaterial_inter.f90
            ./nuclearDatabase_inter.f90
            ./nuclearDataReg_mod.f90
            ./nuclideHandle_inter.f90
            ./nuclideInfo_class.f90
           )

add_unit_tests(
               ./Tests/materialMenu_test.f90
              )


add_integration_tests(
                      ./Tests/aceLibraryRead_iTest.f90
                     )
