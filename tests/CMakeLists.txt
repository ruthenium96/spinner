add_executable(spinner_test
        main.cpp
        time_tests/PerformanceTest.cpp
        unit_tests/group_test.cpp
        unit_tests/consistentModelOptimizationList_tests.cpp
        unit_tests/symmetrizer_tests.cpp
        unit_tests/non_abelian_simplifier_tests.cpp
        matrix_and_spectrum_builders_tests.cpp
        unit_tests/interactions_test.cpp
        non_hamiltonian_operators_tests.cpp
        unit_tests/hamiltonian_operator_tests.cpp
        unit_tests/symbols_tests.cpp
        magnetic_susceptibility_tests.cpp
        integration_tests/spectrum_final_equivalence_tests.cpp)

target_link_libraries(spinner_test gtest gtest_main group space model spectrum common matrix magnetic_susceptibility)