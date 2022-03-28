add_executable(julytest
        ${CMAKE_SOURCE_DIR}/tests/main.cpp
        ${CMAKE_SOURCE_DIR}/tests/PerformanceTest.cpp
        ${CMAKE_SOURCE_DIR}/tests/group_test.cpp
        ${CMAKE_SOURCE_DIR}/tests/symmetrizer_tests.cpp
        ${CMAKE_SOURCE_DIR}/tests/non_abelian_simplifier_tests.cpp
        ${CMAKE_SOURCE_DIR}/tests/matrix_and_spectrum_builders_tests.cpp
        ${CMAKE_SOURCE_DIR}/tests/interactions_test.cpp
        ${CMAKE_SOURCE_DIR}/tests/non_hamiltonian_operators_tests.cpp
        ${CMAKE_SOURCE_DIR}/tests/hamiltonian_operator_tests.cpp
        ${CMAKE_SOURCE_DIR}/tests/symbols_tests.cpp
        ${CMAKE_SOURCE_DIR}/tests/magnetic_susceptibility_tests.cpp
        ${CMAKE_SOURCE_DIR}/tests/spectrum_final_equivalence_tests.cpp)

target_link_libraries(julytest gtest gtest_main group space model spectrum common matrix magnetic_susceptibility)