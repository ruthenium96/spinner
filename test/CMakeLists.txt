add_executable(julytest
        ${CMAKE_SOURCE_DIR}/test/main.cpp
        ${CMAKE_SOURCE_DIR}/test/PerformanceTest.cpp
        ${CMAKE_SOURCE_DIR}/test/group_test.cpp
        ${CMAKE_SOURCE_DIR}/test/symmetrizer_tests.cpp
        ${CMAKE_SOURCE_DIR}/test/non_abelian_simplifier_tests.cpp
        ${CMAKE_SOURCE_DIR}/test/matrix_and_spectrum_builders_tests.cpp
        ${CMAKE_SOURCE_DIR}/test/interactions_test.cpp
        ${CMAKE_SOURCE_DIR}/test/non_hamiltonian_operators_tests.cpp
        ${CMAKE_SOURCE_DIR}/test/hamiltonian_operator_tests.cpp
        ${CMAKE_SOURCE_DIR}/test/symbols_tests.cpp)

target_link_libraries(julytest gtest gtest_main group components space common matrix magnetic_susceptibility)