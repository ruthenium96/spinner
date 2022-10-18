message("Trying to find system googletest...")
find_package(GTest)
if (GTest_FOUND)
    include_directories(${GTest_INCLUDE_DIRS})
    message(" *** System GTest was found")
else ()
    include(FetchContent)
    FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG "release-1.11.0"
    )
    FetchContent_MakeAvailable(googletest)
    include_directories(${CMAKE_SOURCE_DIR}/libs/googletest/googletest/include)
    message(" ***  GTest was downloaded")
endif ()

add_executable(spinner_test
        main.cpp
        time_tests/PerformanceTest.cpp
        integration_tests/spectrum_final_equivalence_tests.cpp
        integration_tests/fitting_theoretical_curves_tests.cpp
        integration_tests/simple_analytical_dependencies_tests.cpp
        unit_tests/group_test.cpp
        unit_tests/consistentModelOptimizationList_tests.cpp
        unit_tests/symmetrizer_tests.cpp
        unit_tests/non_abelian_simplifier_tests.cpp
        unit_tests/interactions_test.cpp
        unit_tests/hamiltonian_operator_tests.cpp
        unit_tests/symbols_tests.cpp
        matrix_and_spectrum_builders_tests.cpp
        non_hamiltonian_operators_tests.cpp
        unit_tests/magnetic_susceptibility_tests.cpp
        )

find_package(Armadillo)
find_package(Eigen3 3.3)

if (ARMADILLO_FOUND)
    target_sources(spinner_test PRIVATE unit_tests/linear_algebra/armaIndividual_tests.cpp)
endif ()

if (Eigen3_FOUND)
    target_sources(spinner_test PRIVATE unit_tests/linear_algebra/eigenIndividual_tests.cpp)
endif ()

if (ARMADILLO_FOUND AND Eigen3_FOUND)
    target_sources(spinner_test PRIVATE unit_tests/linear_algebra/armaEigenComparison_tests.cpp)
endif ()

target_link_libraries(spinner_test gtest gtest_main group space model spectrum common matrix magnetic_susceptibility)