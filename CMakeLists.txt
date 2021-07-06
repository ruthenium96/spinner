cmake_minimum_required(VERSION 3.13)

add_subdirectory(libs)

project(july)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_library(blm
        indexes/blm.h
        indexes/base_blm.h indexes/base_blm.cpp
        indexes/symm_blm.h indexes/symm_blm.cpp
        indexes/P1_blm.h indexes/P1_blm.cpp
        indexes/P2_blm.h indexes/P2_blm.cpp)

add_library(matrices
        matrices/base_matrices.h matrices/base_matrices.cpp)

target_link_libraries(blm ${ARMADILLO_LIBRARIES} wignerSymbols)
target_link_libraries(matrices ${ARMADILLO_LIBRARIES})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable(main main.cpp)

set_target_properties(main PROPERTIES
        CMAKE_CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO)

add_dependencies(main blm matrices)
target_link_libraries(main blm matrices)

option(PACKAGE_TESTS "Build the tests" ON)
if(PACKAGE_TESTS)
    enable_testing()
    include(GoogleTest)
    add_subdirectory(tests)
endif()