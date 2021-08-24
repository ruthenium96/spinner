cmake_minimum_required(VERSION 3.14)

project(july)

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        "release-1.11.0"
)
FetchContent_MakeAvailable(googletest)

FetchContent_Declare(
        wignersymbols
        GIT_REPOSITORY https://github.com/joeydumont/wignerSymbols.git
        GIT_TAG        "9c39f569128fe4cb4835fc852c61f45d0b94aba4"
)
FetchContent_MakeAvailable(wignersymbols)
include_directories(${wignersymbols_SOURCE_DIR}/include/)


find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_library(blm
        indexes/blm.h
        indexes/base_blm.h indexes/base_blm.cpp
        indexes/base_blm_C1.h indexes/base_blm_C1.cpp
        indexes/base_blm_C2.h indexes/base_blm_C2.cpp indexes/quantum_numbers.h)

add_library(matrices
        matrices/base_matrices.h matrices/base_matrices.cpp
        matrices/full.h matrices/base_matrices.h matrices/sum_blocked.h
        matrices/sum_blocked_only_middle.h matrices/sum_s2_blocked.h)

target_link_libraries(blm ${ARMADILLO_LIBRARIES} wignerSymbols ${Boost_LIBRARIES})
target_link_libraries(matrices ${ARMADILLO_LIBRARIES} blm)

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