cmake_minimum_required(VERSION 3.14)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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

add_library(old_blm
        old_code/indexes/blm.h
        old_code/indexes/base_blm.h old_code/indexes/base_blm.cpp
        old_code/indexes/base_blm_C1.h old_code/indexes/base_blm_C1.cpp
        old_code/indexes/base_blm_C2.h old_code/indexes/base_blm_C2.cpp old_code/indexes/quantum_numbers.h)

add_library(old_matrices
        old_code/matrices/base_matrices.h old_code/matrices/base_matrices.cpp
        old_code/matrices/full.h old_code/matrices/base_matrices.h old_code/matrices/sum_blocked.h
        old_code/matrices/sum_blocked_only_middle.h old_code/matrices/sum_s2_blocked.h)

include_directories(${CMAKE_CURRENT_LIST_DIR})
add_executable(components
        main/main_components.cpp 
        components/c2_symmetrizer.cpp)

set_target_properties(components PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(old_blm ${ARMADILLO_LIBRARIES} wignerSymbols ${Boost_LIBRARIES})
target_link_libraries(old_matrices ${ARMADILLO_LIBRARIES} old_blm)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable(old_main old_code/main.cpp)

set_target_properties(old_main PROPERTIES
        CMAKE_CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO)

add_dependencies(old_main old_blm old_matrices)
target_link_libraries(old_main old_blm old_matrices)

option(PACKAGE_TESTS "Build the tests" ON)
if(PACKAGE_TESTS)
    enable_testing()
    include(GoogleTest)
    add_subdirectory(old_code/tests)
endif()