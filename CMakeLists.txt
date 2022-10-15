cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(spinner)

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG "release-1.11.0"
)
FetchContent_MakeAvailable(googletest)
include_directories(${CMAKE_SOURCE_DIR}/libs/googletest/googletest/include)

FetchContent_Declare(
        wignersymbols
        GIT_REPOSITORY https://github.com/joeydumont/wignerSymbols.git
        GIT_TAG "9c39f569128fe4cb4835fc852c61f45d0b94aba4"
)
# TODO: do not build tests of wignersymbols
#FetchContent_MakeAvailable(wignersymbols)
#include_directories(${wignersymbols_SOURCE_DIR}/include/)

FetchContent_Declare(
        stlbfgs
        GIT_REPOSITORY https://github.com/ssloy/stlbfgs.git
        GIT_TAG "0dedf926faf66cf95fffbdc3ff5c3e0067b70296"
)
FetchContent_MakeAvailable(stlbfgs)
include_directories(${stlbfgs_SOURCE_DIR}/)

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(${CMAKE_SOURCE_DIR})

add_subdirectory(src)
add_subdirectory(tests)
