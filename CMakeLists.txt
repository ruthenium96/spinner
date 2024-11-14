cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(spinner)

include(FetchContent)

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Fast hashmap:
FetchContent_Declare(
        emhash
        GIT_REPOSITORY https://github.com/ktprime/emhash.git
        GIT_TAG "2e3edc378373b331afa8bd570d657bf3f1229b86"
)
FetchContent_MakeAvailable(emhash)
include_directories(${emhash_SOURCE_DIR}/)

include_directories(${CMAKE_SOURCE_DIR})

add_subdirectory(src)
add_subdirectory(tests)
