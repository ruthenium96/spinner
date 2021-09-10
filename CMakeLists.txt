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
include_directories(${CMAKE_SOURCE_DIR}/libs/googletest/googletest/include)

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

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(${CMAKE_SOURCE_DIR})

add_library(groups
        groups/Group.h groups/Group.cpp
        groups/Group_Info.h groups/Group_Info.cpp
        )

add_library(space
        entities/BlockProperties.h entities/Entity.h
        entities/Space.h entities/Space.cpp
        entities/Subspace.h entities/Subspace.cpp)

add_library(components
        components/symmetrizer.cpp components/symmetrizer.h
        components/tz_sorter.cpp components/tz_sorter.h
        )

add_executable(main_components
        main/main_components.cpp
        common/LexicographicIndexConverter.cpp common/LexicographicIndexConverter.h)

target_link_libraries(main_components components groups space)

set_target_properties(components PROPERTIES LINKER_LANGUAGE CXX)

add_subdirectory(test)
