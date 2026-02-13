FetchContent_Declare(
        wignersymbols
        GIT_REPOSITORY https://github.com/ruthenium96/wignerSymbols.git
        GIT_TAG "044ddbc084453f6cdb0e01c3d0b3936c13322826"
)
FetchContent_MakeAvailable(wignersymbols)
include_directories(${wignersymbols_SOURCE_DIR}/include/)

add_library(spin_algebra
        MultiplicityDirectSum.cpp MultiplicityDirectSum.h
        GroupAdapter.cpp GroupAdapter.h
        ClebshGordanCalculator.cpp ClebshGordanCalculator.h)

target_link_libraries(spin_algebra wignerSymbols index_converter)