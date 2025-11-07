FetchContent_Declare(
        wignersymbols
        GIT_REPOSITORY https://github.com/ruthenium96/wignerSymbols.git
        GIT_TAG "ee884a90825e788160fd9837f1fba2f4f844b738"
)
FetchContent_MakeAvailable(wignersymbols)
include_directories(${wignersymbols_SOURCE_DIR}/include/)

add_library(spin_algebra
        MultiplicityDirectSum.cpp MultiplicityDirectSum.h
        GroupAdapter.cpp GroupAdapter.h
        RepresentationsMultiplier.cpp RepresentationsMultiplier.h
        ClebshGordanCalculator.cpp ClebshGordanCalculator.h)

target_link_libraries(spin_algebra wignerSymbols index_converter)