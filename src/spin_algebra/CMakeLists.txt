FetchContent_Declare(
        wignersymbols
        GIT_REPOSITORY https://github.com/ruthenium96/wignerSymbols.git
        GIT_TAG "deae28c3b0d92579a193517970e59fa78d9cde9c"
)
FetchContent_MakeAvailable(wignersymbols)
include_directories(${wignersymbols_SOURCE_DIR}/include/)

add_library(spin_algebra
        MultiplicityDirectSum.cpp MultiplicityDirectSum.h
        SSquaredState.cpp SSquaredState.h
        OrderOfSummation.cpp OrderOfSummation.h
        GroupAdapter.cpp GroupAdapter.h
        RepresentationsMultiplier.cpp RepresentationsMultiplier.h
        ClebshGordanCalculator.cpp ClebshGordanCalculator.h
        SSquaredConverter.cpp SSquaredConverter.h)

target_link_libraries(spin_algebra wignerSymbols)