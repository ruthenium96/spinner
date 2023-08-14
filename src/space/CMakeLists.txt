add_library(space
        ../entities/BlockProperties.h ../entities/BlockProperties.cpp
        Space.h Space.cpp
        Subspace.h Subspace.cpp
        optimization/Symmetrizer.cpp optimization/Symmetrizer.h
        optimization/TzSorter.cpp optimization/TzSorter.h
        optimization/PositiveProjectionsEliminator.cpp optimization/PositiveProjectionsEliminator.h
        optimization/NonAbelianSimplifier.cpp optimization/NonAbelianSimplifier.h
        optimization/OptimizedSpaceConstructor.cpp optimization/OptimizedSpaceConstructor.h
        optimization/S2Transformer.cpp optimization/S2Transformer.h)

FetchContent_Declare(
        wignersymbols
        GIT_REPOSITORY https://github.com/ruthenium96/wignerSymbols.git
        GIT_TAG "deae28c3b0d92579a193517970e59fa78d9cde9c"
)
FetchContent_MakeAvailable(wignersymbols)
include_directories(${wignersymbols_SOURCE_DIR}/include/)

target_link_libraries(space wignerSymbols spin_algebra)

target_link_libraries(space data_structures)