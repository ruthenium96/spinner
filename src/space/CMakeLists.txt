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
        GIT_REPOSITORY https://github.com/joeydumont/wignerSymbols.git
        GIT_TAG "9c39f569128fe4cb4835fc852c61f45d0b94aba4"
)
# TODO: do not build tests of wignersymbols
FetchContent_MakeAvailable(wignersymbols)
include_directories(${wignersymbols_SOURCE_DIR}/include/)

target_link_libraries(space wignerSymbols)

target_link_libraries(space data_structures)