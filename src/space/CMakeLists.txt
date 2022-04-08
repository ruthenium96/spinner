add_library(space
        ../entities/BlockProperties.h ../entities/BlockProperties.cpp
        Space.h Space.cpp
        Subspace.h Subspace.cpp
        optimization/Symmetrizer.cpp optimization/Symmetrizer.h
        optimization/TzSorter.cpp optimization/TzSorter.h
        optimization/PositiveProjectionsEliminator.cpp optimization/PositiveProjectionsEliminator.h
        optimization/NonAbelianSimplifier.cpp optimization/NonAbelianSimplifier.h
        optimization/OptimizedSpaceConstructor.cpp optimization/OptimizedSpaceConstructor.h)

target_link_libraries(space data_structures)