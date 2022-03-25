add_library(components
        space/Symmetrizer.cpp space/Symmetrizer.h
        space/TzSorter.cpp space/TzSorter.h
        space/PositiveProjectionsEliminator.cpp space/PositiveProjectionsEliminator.h
        space/NonAbelianSimplifier.cpp space/NonAbelianSimplifier.h
        )

set_target_properties(components PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(components group space common matrix data_structures spectrum)