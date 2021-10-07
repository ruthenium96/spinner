add_library(components
        Symmetrizer.cpp Symmetrizer.h
        TzSorter.cpp TzSorter.h
        NonAbelianSimplifier.cpp NonAbelianSimplifier.h)

set_target_properties(components PROPERTIES LINKER_LANGUAGE CXX)