add_library(components
        space/Symmetrizer.cpp space/Symmetrizer.h
        space/TzSorter.cpp space/TzSorter.h
        space/NonAbelianSimplifier.cpp space/NonAbelianSimplifier.h)

set_target_properties(components PROPERTIES LINKER_LANGUAGE CXX)