add_library(components
        space/Symmetrizer.cpp space/Symmetrizer.h
        space/TzSorter.cpp space/TzSorter.h
        space/PositiveProjectionsEliminator.cpp space/PositiveProjectionsEliminator.h
        space/NonAbelianSimplifier.cpp space/NonAbelianSimplifier.h
        operator/ScalarProduct.cpp operator/ScalarProduct.h
        operator/ConstantOperator.cpp operator/ConstantOperator.h
        operator/Interaction.h
        )

set_target_properties(components PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(components group space common operator matrix data_structures spectrum)