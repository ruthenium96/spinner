add_library(components
        space/Symmetrizer.cpp space/Symmetrizer.h
        space/TzSorter.cpp space/TzSorter.h
        space/NonAbelianSimplifier.cpp space/NonAbelianSimplifier.h
        operator/ScalarProduct.cpp operator/ScalarProduct.h
        operator/ConstantOperator.cpp operator/ConstantOperator.h
        operator/Interaction.h
        matrix/MatrixBuilder.cpp matrix/MatrixBuilder.h
        spectrum/SpectrumBuilder.cpp spectrum/SpectrumBuilder.h)

set_target_properties(components PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(components groups space common operator matrix data_structures spectrum)