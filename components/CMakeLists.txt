add_library(components
        space/Symmetrizer.cpp space/Symmetrizer.h
        space/TzSorter.cpp space/TzSorter.h
        space/NonAbelianSimplifier.cpp space/NonAbelianSimplifier.h
        operator/IsotropicExchangeHamiltonian.cpp operator/IsotropicExchangeHamiltonian.h
        operator/Interaction.h operator/ScalarProduct.h matrix/MatrixBuilder.cpp matrix/MatrixBuilder.h)

set_target_properties(components PROPERTIES LINKER_LANGUAGE CXX)