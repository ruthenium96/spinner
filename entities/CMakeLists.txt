add_library(space
        BlockProperties.h BlockProperties.cpp
        space/Space.h space/Space.cpp space/Subspace.h space/Subspace.cpp space/NewBasisDecomposition.cpp space/NewBasisDecomposition.h)

target_link_libraries(space ${ARMADILLO_LIBRARIES})