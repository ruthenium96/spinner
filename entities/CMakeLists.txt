add_library(space
        BlockProperties.h BlockProperties.cpp
        Space.h Space.cpp Subspace.h Subspace.cpp)

target_link_libraries(space ${ARMADILLO_LIBRARIES})