add_library(space
        BlockProperties.h BlockProperties.cpp
        space/Space.h space/Space.cpp space/Subspace.h space/Subspace.cpp
        space/NewBasisDecomposition_vector_of_maps.cpp
#        space/NewBasisDecomposition_vector_of_sparse_vectors.cpp
#        space/NewBasisDecomposition_sparse_matrix.cpp
        space/NewBasisDecomposition.h )

add_library(operator
        operator/Operator.cpp operator/Operator.h)

add_library(matrix
        matrix/Matrix.cpp matrix/Matrix.h matrix/Submatrix.cpp matrix/Submatrix.h matrix/RawSubmatrixData.cpp matrix/RawSubmatrixData.h)

target_link_libraries(space ${ARMADILLO_LIBRARIES})
target_link_libraries(matrix ${ARMADILLO_LIBRARIES})