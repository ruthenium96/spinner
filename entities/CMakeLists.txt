add_library(space
        BlockProperties.h BlockProperties.cpp
        space/Space.h space/Space.cpp space/Subspace.h space/Subspace.cpp
        space/SubspaceData_vector_of_maps.cpp
#        space/SubspaceData_vector_of_sparse_vectors.cpp
#        space/SubspaceData_sparse_matrix.cpp
        space/SubspaceData.h)

add_library(operator
        operator/Operator.cpp operator/Operator.h)

#find_package(Eigen3 3.3 REQUIRED)

add_library(matrix
        matrix/Matrix.cpp matrix/Matrix.h
        matrix/Submatrix.cpp matrix/Submatrix.h
        matrix/SubmatrixData_armadillo.cpp
#        matrix/SubmatrixData_eigen3.cpp
        matrix/SubmatrixData.h)

target_link_libraries(space ${ARMADILLO_LIBRARIES})
target_link_libraries(matrix ${ARMADILLO_LIBRARIES}
#        Eigen3::Eigen
        )