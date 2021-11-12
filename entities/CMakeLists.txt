add_library(space
        BlockProperties.h BlockProperties.cpp
        space/Space.h space/Space.cpp
        space/Subspace.h space/Subspace.cpp)

add_library(matrix
        matrix/Matrix.cpp matrix/Matrix.h
        matrix/Submatrix.cpp matrix/Submatrix.h)

add_library(spectrum
        spectrum/Spectrum.cpp spectrum/Spectrum.h spectrum/Subspectrum.cpp spectrum/Subspectrum.h)

#find_package(Eigen3 3.3 REQUIRED)

add_library(data_structures
        data_structures/UnitarySpaseMatrix_std-vector_of_std-maps.cpp
        #        space/SubspaceData_vector_of_sparse_vectors.cpp
        #        space/SubspaceData_sparse_matrix.cpp
        data_structures/UnitarySpaseMatrix.h
        data_structures/DenseMatrix_arma-dmat.cpp
        #        matrix/SubmatrixData_eigen3.cpp
        data_structures/DenseMatrix.h
        )

target_link_libraries(data_structures ${ARMADILLO_LIBRARIES}
#        Eigen3::Eigen
        )