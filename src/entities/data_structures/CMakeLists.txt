add_library(data_structures
        UnitarySparseMatrix_std-vector_of_std-maps.cpp
        #        space/SubspaceData_vector_of_sparse_vectors.cpp
        #        space/SubspaceData_sparse_matrix.cpp
        )

find_package(Armadillo)

if (ARMADILLO_FOUND)
    include_directories(${ARMADILLO_INCLUDE_DIRS})
    add_library(arma_structures
            arma/ArmaFactory.cpp arma/ArmaFactory.h
            arma/ArmaMatrix.cpp arma/ArmaMatrix.h
            arma/ArmaVector.cpp arma/ArmaVector.h
            )
    target_link_libraries(arma_structures ${ARMADILLO_LIBRARIES})
    target_link_libraries(data_structures arma_structures)
endif ()

find_package(Eigen3 3.3)

if (Eigen3_FOUND)
    add_library(eigen_structures
            eigen/EigenFactory.cpp eigen/EigenFactory.h
            eigen/EigenMatrix.cpp eigen/EigenMatrix.h
            eigen/EigenVector.cpp eigen/EigenVector.h)
    target_link_libraries(eigen_structures Eigen3::Eigen)
    target_link_libraries(data_structures eigen_structures)
endif ()