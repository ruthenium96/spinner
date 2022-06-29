add_library(data_structures
        UnitarySparseMatrix_std-vector_of_std-maps.cpp
        #        space/SubspaceData_vector_of_sparse_vectors.cpp
        #        space/SubspaceData_sparse_matrix.cpp
        )

message("Trying to find packages of linear algebra...")

find_package(Armadillo)

if (ARMADILLO_FOUND)
    message(" * Armadillo was found")
    add_library(arma_structures
            arma/ArmaFactory.cpp arma/ArmaFactory.h
            arma/ArmaMatrix.cpp arma/ArmaMatrix.h
            arma/ArmaVector.cpp arma/ArmaVector.h
            )
    target_link_libraries(arma_structures ${ARMADILLO_LIBRARIES})
    target_link_libraries(data_structures arma_structures)
else ()
    message(" * Armadillo was NOT found")
endif ()

find_package(Eigen3 3.3)

if (Eigen3_FOUND)
    message(" * Eigen3 was found")
    add_library(eigen_structures
            eigen/EigenFactory.cpp eigen/EigenFactory.h
            eigen/EigenMatrix.cpp eigen/EigenMatrix.h
            eigen/EigenVector.cpp eigen/EigenVector.h)
    target_link_libraries(eigen_structures Eigen3::Eigen)
    target_link_libraries(data_structures eigen_structures)
else ()
    message(" * Eigen3 was NOT found")
endif ()

if (ARMADILLO_FOUND)
    message(" *** Armadillo will be used as default linear algebra package")
    add_library(arma_as_default arma/ArmaAsDefault.cpp)
    target_link_libraries(arma_as_default arma_structures)
    target_link_libraries(data_structures arma_as_default)
elseif (Eigen3_FOUND)
    message(" *** Eigen3 will be used as default linear algebra package")
    add_library(eigen_as_default eigen/EigenAsDefault.cpp)
    target_link_libraries(eigen_as_default eigen_structures)
    target_link_libraries(data_structures eigen_as_default)
else ()
    message(FATAL_ERROR "At least one linear algebra package must be installed")
endif ()