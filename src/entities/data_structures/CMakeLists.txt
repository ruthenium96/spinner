# TODO: implement other sparse matrices:
add_library(data_structures
        hashmap/EmhashSparseSemiunitaryMatrix.cpp
        hashmap/as_default/EmhashSparseTransformAsDefault.cpp
        hashmap/HashmapFactories.cpp
        hashmap/EmhashSparseSymmetricMatrix.cpp
        hashmap/EmhashLogic.cpp
        FactoriesList.cpp)

message("Trying to find packages of linear algebra...")

find_package(Armadillo 10.6)

if (ARMADILLO_FOUND)
    message(" * Armadillo was found")
    add_library(arma_structures
            arma/ArmaFactories.cpp arma/ArmaFactories.h
            arma/ArmaDenseVector.cpp arma/ArmaDenseVector.h
            arma/ArmaSparseDiagonalizableMatrix.cpp
            arma/ArmaDenseSemiunitaryMatrix.cpp arma/ArmaDenseSemiunitaryMatrix.h
            arma/ArmaKrylovDenseSemiunitaryMatrix.cpp arma/ArmaKrylovDenseSemiunitaryMatrix.h
            arma/ArmaDenseDiagonalizableMatrix.cpp arma/ArmaDenseDiagonalizableMatrix.h
            arma/ArmaLogic.cpp arma/ArmaLogic.h
            arma/ArmaSparseSemiunitaryMatrix.cpp arma/ArmaSparseSemiunitaryMatrix.h
            )
    # This include is requested sometimes:
    include_directories(${ARMADILLO_INCLUDE_DIRS})
    target_link_libraries(arma_structures ${ARMADILLO_LIBRARIES})
    target_link_libraries(data_structures arma_structures)
else ()
    message(" * Armadillo was NOT found")
endif ()

find_package(Eigen3 3.3)

if (Eigen3_FOUND)
    message(" * Eigen3 was found")
    add_library(eigen_structures
            eigen/EigenFactories.cpp eigen/EigenFactories.h
            eigen/EigenDenseVector.cpp eigen/EigenDenseVector.h
            eigen/EigenDenseDiagonalizableMatrix.cpp eigen/EigenDenseDiagonalizableMatrix.h
            eigen/EigenDenseSemiunitaryMatrix.cpp eigen/EigenDenseSemiunitaryMatrix.h
            eigen/EigenKrylovDenseSemiunitaryMatrix.cpp eigen/EigenKrylovDenseSemiunitaryMatrix.h
            eigen/EigenSparseDiagonalizableMatrix.cpp eigen/EigenSparseDiagonalizableMatrix.h
            eigen/EigenLogic.cpp eigen/EigenLogic.h)
    target_link_libraries(eigen_structures Eigen3::Eigen)
    target_link_libraries(data_structures eigen_structures)
else ()
    message(" * Eigen3 was NOT found")
endif ()

if (ARMADILLO_FOUND)
    message(" *** Armadillo will be used as default linear algebra package")
    target_sources(data_structures PRIVATE
            arma/as_default/ArmaDenseTransformAndDiagonalizeAsDefault.cpp)
elseif (Eigen3_FOUND)
    message(" *** Eigen3 will be used as default linear algebra package")
    target_sources(data_structures PRIVATE
            eigen/as_default/EigenSymmetricMatrixAsDefault.cpp
            )
else ()
    message(FATAL_ERROR "At least one linear algebra package must be installed")
endif ()