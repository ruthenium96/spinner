# TODO: implement other sparse matrices:
add_library(data_structures
        std/StdSparseMatrix.cpp
        std/StdSparseAsDefault.cpp
        )

message("Trying to find packages of linear algebra...")

find_package(Armadillo)

if (ARMADILLO_FOUND)
    message(" * Armadillo was found")
    add_library(arma_structures
            arma/ArmaDenseFactory.cpp arma/ArmaDenseFactory.h
            arma/ArmaDenseMatrix.cpp arma/ArmaDenseMatrix.h
            arma/ArmaDenseVector.cpp arma/ArmaDenseVector.h
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
            eigen/EigenDenseFactory.cpp eigen/EigenDenseFactory.h
            eigen/EigenDenseMatrix.cpp eigen/EigenDenseMatrix.h
            eigen/EigenDenseVector.cpp eigen/EigenDenseVector.h)
    target_link_libraries(eigen_structures Eigen3::Eigen)
    target_link_libraries(data_structures eigen_structures)
else ()
    message(" * Eigen3 was NOT found")
endif ()

if (ARMADILLO_FOUND)
    message(" *** Armadillo will be used as default linear algebra package")
    target_sources(data_structures PRIVATE arma/ArmaDenseAsDefault.cpp)
elseif (Eigen3_FOUND)
    message(" *** Eigen3 will be used as default linear algebra package")
    target_sources(data_structures PRIVATE eigen/EigenDenseAsDefault.cpp)
else ()
    message(FATAL_ERROR "At least one linear algebra package must be installed")
endif ()