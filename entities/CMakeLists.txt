add_library(space
        BlockProperties.h BlockProperties.cpp
        space/Space.h space/Space.cpp
        space/Subspace.h space/Subspace.cpp)

add_library(operator
        operator/Operator.cpp operator/Operator.h)

add_library(matrix
        matrix/Matrix.cpp matrix/Matrix.h
        matrix/Submatrix.cpp matrix/Submatrix.h
        matrix/SubmatrixData_armadillo.cpp
#        matrix/SubmatrixData_eigen3.cpp
        matrix/SubmatrixData.h)

target_link_libraries(data_structures ${ARMADILLO_LIBRARIES}
#        Eigen3::Eigen
        )