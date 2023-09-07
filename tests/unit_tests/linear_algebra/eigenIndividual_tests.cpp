#include "abstractIndividual_tests.h"
#include "src/entities/data_structures/eigen/EigenFactories.h"

template<>
std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>
createConcreteFactory<quantum::linear_algebra::EigenDenseTransformAndDiagonalizeFactory>() {
    return std::make_shared<quantum::linear_algebra::EigenDenseTransformAndDiagonalizeFactory>();
}

typedef testing::Types<quantum::linear_algebra::EigenDenseTransformAndDiagonalizeFactory> Eigen;
INSTANTIATE_TYPED_TEST_SUITE_P(
    EigenIndividualTests,
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest,
    Eigen);