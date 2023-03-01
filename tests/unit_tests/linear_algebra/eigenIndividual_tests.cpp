#include "abstractIndividual_tests.h"
#include "src/entities/data_structures/eigen/EigenFactories.h"

template<>
std::shared_ptr<quantum::linear_algebra::AbstractSymmetricMatrixFactory>
createConcreteFactory<quantum::linear_algebra::EigenSymmetricMatrixFactory>() {
    return std::make_shared<quantum::linear_algebra::EigenSymmetricMatrixFactory>();
}

typedef testing::Types<quantum::linear_algebra::EigenSymmetricMatrixFactory> Eigen;
INSTANTIATE_TYPED_TEST_SUITE_P(EigenIndividualTests, AbstractFactoryIndividualTest, Eigen);