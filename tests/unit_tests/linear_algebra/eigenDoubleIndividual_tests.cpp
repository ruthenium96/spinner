#include "abstractIndividual_tests.h"
#include "src/entities/data_structures/eigen/EigenFactories.h"

using namespace spinner::linear_algebra;

template<>
std::shared_ptr<AbstractDenseTransformAndDiagonalizeFactory>
createConcreteFactory<std::pair<EigenDenseTransformAndDiagonalizeFactory, double>>() {
    return std::make_shared<EigenDenseTransformAndDiagonalizeFactory>();
}

typedef testing::Types<std::pair<EigenDenseTransformAndDiagonalizeFactory, double>> EigenDouble;
INSTANTIATE_TYPED_TEST_SUITE_P(
    EigenDoubleIndividualTests,
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest,
    EigenDouble);