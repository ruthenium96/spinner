#include "abstractIndividual_tests.h"
#include "src/linalg_structures/eigen/EigenFactories.h"

using namespace spinner::linalg_structures;

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