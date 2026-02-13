#include "abstractIndividual_tests.h"
#include "src/linalg_structures/eigen/EigenFactories.h"

using namespace spinner::linalg_structures;

template<>
std::shared_ptr<AbstractDenseTransformAndDiagonalizeFactory>
createConcreteFactory<std::pair<EigenDenseTransformAndDiagonalizeFactory, float>>() {
    auto factory =
        std::make_shared<EigenDenseTransformAndDiagonalizeFactory>();
    factory->setPrecision(SINGLE);
    return factory;
}

typedef testing::Types<std::pair<EigenDenseTransformAndDiagonalizeFactory, float>> EigenFloat;
INSTANTIATE_TYPED_TEST_SUITE_P(
    EigenFloatIndividualTests,
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest,
    EigenFloat);