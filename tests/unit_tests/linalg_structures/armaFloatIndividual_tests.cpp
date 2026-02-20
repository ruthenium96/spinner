#include "abstractIndividual_tests.h"
#include "src/linalg_structures/arma/ArmaFactories.h"

using namespace spinner::linalg_structures;

template<>
std::shared_ptr<AbstractDenseTransformAndDiagonalizeFactory>
createConcreteFactory<std::pair<armadillo::ArmaDenseTransformAndDiagonalizeFactory, float>>() {
    auto factory =
        std::make_shared<armadillo::ArmaDenseTransformAndDiagonalizeFactory>();
    factory->setPrecision(SINGLE);
    return factory;
};


typedef testing::Types<std::pair<armadillo::ArmaDenseTransformAndDiagonalizeFactory, float>> ArmaFloat;
INSTANTIATE_TYPED_TEST_SUITE_P(
    ArmaFloatIndividualTests,
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest,
    ArmaFloat);
