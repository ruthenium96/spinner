#include "abstractIndividual_tests.h"
#include "src/entities/data_structures/arma/ArmaFactories.h"

using namespace spinner::linear_algebra;

template<>
std::shared_ptr<AbstractDenseTransformAndDiagonalizeFactory>
createConcreteFactory<std::pair<ArmaDenseTransformAndDiagonalizeFactory, float>>() {
    auto factory =
        std::make_shared<ArmaDenseTransformAndDiagonalizeFactory>();
    factory->setPrecision(SINGLE);
    return factory;
};


typedef testing::Types<std::pair<ArmaDenseTransformAndDiagonalizeFactory, float>> ArmaFloat;
INSTANTIATE_TYPED_TEST_SUITE_P(
    ArmaFloatIndividualTests,
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest,
    ArmaFloat);
