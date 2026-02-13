#include "abstractIndividual_tests.h"
#include "src/linalg_structures/arma/ArmaFactories.h"

using namespace spinner::linalg_structures;

template<>
std::shared_ptr<AbstractDenseTransformAndDiagonalizeFactory>
createConcreteFactory<std::pair<ArmaDenseTransformAndDiagonalizeFactory, double>>() {
    auto factory =
        std::make_shared<ArmaDenseTransformAndDiagonalizeFactory>();
    factory->setPrecision(DOUBLE);
    return factory;
};

typedef testing::Types<std::pair<ArmaDenseTransformAndDiagonalizeFactory, double>> ArmaDouble;
INSTANTIATE_TYPED_TEST_SUITE_P(
    ArmaDoubleIndividualTests,
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest,
    ArmaDouble);