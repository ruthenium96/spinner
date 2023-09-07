#include "abstractIndividual_tests.h"
#include "src/entities/data_structures/arma/ArmaFactories.h"

template<>
std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>
createConcreteFactory<quantum::linear_algebra::ArmaDenseTransformAndDiagonalizeFactory>() {
    return std::make_shared<quantum::linear_algebra::ArmaDenseTransformAndDiagonalizeFactory>();
};

typedef testing::Types<quantum::linear_algebra::ArmaDenseTransformAndDiagonalizeFactory> Arma;
INSTANTIATE_TYPED_TEST_SUITE_P(
    ArmaIndividualTests,
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest,
    Arma);
