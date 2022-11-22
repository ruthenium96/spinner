#include "abstractIndividual_tests.h"
#include "src/entities/data_structures/arma/ArmaDenseFactory.h"

template<>
std::shared_ptr<quantum::linear_algebra::AbstractDenseFactory>
createConcreteFactory<quantum::linear_algebra::ArmaDenseFactory>() {
    return std::make_shared<quantum::linear_algebra::ArmaDenseFactory>();
};

typedef testing::Types<quantum::linear_algebra::ArmaDenseFactory> Arma;
INSTANTIATE_TYPED_TEST_SUITE_P(ArmaIndividualTests, AbstractFactoryIndividualTest, Arma);
