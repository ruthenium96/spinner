#include "abstractIndividual_tests.h"
#include "src/entities/data_structures/arma/ArmaFactory.h"

template<>
std::shared_ptr<quantum::linear_algebra::AbstractFactory>
createConcreteFactory<quantum::linear_algebra::ArmaFactory>() {
    return std::make_shared<quantum::linear_algebra::ArmaFactory>();
};

typedef testing::Types<quantum::linear_algebra::ArmaFactory> Arma;
INSTANTIATE_TYPED_TEST_SUITE_P(ArmaIndividualTests, AbstractFactoryIndividualTest, Arma);
