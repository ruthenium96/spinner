#include "abstractIndividual_tests.h"
#include "src/entities/data_structures/arma/ArmaFactories.h"

template<>
std::shared_ptr<quantum::linear_algebra::AbstractSymmetricMatrixFactory>
createConcreteFactory<quantum::linear_algebra::ArmaSymmetricMatrixFactory>() {
    return std::make_shared<quantum::linear_algebra::ArmaSymmetricMatrixFactory>();
};

typedef testing::Types<quantum::linear_algebra::ArmaSymmetricMatrixFactory> Arma;
INSTANTIATE_TYPED_TEST_SUITE_P(ArmaIndividualTests, AbstractFactoryIndividualTest, Arma);
