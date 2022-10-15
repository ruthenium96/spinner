#include "abstractIndividual_tests.h"
#include "src/entities/data_structures/eigen/EigenFactory.h"

template<>
std::shared_ptr<quantum::linear_algebra::AbstractFactory>
createConcreteFactory<quantum::linear_algebra::EigenFactory>() {
    return std::make_shared<quantum::linear_algebra::EigenFactory>();
}

typedef testing::Types<quantum::linear_algebra::EigenFactory> Eigen;
INSTANTIATE_TYPED_TEST_SUITE_P(EigenIndividualTests, AbstractFactoryIndividualTest, Eigen);