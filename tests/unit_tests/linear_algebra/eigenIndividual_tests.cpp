#include "abstractIndividual_tests.h"
#include "src/entities/data_structures/eigen/EigenDenseFactory.h"

template<>
std::shared_ptr<quantum::linear_algebra::AbstractDenseFactory>
createConcreteFactory<quantum::linear_algebra::EigenDenseFactory>() {
    return std::make_shared<quantum::linear_algebra::EigenDenseFactory>();
}

typedef testing::Types<quantum::linear_algebra::EigenDenseFactory> Eigen;
INSTANTIATE_TYPED_TEST_SUITE_P(EigenIndividualTests, AbstractFactoryIndividualTest, Eigen);