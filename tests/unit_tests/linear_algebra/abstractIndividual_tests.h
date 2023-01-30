#ifndef SPINNER_ABSTRACT_INDIVIDUAL_TESTS_H
#define SPINNER_ABSTRACT_INDIVIDUAL_TESTS_H

#include "gtest/gtest.h"
#include "src/entities/data_structures/AbstractFactories.h"

template<class T>
std::shared_ptr<quantum::linear_algebra::AbstractSymmetricMatrixFactory> createConcreteFactory();

template<class T>
class AbstractFactoryIndividualTest: public testing::Test {
  protected:
    AbstractFactoryIndividualTest() : factory_(createConcreteFactory<T>()) {}
    std::shared_ptr<quantum::linear_algebra::AbstractSymmetricMatrixFactory> const factory_;
};

TYPED_TEST_SUITE_P(AbstractFactoryIndividualTest);

// TODO: implement tests for abstract factory, abstract matrix and abstract vector
TYPED_TEST_P(AbstractFactoryIndividualTest, NonNullptrObjects) {
    EXPECT_FALSE(this->factory_->createDenseSymmetricMatrix(0) == nullptr);
}

REGISTER_TYPED_TEST_SUITE_P(AbstractFactoryIndividualTest, NonNullptrObjects);

#endif  //SPINNER_ABSTRACT_INDIVIDUAL_TESTS_H
