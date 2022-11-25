#ifndef SPINNER_ABSTRACT_INDIVIDUAL_TESTS_H
#define SPINNER_ABSTRACT_INDIVIDUAL_TESTS_H

#include "gtest/gtest.h"
#include "src/entities/data_structures/AbstractDenseFactory.h"

template<class T>
std::shared_ptr<quantum::linear_algebra::AbstractDenseFactory> createConcreteFactory();

template<class T>
class AbstractFactoryIndividualTest: public testing::Test {
  protected:
    AbstractFactoryIndividualTest() : factory_(createConcreteFactory<T>()) {}
    std::shared_ptr<quantum::linear_algebra::AbstractDenseFactory> const factory_;
};

TYPED_TEST_SUITE_P(AbstractFactoryIndividualTest);

// TODO: implement tests for abstract factory, abstract matrix and abstract vector
TYPED_TEST_P(AbstractFactoryIndividualTest, NonNullptrObjects) {
    EXPECT_FALSE(this->factory_->createMatrix(0, 0) == nullptr);
    EXPECT_FALSE(this->factory_->createVector() == nullptr);
}

REGISTER_TYPED_TEST_SUITE_P(AbstractFactoryIndividualTest, NonNullptrObjects);

#endif  //SPINNER_ABSTRACT_INDIVIDUAL_TESTS_H
