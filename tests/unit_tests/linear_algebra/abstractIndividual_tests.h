#ifndef SPINNER_ABSTRACT_INDIVIDUAL_TESTS_H
#define SPINNER_ABSTRACT_INDIVIDUAL_TESTS_H

#include "gtest/gtest.h"
#include "src/entities/data_structures/AbstractFactories.h"
#include "tests/tools/GenerateSameMatrix.h"

template<class T>
std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>
createConcreteFactory();

template<class T>
class AbstractDenseTransformAndDiagonalizeFactoryIndividualTest: public testing::Test {
  protected:
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest() :
        factory_(createConcreteFactory<T>()) {}
    std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory> const
        factory_;
};

TYPED_TEST_SUITE_P(AbstractDenseTransformAndDiagonalizeFactoryIndividualTest);

// TODO: implement tests for abstract factory, abstract matrix and abstract vector
TYPED_TEST_P(AbstractDenseTransformAndDiagonalizeFactoryIndividualTest, NonNullptrObjects) {
    EXPECT_FALSE(this->factory_->createDenseDiagonalizableMatrix(0) == nullptr);
}

TYPED_TEST_P(
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest,
    unitary_transformation_and_unitary_transformation_and_return_main_diagonal_equivalence) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 2; size < 100; ++size) {
        auto denseDiagonalizableMatrix =
            generateDenseDiagonalizableMatrix(size, this->factory_, dist, rng);
        auto denseUnitaryMatrix = generateDenseUnitaryMatrix(size, this->factory_, dist, rng);

        auto denseDiagonalizableMatrixTransformedMatrix =
            denseUnitaryMatrix->unitaryTransform(denseDiagonalizableMatrix);
        auto denseDiagonalizableMatrixTransformedMainDiagonal =
            denseUnitaryMatrix->unitaryTransformAndReturnMainDiagonal(denseDiagonalizableMatrix);

        // check equality:
        for (size_t i = 0; i < size; ++i) {
            double epsilon =
                std::abs(denseDiagonalizableMatrixTransformedMainDiagonal->at(i) * 1e-6);
            EXPECT_NEAR(
                denseDiagonalizableMatrixTransformedMainDiagonal->at(i),
                denseDiagonalizableMatrixTransformedMatrix->at(i, i),
                epsilon);
        }
    }
}

TYPED_TEST_P(
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest,
    krylovDiagonalizeValuesAndDiagonalizeValues) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-100, +100);
    
    for (size_t size = 8; size <= 16; size*=2) {
        auto matrix = generateSparseDiagonalizableMatrix(size, this->factory_, dist, rng);
        auto seed_vector = generateOrthDenseVector(size, this->factory_);

        auto krylov_eigenvalues = std::move(matrix->krylovDiagonalizeValues(seed_vector, size).eigenvalues);
        auto exact_eigenvalues = matrix->diagonalizeValues();
        ASSERT_EQ(krylov_eigenvalues->size(), exact_eigenvalues->size());
        for (size_t i = 0; i < krylov_eigenvalues->size(); ++i) {
            double range = std::abs(krylov_eigenvalues->at(i) * 1e-1);
            EXPECT_NEAR(krylov_eigenvalues->at(i), exact_eigenvalues->at(i), range);
        }
    }
}

TYPED_TEST_P(
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest,
    randomUnitVectorsAreUnit) {
    
    for (size_t size = 1; size <= 2048; size*=2) {
        auto vector = this->factory_->createRandomUnitVector(size);

        ASSERT_EQ(vector->size(), size);
        double norm = 0.0;
        for (size_t i = 0; i < size; ++i) {
            double value = vector->at(i);
            norm += value * value; 
        }
        EXPECT_NEAR(norm, 1.0, 1e-6);
    }
}

TYPED_TEST_P(
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest,
    randomUnitVectorsAreUnbiased) {
    for (size_t size = 1024; size <= 65536; size*=2) {
        auto vector = this->factory_->createRandomUnitVector(size);

        ASSERT_EQ(vector->size(), size);
        double sum = 0.0;
        for (size_t i = 0; i < size; ++i) {
            double value = vector->at(i);
            sum += value; 
        }
        sum /= size;
        EXPECT_NEAR(sum, 0.0, 5e-3);
    }
}

REGISTER_TYPED_TEST_SUITE_P(
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest,
    NonNullptrObjects,
    unitary_transformation_and_unitary_transformation_and_return_main_diagonal_equivalence,
    krylovDiagonalizeValuesAndDiagonalizeValues,
    randomUnitVectorsAreUnit,
    randomUnitVectorsAreUnbiased);

#endif  //SPINNER_ABSTRACT_INDIVIDUAL_TESTS_H
