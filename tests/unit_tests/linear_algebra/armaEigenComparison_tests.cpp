#include <random>

#include "gtest/gtest.h"
#include "src/entities/data_structures/arma/ArmaFactory.h"
#include "src/entities/data_structures/eigen/EigenFactory.h"

void multiplyColumnByMinusOne(
    std::unique_ptr<quantum::linear_algebra::AbstractMatrix>& rhs,
    size_t column) {
    for (size_t i = 0; i < rhs->size_rows(); ++i) {
        rhs->assign_to_position(-rhs->at(i, column), i, column);
    }
}

void makeUnitaryMatrixSame(
    const std::unique_ptr<quantum::linear_algebra::AbstractMatrix>& lhs,
    std::unique_ptr<quantum::linear_algebra::AbstractMatrix>& rhs) {
    for (size_t column = 0; column < lhs->size_cols(); ++column) {
        if (std::abs(lhs->at(0, column) - (-rhs->at(0, column))) < 1e-6) {
            multiplyColumnByMinusOne(rhs, column);
        }
    }
}

// TODO: I guess, it will be better to return std::vector<> or something like that

std::pair<
    std::unique_ptr<quantum::linear_algebra::AbstractMatrix>,
    std::unique_ptr<quantum::linear_algebra::AbstractMatrix>>
generatePairSymmetricMatrices(
    size_t size,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng) {
    auto armaFactory = std::make_shared<quantum::linear_algebra::ArmaFactory>();
    auto eigenFactory = std::make_shared<quantum::linear_algebra::EigenFactory>();

    // construct identical matrices:
    auto armaMatrix = armaFactory->createMatrix();
    auto eigenMatrix = eigenFactory->createMatrix();
    armaMatrix->resize(size, size);
    eigenMatrix->resize(size, size);

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            double value = dist(rng);
            armaMatrix->assign_to_position(value, i, j);
            armaMatrix->assign_to_position(value, j, i);
            eigenMatrix->assign_to_position(value, i, j);
            eigenMatrix->assign_to_position(value, j, i);
        }
    }
    std::pair<
        std::unique_ptr<quantum::linear_algebra::AbstractMatrix>,
        std::unique_ptr<quantum::linear_algebra::AbstractMatrix>>
        answer;
    answer.first = std::move(armaMatrix);
    answer.second = std::move(eigenMatrix);

    return answer;
}

std::pair<
    std::unique_ptr<quantum::linear_algebra::AbstractMatrix>,
    std::unique_ptr<quantum::linear_algebra::AbstractMatrix>>
generatePairUnitaryMatrix(
    size_t size,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng) {
    // construct symmetrical matrix:
    auto [armaSymmetricalMatrix, eigenSymmetricalMatrix] =
        generatePairSymmetricMatrices(size, dist, rng);

    // construct unitary matrix as eigenvectors matrix:
    auto armaUnitaryMatrix = armaSymmetricalMatrix->diagonalizeValuesVectors().eigenvectors;
    auto eigenUnitaryMatrix = eigenSymmetricalMatrix->diagonalizeValuesVectors().eigenvectors;

    makeUnitaryMatrixSame(armaUnitaryMatrix, eigenUnitaryMatrix);

    std::pair<
        std::unique_ptr<quantum::linear_algebra::AbstractMatrix>,
        std::unique_ptr<quantum::linear_algebra::AbstractMatrix>>
        answer;
    answer.first = std::move(armaUnitaryMatrix);
    answer.second = std::move(eigenUnitaryMatrix);

    return answer;
}

TEST(linearAlgebraFactories, throw_combination_of_different_objects) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    auto [armaUnitaryMatrix, eigenUnitaryMatrix] = generatePairUnitaryMatrix(10, dist, rng);
    auto [armaSymmetricMatrix, eigenSymmetricMatrix] = generatePairSymmetricMatrices(10, dist, rng);

    // unitary_transform
    EXPECT_ANY_THROW(eigenUnitaryMatrix->unitary_transform(armaSymmetricMatrix));
    EXPECT_ANY_THROW(armaUnitaryMatrix->unitary_transform(eigenSymmetricMatrix));
    EXPECT_NO_THROW(eigenUnitaryMatrix->unitary_transform(eigenSymmetricMatrix));
    EXPECT_NO_THROW(armaUnitaryMatrix->unitary_transform(armaSymmetricMatrix));

    auto armaEigenVector = armaSymmetricMatrix->diagonalizeValues();
    auto eigenEigenVector = eigenSymmetricMatrix->diagonalizeValues();

    // concatenate_with
    EXPECT_ANY_THROW(eigenEigenVector->concatenate_with(armaEigenVector));
    EXPECT_ANY_THROW(armaEigenVector->concatenate_with(eigenEigenVector));
    EXPECT_NO_THROW(eigenEigenVector->concatenate_with(eigenEigenVector));
    EXPECT_NO_THROW(armaEigenVector->concatenate_with(armaEigenVector));

    // dot
    EXPECT_ANY_THROW(eigenEigenVector->dot(armaEigenVector));
    EXPECT_ANY_THROW(armaEigenVector->dot(eigenEigenVector));
    EXPECT_NO_THROW(eigenEigenVector->dot(eigenEigenVector));
    EXPECT_NO_THROW(armaEigenVector->dot(armaEigenVector));

    // element_wise_multiplication
    EXPECT_ANY_THROW(eigenEigenVector->element_wise_multiplication(armaEigenVector));
    EXPECT_ANY_THROW(armaEigenVector->element_wise_multiplication(eigenEigenVector));
    EXPECT_NO_THROW(eigenEigenVector->element_wise_multiplication(eigenEigenVector));
    EXPECT_NO_THROW(armaEigenVector->element_wise_multiplication(armaEigenVector));
}

// Check if different linear algebra packages make the same eigen decomposition
TEST(linearAlgebraFactories, eigendecomposition) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 2; size < 100; ++size) {
        // construct identical symmetrical matrix:
        auto [armaMatrix, eigenMatrix] = generatePairSymmetricMatrices(size, dist, rng);
        // only-values-eigendecomposition:
        {
            // decomposition
            auto armaOnlyValueVector = armaMatrix->diagonalizeValues();
            auto eigenOnlyValueVector = eigenMatrix->diagonalizeValues();
            // check equality:
            for (size_t i = 0; i < size; ++i) {
                // TODO: epsilon
                EXPECT_NEAR(armaOnlyValueVector->at(i), eigenOnlyValueVector->at(i), 1e-6);
            }
        }
        // eigendecomposition-with-eigenvectors
        {
            // decomposition
            auto armaEigenCouple = armaMatrix->diagonalizeValuesVectors();
            auto eigenEigenCouple = eigenMatrix->diagonalizeValuesVectors();
            makeUnitaryMatrixSame(armaEigenCouple.eigenvectors, eigenEigenCouple.eigenvectors);
            // check equality:
            for (size_t i = 0; i < size; ++i) {
                // TODO: epsilon
                EXPECT_NEAR(
                    armaEigenCouple.eigenvalues->at(i),
                    eigenEigenCouple.eigenvalues->at(i),
                    1e-6);
                for (size_t j = 0; j < size; ++j) {
                    // TODO: epsilon
                    EXPECT_NEAR(
                        armaEigenCouple.eigenvectors->at(i, j),
                        eigenEigenCouple.eigenvectors->at(i, j),
                        1e-6);
                }
            }
        }
    }
}

// Check if different linear algebra packages make the same unitary transformation
TEST(linearAlgebraFactories, unitary_transformation) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 2; size < 100; ++size) {
        auto [armaUnitaryMatrix, eigenUnitaryMatrix] = generatePairUnitaryMatrix(size, dist, rng);
        auto [armaSymmetricMatrix, eigenSymmetricMatrix] =
            generatePairSymmetricMatrices(size, dist, rng);

        auto armaSymmetricMatrixTransformed =
            armaUnitaryMatrix->unitary_transform(armaSymmetricMatrix);
        auto eigenSymmetricMatrixTransformed =
            eigenUnitaryMatrix->unitary_transform(eigenSymmetricMatrix);

        // check equality:
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < size; ++j) {
                // TODO: epsilon
                EXPECT_NEAR(
                    armaSymmetricMatrixTransformed->at(i, j),
                    eigenSymmetricMatrixTransformed->at(i, j),
                    1e-6);
            }
        }
    }
}