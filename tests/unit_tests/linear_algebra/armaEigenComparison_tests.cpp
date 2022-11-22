#include <random>

#include "gtest/gtest.h"
#include "tests/tools/GenerateSameDenseMatrix.h"

TEST(linearAlgebraFactories, throw_combination_of_different_objects) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    auto unitaryMatrices = generateUnitaryMatrix(10, constructAllDenseFactories(), dist, rng);
    auto symmetricMatrices = generateSymmetricMatrices(10, constructAllDenseFactories(), dist, rng);

    auto armaUnitaryMatrix = std::move(unitaryMatrices[0]);
    auto eigenUnitaryMatrix = std::move(unitaryMatrices[1]);
    auto armaSymmetricMatrix = std::move(symmetricMatrices[0]);
    auto eigenSymmetricMatrix = std::move(symmetricMatrices[1]);
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
        auto matrices = generateSymmetricMatrices(size, constructAllDenseFactories(), dist, rng);
        auto armaMatrix = std::move(matrices[0]);
        auto eigenMatrix = std::move(matrices[1]);
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
        auto symmetricMatrices =
            generateUnitaryMatrix(size, constructAllDenseFactories(), dist, rng);
        auto unitaryMatrices =
            generateSymmetricMatrices(size, constructAllDenseFactories(), dist, rng);

        auto armaSymmetricMatrixTransformed =
            unitaryMatrices[0]->unitary_transform(symmetricMatrices[0]);
        auto eigenSymmetricMatrixTransformed =
            unitaryMatrices[1]->unitary_transform(symmetricMatrices[1]);

        // check equality:
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < size; ++j) {
                // TODO: epsilon
                EXPECT_NEAR(
                    armaSymmetricMatrixTransformed->at(i, j),
                    eigenSymmetricMatrixTransformed->at(i, j),
                    1e-5);
            }
        }
    }
}