#include <random>

#include "gtest/gtest.h"
#include "src/entities/data_structures/arma/ArmaFactory.h"
#include "src/entities/data_structures/eigen/EigenFactory.h"

// TODO: tests for inconsistency of different packages

TEST(linearAlgebraFactories, eigendecomposition) {
    auto armaFactory = std::make_shared<quantum::linear_algebra::ArmaFactory>();
    auto eigenFactory = std::make_shared<quantum::linear_algebra::EigenFactory>();

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 2; size < 100; ++size) {
        // construct identical matrix:
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
                        std::abs(armaEigenCouple.eigenvectors->at(i, j)),
                        std::abs(eigenEigenCouple.eigenvectors->at(i, j)),
                        1e-6);
                }
            }
        }
    }
}