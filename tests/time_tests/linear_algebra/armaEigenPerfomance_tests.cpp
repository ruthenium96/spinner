#include <chrono>
#include <random>

#include "gtest/gtest.h"
#include "tests/tools/GenerateSameDenseMatrix.h"
#include "tests/tools/MeanAndDeviation.h"

// Compare time of eigendecomposition of the same matrix:
TEST(DenseFactoriesPerfomance, eigendecomposition) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 2; size <= 2048; size *= 2) {
        size_t cycles = 2048 * 8 / size;
        std::cout << "SIZE = " << size << std::endl;
        // construct identical symmetrical matrix:
        auto matrices = generateSymmetricMatrices(size, constructAllDenseFactories(), dist, rng);
        auto armaMatrix = std::move(matrices[0]);
        auto eigenMatrix = std::move(matrices[1]);
        // only-values-eigendecomposition:
        std::function<void(void)> armaOnlyValuesDecomposition = [&]() {
            armaMatrix->diagonalizeValues();
        };
        std::function<void(void)> eigenOnlyValuesDecomposition = [&]() {
            eigenMatrix->diagonalizeValues();
        };
        // eigendecomposition-with-eigenvectors
        std::function<void(void)> armaBothDecomposition = [&]() {
            armaMatrix->diagonalizeValuesVectors();
        };
        std::function<void(void)> eigenBothDecomposition = [&]() {
            eigenMatrix->diagonalizeValuesVectors();
        };

        std::cout << "Armadillo:" << std::endl;
        std::cout << "Only eigenvalues:                  ";
        PerformanceTest(armaOnlyValuesDecomposition, cycles);
        std::cout << "Both eigenvalues and eigenvectors: ";
        PerformanceTest(armaBothDecomposition, cycles);

        std::cout << "Eigen:" << std::endl;
        std::cout << "Only eigenvalues:                  ";
        PerformanceTest(eigenOnlyValuesDecomposition, cycles);
        std::cout << "Both eigenvalues and eigenvectors: ";
        PerformanceTest(eigenBothDecomposition, cycles);
        std::cout << std::endl;
    }
}

// Compare time of unitary transformation of the same matrix:
TEST(DenseFactoriesPerfomance, unitary_transformation) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 2; size <= 2048; size *= 2) {
        size_t cycles = 2048 * 16 / size;
        std::cout << "SIZE = " << size << std::endl;
        auto symmetricMatrices =
            generateUnitaryMatrix(size, constructAllDenseFactories(), dist, rng);
        auto unitaryMatrices =
            generateSymmetricMatrices(size, constructAllDenseFactories(), dist, rng);

        // unitary transformation
        std::function<void(void)> armaUnitaryTransformation = [&]() {
            unitaryMatrices[0]->unitary_transform(symmetricMatrices[0]);
        };
        std::function<void(void)> eigenUnitaryTransformation = [&]() {
            unitaryMatrices[1]->unitary_transform(symmetricMatrices[1]);
        };
        std::cout << "Armadillo:" << std::endl;
        PerformanceTest(armaUnitaryTransformation, cycles);

        std::cout << "Eigen:" << std::endl;
        PerformanceTest(eigenUnitaryTransformation, cycles);
        std::cout << std::endl;
    }
}