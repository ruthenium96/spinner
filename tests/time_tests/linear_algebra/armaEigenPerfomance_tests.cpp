#include <chrono>
#include <random>

#include "gtest/gtest.h"
#include "tests/tools/AllSymmetricMatrixFactories.h"
#include "tests/tools/GenerateSameMatrix.h"
#include "tests/tools/MeanAndDeviation.h"

// Compare time of eigendecomposition of the same matrix:
TEST(DenseFactoriesPerfomance, eigendecomposition_dense) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 2; size <= 2048; size *= 2) {
        size_t cycles = 2048 * 8 / size;
        std::cout << "SIZE = " << size << std::endl;
        // construct identical dense symmetrical matrix:
        auto matrices = generateDenseDiagonalizableMatrices(
            size,
            constructAllDenseTransformAndDiagonalizeFactories(),
            dist,
            rng);
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

TEST(DenseFactoriesPerfomance, eigendecomposition_sparse) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 2; size <= 2048; size *= 2) {
        size_t cycles = 2048 * 8 / size;
        std::cout << "SIZE = " << size << std::endl;
        // construct identical sparse symmetrical matrix:
        auto matrices = generateSparseDiagonalizableMatrices(
            size,
            constructAllDenseTransformAndDiagonalizeFactories(),
            dist,
            rng);
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
TEST(DenseFactoriesPerfomance, unitary_transformation_dense) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 2; size <= 2048; size *= 2) {
        size_t cycles = 2048 * 16 / size;
        std::cout << "SIZE = " << size << std::endl;
        auto denseDiagonalizableMatrices = generateDenseDiagonalizableMatrices(
            size,
            constructAllDenseTransformAndDiagonalizeFactories(),
            dist,
            rng);
        auto denseUnitaryMatrices = generateDenseUnitaryMatrices(
            size,
            constructAllDenseTransformAndDiagonalizeFactories(),
            dist,
            rng);

        // unitary transformation
        std::function<void(void)> armaUnitaryTransformation = [&]() {
            denseUnitaryMatrices[0]->unitaryTransform(denseDiagonalizableMatrices[0]);
        };
        std::function<void(void)> eigenUnitaryTransformation = [&]() {
            denseUnitaryMatrices[1]->unitaryTransform(denseDiagonalizableMatrices[1]);
        };
        std::cout << "Armadillo:" << std::endl;
        PerformanceTest(armaUnitaryTransformation, cycles);

        std::cout << "Eigen:" << std::endl;
        PerformanceTest(eigenUnitaryTransformation, cycles);
        std::cout << std::endl;
    }
}

// Compare time of unitary transformation (main diagonal only) of the same matrix:
TEST(DenseFactoriesPerfomance, unitary_transformation_and_return_main_diagonal_dense) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 2; size <= 2048; size *= 2) {
        size_t cycles = 2048 * 16 / size;
        std::cout << "SIZE = " << size << std::endl;
        auto denseDiagonalizableMatrices = generateDenseDiagonalizableMatrices(
            size,
            constructAllDenseTransformAndDiagonalizeFactories(),
            dist,
            rng);
        auto denseUnitaryMatrices = generateDenseUnitaryMatrices(
            size,
            constructAllDenseTransformAndDiagonalizeFactories(),
            dist,
            rng);

        // unitary transformation
        std::function<void(void)> armaUnitaryTransformation = [&]() {
            denseUnitaryMatrices[0]->unitaryTransformAndReturnMainDiagonal(
                denseDiagonalizableMatrices[0]);
        };
        std::function<void(void)> eigenUnitaryTransformation = [&]() {
            denseUnitaryMatrices[1]->unitaryTransformAndReturnMainDiagonal(
                denseDiagonalizableMatrices[1]);
        };
        std::cout << "Armadillo:" << std::endl;
        PerformanceTest(armaUnitaryTransformation, cycles);

        std::cout << "Eigen:" << std::endl;
        PerformanceTest(eigenUnitaryTransformation, cycles);
        std::cout << std::endl;
    }
}

// Compare time of unitary transformation (main diagonal only) of the same sparse matrix:
TEST(DenseFactoriesPerfomance, unitary_transformation_and_return_main_diagonal_sparse) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 2; size <= 2048; size *= 2) {
        size_t cycles = 2048 * 16 / size;
        std::cout << "SIZE = " << size << std::endl;
        auto sparseDiagonalizableMatrices = generateSparseDiagonalizableMatrices(
            size,
            constructAllDenseTransformAndDiagonalizeFactories(),
            dist,
            rng);
        auto denseUnitaryMatrix = generateDenseUnitaryMatrices(
            size,
            constructAllDenseTransformAndDiagonalizeFactories(),
            dist,
            rng);

        // unitary transformation
        std::function<void(void)> armaUnitaryTransformation = [&]() {
            denseUnitaryMatrix[0]->unitaryTransformAndReturnMainDiagonal(
                sparseDiagonalizableMatrices[0]);
        };
        std::function<void(void)> eigenUnitaryTransformation = [&]() {
            denseUnitaryMatrix[1]->unitaryTransformAndReturnMainDiagonal(
                sparseDiagonalizableMatrices[1]);
        };
        std::cout << "Armadillo:" << std::endl;
        PerformanceTest(armaUnitaryTransformation, cycles);

        std::cout << "Eigen:" << std::endl;
        PerformanceTest(eigenUnitaryTransformation, cycles);
        std::cout << std::endl;
    }
}
