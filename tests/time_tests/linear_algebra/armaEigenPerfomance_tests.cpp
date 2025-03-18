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
        auto armaDoubleMatrix = std::move(matrices[0]);
        auto armaSingleMatrix = std::move(matrices[1]);
        auto eigenDoubleMatrix = std::move(matrices[2]);
        auto eigenSingleMatrix = std::move(matrices[3]);
        // only-values-eigendecomposition:
        std::function<void(void)> armaDoubleOnlyValuesDecomposition = [&]() {
            armaDoubleMatrix->diagonalizeValues();
        };
        std::function<void(void)> armaSingleOnlyValuesDecomposition = [&]() {
            armaSingleMatrix->diagonalizeValues();
        };
        std::function<void(void)> eigenDoubleOnlyValuesDecomposition = [&]() {
            eigenDoubleMatrix->diagonalizeValues();
        };
        std::function<void(void)> eigenSingleOnlyValuesDecomposition = [&]() {
            eigenSingleMatrix->diagonalizeValues();
        };
        // eigendecomposition-with-eigenvectors
        std::function<void(void)> armaDoubleBothDecomposition = [&]() {
            armaDoubleMatrix->diagonalizeValuesVectors();
        };
        std::function<void(void)> armaSingleBothDecomposition = [&]() {
            armaSingleMatrix->diagonalizeValuesVectors();
        };
        std::function<void(void)> eigenDoubleBothDecomposition = [&]() {
            eigenDoubleMatrix->diagonalizeValuesVectors();
        };
        std::function<void(void)> eigenSingleBothDecomposition = [&]() {
            eigenSingleMatrix->diagonalizeValuesVectors();
        };

        std::cout << "Armadillo double precision:" << std::endl;
        std::cout << "Only eigenvalues:                  ";
        PerformanceTest(armaDoubleOnlyValuesDecomposition, cycles);
        std::cout << "Both eigenvalues and eigenvectors: ";
        PerformanceTest(armaDoubleBothDecomposition, cycles);

        std::cout << "Armadillo single precision:" << std::endl;
        std::cout << "Only eigenvalues:                  ";
        PerformanceTest(armaSingleOnlyValuesDecomposition, cycles);
        std::cout << "Both eigenvalues and eigenvectors: ";
        PerformanceTest(armaSingleBothDecomposition, cycles);

        std::cout << "Eigen double precision:" << std::endl;
        std::cout << "Only eigenvalues:                  ";
        PerformanceTest(eigenDoubleOnlyValuesDecomposition, cycles);
        std::cout << "Both eigenvalues and eigenvectors: ";
        PerformanceTest(eigenDoubleBothDecomposition, cycles);

        std::cout << "Eigen single precision:" << std::endl;
        std::cout << "Only eigenvalues:                  ";
        PerformanceTest(eigenSingleOnlyValuesDecomposition, cycles);
        std::cout << "Both eigenvalues and eigenvectors: ";
        PerformanceTest(eigenSingleBothDecomposition, cycles);
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
        auto armaDoubleMatrix = std::move(matrices[0]);
        auto armaSingleMatrix = std::move(matrices[1]);
        auto eigenDoubleMatrix = std::move(matrices[2]);
        auto eigenSingleMatrix = std::move(matrices[3]);
        // only-values-eigendecomposition:
        std::function<void(void)> armaDoubleOnlyValuesDecomposition = [&]() {
            armaDoubleMatrix->diagonalizeValues();
        };
        std::function<void(void)> armaSingleOnlyValuesDecomposition = [&]() {
            armaSingleMatrix->diagonalizeValues();
        };
        std::function<void(void)> eigenDoubleOnlyValuesDecomposition = [&]() {
            eigenDoubleMatrix->diagonalizeValues();
        };
        std::function<void(void)> eigenSingleOnlyValuesDecomposition = [&]() {
            eigenSingleMatrix->diagonalizeValues();
        };
        // eigendecomposition-with-eigenvectors
        std::function<void(void)> armaDoubleBothDecomposition = [&]() {
            armaDoubleMatrix->diagonalizeValuesVectors();
        };
        std::function<void(void)> armaSingleBothDecomposition = [&]() {
            armaSingleMatrix->diagonalizeValuesVectors();
        };
        std::function<void(void)> eigenDoubleBothDecomposition = [&]() {
            eigenDoubleMatrix->diagonalizeValuesVectors();
        };
        std::function<void(void)> eigenSingleBothDecomposition = [&]() {
            eigenSingleMatrix->diagonalizeValuesVectors();
        };

        std::cout << "Armadillo double precision:" << std::endl;
        std::cout << "Only eigenvalues:                  ";
        PerformanceTest(armaDoubleOnlyValuesDecomposition, cycles);
        std::cout << "Both eigenvalues and eigenvectors: ";
        PerformanceTest(armaDoubleBothDecomposition, cycles);

        std::cout << "Armadillo single precision:" << std::endl;
        std::cout << "Only eigenvalues:                  ";
        PerformanceTest(armaSingleOnlyValuesDecomposition, cycles);
        std::cout << "Both eigenvalues and eigenvectors: ";
        PerformanceTest(armaSingleBothDecomposition, cycles);

        std::cout << "Eigen double precision:" << std::endl;
        std::cout << "Only eigenvalues:                  ";
        PerformanceTest(eigenDoubleOnlyValuesDecomposition, cycles);
        std::cout << "Both eigenvalues and eigenvectors: ";
        PerformanceTest(eigenDoubleBothDecomposition, cycles);

        std::cout << "Eigen single precision:" << std::endl;
        std::cout << "Only eigenvalues:                  ";
        PerformanceTest(eigenSingleOnlyValuesDecomposition, cycles);
        std::cout << "Both eigenvalues and eigenvectors: ";
        PerformanceTest(eigenSingleBothDecomposition, cycles);
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
        std::function<void(void)> armaDoubleUnitaryTransformation = [&]() {
            denseUnitaryMatrices[0]->unitaryTransform(denseDiagonalizableMatrices[0]);
        };
        std::function<void(void)> armaSingleUnitaryTransformation = [&]() {
            denseUnitaryMatrices[1]->unitaryTransform(denseDiagonalizableMatrices[1]);
        };
        std::function<void(void)> eigenDoubleUnitaryTransformation = [&]() {
            denseUnitaryMatrices[2]->unitaryTransform(denseDiagonalizableMatrices[2]);
        };
        std::function<void(void)> eigenSingleUnitaryTransformation = [&]() {
            denseUnitaryMatrices[3]->unitaryTransform(denseDiagonalizableMatrices[3]);
        };
        std::cout << "Armadillo double precision:" << std::endl;
        PerformanceTest(armaDoubleUnitaryTransformation, cycles);

        std::cout << "Armadillo single precision:" << std::endl;
        PerformanceTest(armaSingleUnitaryTransformation, cycles);

        std::cout << "Eigen double precision:" << std::endl;
        PerformanceTest(eigenDoubleUnitaryTransformation, cycles);

        std::cout << "Eigen single precision:" << std::endl;
        PerformanceTest(eigenSingleUnitaryTransformation, cycles);
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
        std::function<void(void)> armaDoubleUnitaryTransformation = [&]() {
            denseUnitaryMatrices[0]->unitaryTransformAndReturnMainDiagonal(
                denseDiagonalizableMatrices[0]);
        };
        std::function<void(void)> armaSingleUnitaryTransformation = [&]() {
            denseUnitaryMatrices[1]->unitaryTransformAndReturnMainDiagonal(
                denseDiagonalizableMatrices[1]);
        };
        std::function<void(void)> eigenDoubleUnitaryTransformation = [&]() {
            denseUnitaryMatrices[2]->unitaryTransformAndReturnMainDiagonal(
                denseDiagonalizableMatrices[2]);
        };
        std::function<void(void)> eigenSingleUnitaryTransformation = [&]() {
            denseUnitaryMatrices[3]->unitaryTransformAndReturnMainDiagonal(
                denseDiagonalizableMatrices[3]);
        };
        std::cout << "Armadillo double precision:" << std::endl;
        PerformanceTest(armaDoubleUnitaryTransformation, cycles);

        std::cout << "Armadillo single precision:" << std::endl;
        PerformanceTest(armaSingleUnitaryTransformation, cycles);

        std::cout << "Eigen double precision:" << std::endl;
        PerformanceTest(eigenDoubleUnitaryTransformation, cycles);

        std::cout << "Eigen single precision:" << std::endl;
        PerformanceTest(eigenSingleUnitaryTransformation, cycles);
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
        std::function<void(void)> armaDoubleUnitaryTransformation = [&]() {
            denseUnitaryMatrix[0]->unitaryTransformAndReturnMainDiagonal(
                sparseDiagonalizableMatrices[0]);
        };
        std::function<void(void)> armaSingleUnitaryTransformation = [&]() {
            denseUnitaryMatrix[1]->unitaryTransformAndReturnMainDiagonal(
                sparseDiagonalizableMatrices[1]);
        };
        std::function<void(void)> eigenDoubleUnitaryTransformation = [&]() {
            denseUnitaryMatrix[2]->unitaryTransformAndReturnMainDiagonal(
                sparseDiagonalizableMatrices[2]);
        };
        std::function<void(void)> eigenSingleUnitaryTransformation = [&]() {
            denseUnitaryMatrix[3]->unitaryTransformAndReturnMainDiagonal(
                sparseDiagonalizableMatrices[3]);
        };
        std::cout << "Armadillo double precision:" << std::endl;
        PerformanceTest(armaDoubleUnitaryTransformation, cycles);

        std::cout << "Armadillo single precision:" << std::endl;
        PerformanceTest(armaSingleUnitaryTransformation, cycles);

        std::cout << "Eigen double precision:" << std::endl;
        PerformanceTest(eigenDoubleUnitaryTransformation, cycles);

        std::cout << "Eigen single precision:" << std::endl;
        PerformanceTest(eigenSingleUnitaryTransformation, cycles);
        std::cout << std::endl;
    }
}
