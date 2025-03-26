#ifndef SPINNER_GENERATESAMEMATRIX_H
#define SPINNER_GENERATESAMEMATRIX_H

#include <random>
#include <vector>

#include "src/entities/data_structures/AbstractFactories.h"

std::unique_ptr<quantum::linear_algebra::AbstractDiagonalizableMatrix>
generateSparseDiagonalizableMatrix(
    size_t size,
    std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory> factory,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng);

std::unique_ptr<quantum::linear_algebra::AbstractDiagonalizableMatrix>
generateDenseDiagonalizableMatrix(
    size_t size,
    std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory> factory,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng);

std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix> generateDenseUnitaryMatrix(
    size_t size,
    std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory> factory,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng);

std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>
generateOrthDenseVector(
    size_t size,
    std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory> factory);

std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDiagonalizableMatrix>>
generateDenseDiagonalizableMatrices(
    size_t size,
    const std::vector<
        std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>>&
        factories,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng);

std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDiagonalizableMatrix>>
generateSparseDiagonalizableMatrices(
    size_t size,
    const std::vector<
        std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>>&
        factories,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng);

std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
generateDenseUnitaryMatrices(
    size_t size,
    const std::vector<
        std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>>&
        factories,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng);

std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>
generateOrthDenseVectors(
    size_t size,
    const std::vector<
        std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>>&
        factories);

void makeUnitaryMatrixSame(
    const std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>& lhs,
    std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>& rhs);

#endif  //SPINNER_GENERATESAMEMATRIX_H
