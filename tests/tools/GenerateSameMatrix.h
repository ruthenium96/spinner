#ifndef SPINNER_GENERATESAMEMATRIX_H
#define SPINNER_GENERATESAMEMATRIX_H

#include <random>
#include <vector>

#include "src/entities/data_structures/AbstractFactories.h"

std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDiagonalizableMatrix>>
generateSymmetricMatrices(
    size_t size,
    const std::vector<
        std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>>&
        factories,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng);

std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDiagonalizableMatrix>>
generateSparseSymmetricMatrices(
    size_t size,
    const std::vector<
        std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>>&
        factories,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng);

std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
generateUnitaryMatrix(
    size_t size,
    const std::vector<
        std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>>&
        factories,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng);

void makeUnitaryMatrixSame(
    const std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>& lhs,
    std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>& rhs);

#endif  //SPINNER_GENERATESAMEMATRIX_H
