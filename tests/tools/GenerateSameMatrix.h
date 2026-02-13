#ifndef SPINNER_GENERATESAMEMATRIX_H
#define SPINNER_GENERATESAMEMATRIX_H

#include <random>
#include <vector>

#include "src/linalg_structures/AbstractFactories.h"

std::unique_ptr<spinner::linalg_structures::AbstractDiagonalizableMatrix>
generateSparseDiagonalizableMatrix(
    size_t size,
    std::shared_ptr<spinner::linalg_structures::AbstractDenseTransformAndDiagonalizeFactory> factory,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng);

std::unique_ptr<spinner::linalg_structures::AbstractDiagonalizableMatrix>
generateDenseDiagonalizableMatrix(
    size_t size,
    std::shared_ptr<spinner::linalg_structures::AbstractDenseTransformAndDiagonalizeFactory> factory,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng);

std::unique_ptr<spinner::linalg_structures::AbstractDenseSemiunitaryMatrix> generateDenseUnitaryMatrix(
    size_t size,
    std::shared_ptr<spinner::linalg_structures::AbstractDenseTransformAndDiagonalizeFactory> factory,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng);

std::unique_ptr<spinner::linalg_structures::AbstractDenseVector>
generateOrthDenseVector(
    size_t size,
    std::shared_ptr<spinner::linalg_structures::AbstractDenseTransformAndDiagonalizeFactory> factory);

std::vector<std::unique_ptr<spinner::linalg_structures::AbstractDiagonalizableMatrix>>
generateDenseDiagonalizableMatrices(
    size_t size,
    const std::vector<
        std::shared_ptr<spinner::linalg_structures::AbstractDenseTransformAndDiagonalizeFactory>>&
        factories,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng);

std::vector<std::unique_ptr<spinner::linalg_structures::AbstractDiagonalizableMatrix>>
generateSparseDiagonalizableMatrices(
    size_t size,
    const std::vector<
        std::shared_ptr<spinner::linalg_structures::AbstractDenseTransformAndDiagonalizeFactory>>&
        factories,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng);

std::vector<std::unique_ptr<spinner::linalg_structures::AbstractDenseSemiunitaryMatrix>>
generateDenseUnitaryMatrices(
    size_t size,
    const std::vector<
        std::shared_ptr<spinner::linalg_structures::AbstractDenseTransformAndDiagonalizeFactory>>&
        factories,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng);

std::vector<std::unique_ptr<spinner::linalg_structures::AbstractDenseVector>>
generateOrthDenseVectors(
    size_t size,
    const std::vector<
        std::shared_ptr<spinner::linalg_structures::AbstractDenseTransformAndDiagonalizeFactory>>&
        factories);

void makeUnitaryMatrixSame(
    const std::unique_ptr<spinner::linalg_structures::AbstractDenseSemiunitaryMatrix>& lhs,
    std::unique_ptr<spinner::linalg_structures::AbstractDenseSemiunitaryMatrix>& rhs);

#endif  //SPINNER_GENERATESAMEMATRIX_H
