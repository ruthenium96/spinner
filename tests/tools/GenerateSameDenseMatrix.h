#ifndef SPINNER_GENERATESAMEDENSEMATRIX_H
#define SPINNER_GENERATESAMEDENSEMATRIX_H

#include <random>
#include <vector>

#include "src/entities/data_structures/AbstractDenseFactory.h"
#include "tests/tools/AllDenseFactories.h"

std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseMatrix>>
generateSymmetricMatrices(
    size_t size,
    const std::vector<std::shared_ptr<quantum::linear_algebra::AbstractDenseFactory>>& factories,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng);

std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseMatrix>> generateUnitaryMatrix(
    size_t size,
    const std::vector<std::shared_ptr<quantum::linear_algebra::AbstractDenseFactory>>& factories,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng);

void makeUnitaryMatrixSame(
    const std::unique_ptr<quantum::linear_algebra::AbstractDenseMatrix>& lhs,
    std::unique_ptr<quantum::linear_algebra::AbstractDenseMatrix>& rhs);

#endif  //SPINNER_GENERATESAMEDENSEMATRIX_H