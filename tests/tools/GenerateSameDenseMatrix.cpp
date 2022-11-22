#include "GenerateSameDenseMatrix.h"

std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseMatrix>>
generateSymmetricMatrices(
    size_t size,
    const std::vector<std::shared_ptr<quantum::linear_algebra::AbstractDenseFactory>>& factories,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng) {
    std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseMatrix>> answer;
    answer.reserve(factories.size());

    // create empty matrices:
    for (const auto& factory : factories) {
        answer.emplace_back(factory->createMatrix());
    }

    // resize it:
    for (auto& matrix : answer) {
        matrix->resize(size, size);
    }

    // fill it with identical values:
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            double value = dist(rng);
            for (auto& matrix : answer) {
                matrix->assign_to_position(value, i, j);
                matrix->assign_to_position(value, j, i);
            }
        }
    }

    return answer;
}

void multiplyColumnByMinusOne(
    std::unique_ptr<quantum::linear_algebra::AbstractDenseMatrix>& rhs,
    size_t column) {
    for (size_t i = 0; i < rhs->size_rows(); ++i) {
        rhs->assign_to_position(-rhs->at(i, column), i, column);
    }
}

void makeUnitaryMatrixSame(
    const std::unique_ptr<quantum::linear_algebra::AbstractDenseMatrix>& lhs,
    std::unique_ptr<quantum::linear_algebra::AbstractDenseMatrix>& rhs) {
    for (size_t column = 0; column < lhs->size_cols(); ++column) {
        if (std::abs(lhs->at(0, column) - (-rhs->at(0, column))) < 1e-6) {
            multiplyColumnByMinusOne(rhs, column);
        }
    }
}

std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseMatrix>> generateUnitaryMatrix(
    size_t size,
    const std::vector<std::shared_ptr<quantum::linear_algebra::AbstractDenseFactory>>& factories,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng) {
    // construct symmetrical matrix:
    auto symmetricMatrices = generateSymmetricMatrices(size, factories, dist, rng);

    auto unitaryMatrices =
        std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseMatrix>>(
            symmetricMatrices.size());

    // construct unitary matrix as eigenvectors matrix:
    for (size_t i = 0; i < symmetricMatrices.size(); ++i) {
        unitaryMatrices[i] = symmetricMatrices[i]->diagonalizeValuesVectors().eigenvectors;
    }

    for (size_t i = 1; i < unitaryMatrices.size(); ++i) {
        makeUnitaryMatrixSame(unitaryMatrices[0], unitaryMatrices[i]);
    }

    return unitaryMatrices;
}
