#include "GenerateSameDenseMatrix.h"

std::vector<std::unique_ptr<quantum::linear_algebra::AbstractSymmetricMatrix>>
generateSymmetricMatrices(
    size_t size,
    const std::vector<std::shared_ptr<quantum::linear_algebra::AbstractSymmetricMatrixFactory>>&
        factories,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng) {
    std::vector<std::unique_ptr<quantum::linear_algebra::AbstractSymmetricMatrix>> answer;
    answer.reserve(factories.size());

    // create zero matrices:
    for (const auto& factory : factories) {
        answer.emplace_back(factory->createSymmetricMatrix(size));
    }

    // fill it with identical values:
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            double value = dist(rng);
            for (auto& matrix : answer) {
                matrix->add_to_position(value, i, j);
                matrix->add_to_position(value, j, i);
            }
        }
    }

    return answer;
}

std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
generateUnitaryMatrix(
    size_t size,
    const std::vector<std::shared_ptr<quantum::linear_algebra::AbstractSymmetricMatrixFactory>>&
        factories,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng) {
    // construct symmetrical matrix:
    auto symmetricMatrices = generateSymmetricMatrices(size, factories, dist, rng);

    auto unitaryMatrices =
        std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>(
            symmetricMatrices.size());

    // construct unitary matrix as eigenvectors matrix:
    for (size_t i = 0; i < symmetricMatrices.size(); ++i) {
        unitaryMatrices[i] = symmetricMatrices[i]->diagonalizeValuesVectors().eigenvectors;
    }

    return unitaryMatrices;
}
