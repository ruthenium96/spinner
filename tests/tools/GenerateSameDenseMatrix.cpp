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
        answer.emplace_back(factory->createDenseSymmetricMatrix(size));
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

std::vector<std::unique_ptr<quantum::linear_algebra::AbstractSymmetricMatrix>>
generateSparseSymmetricMatrices(
    size_t size,
    const std::vector<std::shared_ptr<quantum::linear_algebra::AbstractSymmetricMatrixFactory>>&
        factories,
    std::uniform_real_distribution<double> dist,
    std::mt19937 rng) {
    std::vector<std::unique_ptr<quantum::linear_algebra::AbstractSymmetricMatrix>> answer;
    answer.reserve(factories.size());

    // create zero matrices:
    for (const auto& factory : factories) {
        answer.emplace_back(factory->createSparseSymmetricMatrix(size));
    }

    size_t numberOfElementsPerRow = log(size) * (log(size) - 1);

    std::uniform_int_distribution<size_t> distOfColumn(0, size - 1);

    // fill it with identical values:
    for (size_t i = 0; i < size; ++i) {
        size_t addedValues = 0;
        while (addedValues < numberOfElementsPerRow) {
            size_t j = distOfColumn(rng);
            if (answer[0]->at(i, j) == 0) {
                double value = dist(rng);
                for (auto& matrix : answer) {
                    matrix->add_to_position(value, i, j);
                    matrix->add_to_position(value, j, i);
                }
                addedValues++;
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
