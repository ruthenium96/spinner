#include <random>

#include "gtest/gtest.h"
#include "tests/tools/AllSymmetricMatrixFactories.h"
#include "tests/tools/GenerateSameMatrix.h"

inline int sign_f(double a, double b)
{
    return -2 * (std::signbit(a) ^ std::signbit(b)) + 1;
}

TEST(linearAlgebraFactories, throw_combination_of_different_objects) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    auto denseUnitaryMatrices = generateDenseUnitaryMatrices(
        10,
        constructAllDenseTransformAndDiagonalizeFactories(),
        dist,
        rng);
    auto denseDiagonalizableMatrices = generateDenseDiagonalizableMatrices(
        10,
        constructAllDenseTransformAndDiagonalizeFactories(),
        dist,
        rng);

    for (int i = 0; i < denseDiagonalizableMatrices.size(); ++i) {
        auto& denseDiagonalizableMatrix_i = denseDiagonalizableMatrices[i];
        for (int j = 0; j < denseUnitaryMatrices.size(); ++j) {
            auto& denseUnitatyMatrix_j = denseUnitaryMatrices[j];
            if (i != j) {
                EXPECT_ANY_THROW(denseUnitatyMatrix_j->
                    getUnitaryTransformer()->calculateUnitaryTransformationOfMatrix(denseDiagonalizableMatrix_i));
            } else {
                EXPECT_NO_THROW(denseUnitatyMatrix_j->
                    getUnitaryTransformer()->calculateUnitaryTransformationOfMatrix(denseDiagonalizableMatrix_i));
            }
        }
    }

    std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>> denseVectors;
    for (int i = 0; i < denseUnitaryMatrices.size(); ++i) {
        denseVectors.emplace_back(denseDiagonalizableMatrices[i]->diagonalizeValues());
    }

    for (int i = 0; i < denseVectors.size(); ++i) {
        auto& vector_i = denseVectors[i];
        for (int j = 0; j < denseVectors.size(); ++j) {
            auto& vector_j = denseVectors[j];
            if (i != j) {
                EXPECT_ANY_THROW(vector_i->concatenate_with(vector_j));
                EXPECT_ANY_THROW(vector_i->dot(vector_j));
                EXPECT_ANY_THROW(vector_i->element_wise_multiplication(vector_j));
            } else {
                EXPECT_NO_THROW(vector_i->concatenate_with(vector_j));
                EXPECT_NO_THROW(vector_i->dot(vector_j));
                EXPECT_NO_THROW(vector_i->element_wise_multiplication(vector_j));
            }
        }
    }
}

// Check if different linear algebra packages make the same eigen decomposition
TEST(linearAlgebraFactories, eigendecomposition) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 2; size < 100; ++size) {
        // construct identical dense diagonalizable matrix:
        auto matrices = generateDenseDiagonalizableMatrices(
            size,
            constructAllDenseTransformAndDiagonalizeFactories(),
            dist,
            rng);
        // only-values-eigendecomposition:
        {
            // decomposition
            std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>> denseVectors;
            for (const auto& matrix : matrices) {
                denseVectors.emplace_back(matrix->diagonalizeValues());
            }
            // check equality:
            for (size_t i = 0; i < denseVectors.size(); ++i) {
                const auto& denseVector_i = denseVectors[i];
                for (size_t j = 0; j < denseVectors.size(); ++j) {
                    if (i == j) {
                        continue;
                    }
                    const auto& denseVector_j = denseVectors[j];
                    for (size_t k = 0; k < size; ++k) {
                        double epsilon = std::abs(denseVector_i->at(k) * 5e-3);
                        EXPECT_NEAR(denseVector_i->at(k), denseVector_j->at(k), epsilon);
                    }
                }
            }
        }
        // eigendecomposition-with-eigenvectors
        {
            // decomposition
            std::vector<quantum::linear_algebra::EigenCouple> denseEigenCouples;
            for (const auto& matrix : matrices) {
                denseEigenCouples.emplace_back(matrix->diagonalizeValuesVectors());
            }
            // check equality:
            for (size_t i = 0; i < denseEigenCouples.size(); ++i) {
                auto& eigencouple_i = denseEigenCouples[i];
                for (size_t j = 0; j < denseEigenCouples.size(); ++j) {
                    if (i == j) {
                        continue;
                    }
                    auto& eigencouple_j = denseEigenCouples[j];
                    for (size_t k = 0; k < size; ++k) {
                        double epsilon_values = std::abs(eigencouple_i.eigenvalues->at(k) * 5e-3);
                        EXPECT_NEAR(
                            eigencouple_i.eigenvalues->at(k),
                            eigencouple_j.eigenvalues->at(k),
                            epsilon_values);
                        int sign;
                        for (size_t l = 0; l < size; ++l) {
                            if (l == 0) {
                                sign = sign_f(eigencouple_i.eigenvectors->at(k, l),
                                              eigencouple_j.eigenvectors->at(k, l));
                            }
                            // TODO: epsilon
                            EXPECT_NEAR(
                                eigencouple_i.eigenvectors->at(k, l),
                                sign * eigencouple_j.eigenvectors->at(k, l),
                                5e-4);
                        }
                    }
                }
            }
        }
    }
}

TEST(linearAlgebraFactories, krylov_eigendecomposition) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 8; size <= 16; size*=2) {
        // construct identical dense diagonalizable matrix:
        auto matrices = generateSparseDiagonalizableMatrices(
            size,
            constructAllDenseTransformAndDiagonalizeFactories(),
            dist,
            rng);
        auto orth_vectors = generateOrthDenseVectors(
            size, 
            constructAllDenseTransformAndDiagonalizeFactories());
        // only-values-eigendecomposition:
        {
            // decomposition
            std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>> denseVectorsEnergy;
            std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>> denseVectorsWeights;
            for (int i = 0; i < matrices.size(); ++i) {
                const auto& matrix = matrices[i];
                const auto& orth_vector = orth_vectors[i];
                auto pair = matrix->krylovDiagonalizeValues(orth_vector, size);
                denseVectorsEnergy.emplace_back(std::move(pair.eigenvalues));
                denseVectorsWeights.emplace_back(std::move(pair.ftlm_weights_of_states));
            }
            // check equality:
            for (size_t i = 0; i < denseVectorsEnergy.size(); ++i) {
                const auto& denseVectorEnergy_i = denseVectorsEnergy[i];
                const auto& denseVectorsWeights_i = denseVectorsWeights[i];
                for (size_t j = 0; j < denseVectorsEnergy.size(); ++j) {
                    if (i == j) {
                        continue;
                    }
                    const auto& denseVectorEnergy_j = denseVectorsEnergy[j];
                    const auto& denseVectorsWeights_j = denseVectorsWeights[j];
                    for (size_t k = 0; k < size; ++k) {
                        double epsilonEnergy = std::abs(denseVectorEnergy_i->at(k) * 5e-3);
                        EXPECT_NEAR(denseVectorEnergy_i->at(k), denseVectorEnergy_j->at(k), epsilonEnergy);
                        EXPECT_NEAR(denseVectorsWeights_i->at(k), denseVectorsWeights_j->at(k), 1e-3);
                    }
                }
            }
        }
        // eigendecomposition-with-eigenvectors
        {
            // decomposition
            std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>> denseVectorsEnergy;
            std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>> denseVectorsWeights;
            std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>> denseUnitaryMatrices;
            for (int i = 0; i < matrices.size(); ++i) {
                const auto& matrix = matrices[i];
                const auto& orth_vector = orth_vectors[i];
                auto triple = matrix->krylovDiagonalizeValuesVectors(orth_vector, size);
                denseVectorsEnergy.emplace_back(std::move(triple.eigenvalues));
                denseVectorsWeights.emplace_back(std::move(triple.ftlm_weights_of_states));
                denseUnitaryMatrices.emplace_back(std::move(triple.eigenvectors));
            }
            // check equality:
            for (size_t i = 0; i < denseVectorsEnergy.size(); ++i) {
                const auto& denseVectorEnergy_i = denseVectorsEnergy[i];
                const auto& denseVectorsWeights_i = denseVectorsWeights[i];
                for (size_t j = 0; j < denseVectorsEnergy.size(); ++j) {
                    if (i == j) {
                        continue;
                    }
                    const auto& denseVectorEnergy_j = denseVectorsEnergy[j];
                    const auto& denseVectorsWeights_j = denseVectorsWeights[j];
                    for (size_t k = 0; k < size; ++k) {
                        double epsilonEnergy = std::abs(denseVectorEnergy_i->at(k) * 5e-3);
                        EXPECT_NEAR(denseVectorEnergy_i->at(k), denseVectorEnergy_j->at(k), epsilonEnergy);
                        EXPECT_NEAR(denseVectorsWeights_i->at(k), denseVectorsWeights_j->at(k), 1e-3);
                        int sign;
                        for (size_t l = 0; l < size; ++l) {
                            if (l == 0) {
                                sign = sign_f(denseUnitaryMatrices[i]->at(k, l),
                                denseUnitaryMatrices[j]->at(k, l));
                            }
                            // TODO: epsilon
                            EXPECT_NEAR(
                                denseUnitaryMatrices[i]->at(k, l),
                                sign * denseUnitaryMatrices[j]->at(k, l),
                                5e-3);
                        }
                    }
                }
            }
        }
    }
}

// Check if different linear algebra packages make the same unitary transformation
TEST(linearAlgebraFactories, unitary_transformation_and_return_main_diagonal) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 2; size < 100; ++size) {
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

        // unitary transformation and returning main diagonal:
        std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>
            denseTransformedMainDiagonal;
        for (size_t i = 0; i < denseDiagonalizableMatrices.size(); ++i) {
            const auto& denseDiagonalizableMatrix = denseDiagonalizableMatrices[i];
            const auto& denseUnitaryMatrix = denseUnitaryMatrices[i];

            denseTransformedMainDiagonal.emplace_back(
                denseUnitaryMatrix->getUnitaryTransformer()->calculateUnitaryTransformationOfMatrix(denseDiagonalizableMatrix));
        }

        // check equality:
        for (size_t i = 0; i < denseTransformedMainDiagonal.size(); ++i) {
            const auto& transformedMainDiagonal_i = denseTransformedMainDiagonal[i];
            for (size_t j = 0; j < denseTransformedMainDiagonal.size(); ++j) {
                if (i == j) {
                    continue;
                }
                const auto& transformedMainDiagonal_j = denseTransformedMainDiagonal[j];
                for (size_t k = 0; k < size; ++k) {
                    double epsilon = std::max(
                        std::abs(transformedMainDiagonal_i->at(k) * 1e-4),
                        5e-3);
                    EXPECT_NEAR(
                        transformedMainDiagonal_i->at(k),
                        transformedMainDiagonal_j->at(k),
                        epsilon);
                }
            }
        }
    }
}