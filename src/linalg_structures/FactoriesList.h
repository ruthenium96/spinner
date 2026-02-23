#ifndef SPINNER_FACTORIESLIST_H
#define SPINNER_FACTORIESLIST_H

#include <cstdint>
#include <vector>
#include "AbstractFactories.h"

// TODO: refactor the whole project and pass FactoriesList only to some controlling class.
//  in case of space/optimization it can be done by means of replacing move_vector_from to
//  move_vector_to function. In other cases it can be done through refactoring of constructors.

namespace spinner::linalg_structures {
class FactoriesList {
  public:
    explicit FactoriesList(
        std::shared_ptr<AbstractDenseTransformAndDiagonalizeFactory> symmetricMatrixFactory =
            AbstractDenseTransformAndDiagonalizeFactory::defaultFactory(),
        std::shared_ptr<AbstractSparseTransformFactory> sparseMatrix =
            AbstractSparseTransformFactory::defaultSparseFactory());

    std::unique_ptr<AbstractDiagonalizableMatrix>
    createDenseDiagonalizableMatrix(uint32_t size) const;
    std::unique_ptr<AbstractDiagonalizableMatrix>
    createSparseDiagonalizableMatrix(uint32_t size) const;
    /**
     * @brief Creates a given number of random unit vectors (using AbstractDenseTransformAndDiagonalizeFactory).
     * For each vector 1) independently generates size_of_vector elements 
     * according to some distribution, 2) Euclid-normalize these elements.
     * @param size_of_vector  Dimension of each vector.
     * @param number_of_vectors How many vectors to create.
     * @return std::vector of unique_ptrs to AbstractDenseVector, each of unit norm.
     */
    std::vector<std::unique_ptr<AbstractDenseVector>> createRandomUnitVectors(uint32_t size_of_vector, uint32_t number_of_vectors) const;
    /**
     * @brief Creates an empty dense vector (using AbstractDenseTransformAndDiagonalizeFactory).
     * @return std::unique_ptr<AbstractDenseVector>
     */
    std::unique_ptr<AbstractDenseVector> createVector() const;
    /**
     * @brief Creates a sparse semi‑unitary matrix (using AbstractSparseTransformFactory).
     * @param cols Number of columns (vectors).
     * @param rows Number of rows.
     * @return std::unique_ptr<AbstractSparseSemiunitaryMatrix>
     */
    std::unique_ptr<AbstractSparseSemiunitaryMatrix>
    createSparseSemiunitaryMatrix(uint32_t cols, uint32_t rows) const;
    /**
     * @brief Creates a sparse symmetric matrix (using AbstractSparseTransformFactory).
     * @param size Matrix dimension.
     * @return std::unique_ptr<AbstractSymmetricMatrix>
     */
    std::unique_ptr<AbstractSymmetricMatrix> createSparseSymmetricMatrix(uint32_t size) const;

  private:
    std::shared_ptr<AbstractDenseTransformAndDiagonalizeFactory> denseFactory_;
    std::shared_ptr<AbstractSparseTransformFactory> sparseFactory_;
};
} // namespace spinner::linalg_structures
#endif  //SPINNER_FACTORIESLIST_H
