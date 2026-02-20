#ifndef SPINNER_ABSTRACTFACTORIES_H
#define SPINNER_ABSTRACTFACTORIES_H

#include <vector>
#include "AbstractDenseSemiunitaryMatrix.h"
#include "AbstractDenseVector.h"
#include "AbstractDiagonalizableMatrix.h"
#include "AbstractSparseSemiunitaryMatrix.h"
#include "AbstractSymmetricMatrix.h"

namespace spinner::linalg_structures {

enum Precision { SINGLE, DOUBLE };

class AbstractDenseTransformAndDiagonalizeFactory {
  public:
    static std::shared_ptr<AbstractDenseTransformAndDiagonalizeFactory> defaultFactory();
    virtual std::unique_ptr<AbstractDiagonalizableMatrix>
    createDenseDiagonalizableMatrix(uint32_t size) = 0;
    virtual std::unique_ptr<AbstractDiagonalizableMatrix>
    createSparseDiagonalizableMatrix(uint32_t size) = 0;
    virtual std::unique_ptr<AbstractDenseSemiunitaryMatrix>
    createDenseSemiunitaryMatrix(uint32_t cols, uint32_t rows) = 0;
    /**
     * @brief Creates a given number of random unit vectors.
     * For each vector 1) independently generates size_of_vector elements 
     * according to some distribution, 2) Euclid-normalize these elements.
     * @param size_of_vector  Dimension of each vector.
     * @param number_of_vectors How many vectors to create.
     * @return std::vector of unique_ptrs to AbstractDenseVector, each of unit norm.
     */
    virtual std::vector<std::unique_ptr<AbstractDenseVector>> createRandomUnitVectors(uint32_t size_of_vector, uint32_t number_of_vectors) = 0;
    /**
     * @brief Creates an empty dense vector.
     * @return std::unique_ptr<AbstractDenseVector>
     */
    virtual std::unique_ptr<AbstractDenseVector> createVector() = 0;

    /// Virtual destructor.
    ~AbstractDenseTransformAndDiagonalizeFactory() = default;
  private:
    Precision precision_ = Precision::DOUBLE;
  public:
    void setPrecision(Precision precision) {
        precision_ = precision;
    }
    Precision getPrecision() const {
        return precision_;
    }
};

/**
 * @class AbstractSparseTransformFactory
 * @brief Factory for sparse matrix structures tailored to incremental assembly.
 *
 * This factory produces sparse semi‑unitary and sparse symmetric matrices that
 * are built incrementally via `add_to_position()` calls. For sparse semi‑unitary
 * matrices, the key operation `move_vector_from()` efficiently transfers whole
 * rows between matrices, which is essential for the performance.
 *
 * @note These structures deliberately avoid using standard sparse storage
 *       formats (like CSC or CSR) because implementing `add_to_position` and
 *       `move_vector_from` efficiently in those formats would be impractical.
 *       Instead, the concrete implementations use custom data layouts optimized
 *       for incremental construction and row moves (e.g. hashmaps). Consequently,
 *       they are provided by a separate factory, distinct from the dense linear 
 *       algebra backend factory.
 */
class AbstractSparseTransformFactory {
  public:
    /**
     * @brief Returns a default (concrete) factory instance.
     *
     * The exact type depends on the CMake configuration.
     *
     * @return std::shared_ptr<AbstractSparseTransformFactory>
     */
    static std::shared_ptr<AbstractSparseTransformFactory> defaultSparseFactory();
    /**
     * @brief Creates a sparse semi‑unitary matrix.
     * @param cols Number of columns (vectors).
     * @param rows Number of rows.
     * @return std::unique_ptr<AbstractSparseSemiunitaryMatrix>
     */
    virtual std::unique_ptr<AbstractSparseSemiunitaryMatrix>
    createSparseSemiunitaryMatrix(uint32_t cols, uint32_t rows) = 0;
    /**
     * @brief Creates a sparse symmetric matrix.
     * @param size Matrix dimension.
     * @return std::unique_ptr<AbstractSymmetricMatrix>
     */
    virtual std::unique_ptr<AbstractSymmetricMatrix> createSparseSymmetricMatrix(uint32_t size) = 0;

    /// Virtual destructor.
    ~AbstractSparseTransformFactory() = default;
};

} // namespace spinner::linalg_structures

#endif  //SPINNER_ABSTRACTFACTORIES_H
